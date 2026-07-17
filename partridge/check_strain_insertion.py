"""
For each candidate insertion breakpoint in a collect_evidence.py evidence
file, check whether the insertion already exists in a second (query) genome
rather than being novel/somatic — e.g. checking a C57BL/6-called breakpoint
against a CD-1 genome, to rule out a strain-polymorphic insertion that's
simply absent from the C57BL/6 reference.

Method, per breakpoint:
  1. Extract a window around the breakpoint from the reference fasta.
  2. Extract the corresponding window from the query fasta — either the same
     coordinates directly (--same-scaffold) or, by default, by locating the
     orthologous locus with minimap2.
  3. Align the reference window to the query window directly (a coarse,
     first-pass presence/absence signal on its own).
  4. Align the breakpoint's RIGHT_CLIP/LEFT_CLIP sequences (pure insertion
     sequence, per collect_evidence.py) against both windows. A real hit in
     the query window (and none in the reference window, as a sanity check)
     is evidence the insertion pre-exists in the query genome.

This reports alignment metrics for manual review — it does not make the
present/absent call for you (a "note" column applies a simple, overridable
threshold as a starting point only).

Multiple evidence files (e.g. one per sample) can be checked in a single run
via --evidence-dir — output rows gain a "sample" column (each evidence
file's basename with its extension stripped), and the companion windows.fa
gets sample-prefixed headers, so results from every sample land in one
combined tsv/fasta pair.

Dependencies beyond this repo's usual env.txt: biopython>=1.80 (pairwise
local alignment), and the minimap2 binary on PATH (only needed unless
--same-scaffold is given).
"""

import argparse
import glob
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field

import pandas as pd
import pysam
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

SECTION_NAMES = {"RIGHT_CLIP", "LEFT_CLIP", "LEFT_MATE", "RIGHT_MATE"}
HEADER_RE = re.compile(
    r"^@(?P<seqname>[^:]+):(?P<pos1>\d+)-(?P<pos2>\d+) \((?P<strand>[+-])\)$"
)

logger = logging.getLogger("check_strain_insertion")


@dataclass
class Breakpoint:
    seqname: str
    pos1: int
    pos2: int
    strand: str
    right_clips: list = field(default_factory=list)
    left_clips: list = field(default_factory=list)

    @property
    def name(self):
        return f"{self.seqname}:{self.pos1}"

    @property
    def longest_right_clip(self):
        return max(self.right_clips, key=len, default="")

    @property
    def longest_left_clip(self):
        return max(self.left_clips, key=len, default="")


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_evidence_fa(path):
    """Yield one Breakpoint per entry in a collect_evidence.py output file."""
    bp = None
    section = None
    with open_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            header = HEADER_RE.match(line)
            if header:
                if bp is not None:
                    yield bp
                bp = Breakpoint(
                    seqname=header["seqname"],
                    pos1=int(header["pos1"]),
                    pos2=int(header["pos2"]),
                    strand=header["strand"],
                )
                section = None
                continue
            if line[1:] in SECTION_NAMES if line.startswith(">") else False:
                section = line[1:]
                continue
            if line.startswith(">"):
                # LEFT_CLIP data line: ">" + right-alignment padding + sequence.
                seq = line[1:].lstrip(" ")
                if bp is not None and section == "LEFT_CLIP" and seq:
                    bp.left_clips.append(seq)
                continue
            # plain sequence line (RIGHT_CLIP, or a mate section we don't use here)
            if bp is not None and section == "RIGHT_CLIP" and line != "None":
                bp.right_clips.append(line)
    if bp is not None:
        yield bp


def sample_name(path):
    """Derive a sample name from an evidence file path by stripping its
    extension(s), e.g. '10.3C.fa.gz' -> '10.3C', 'sample1.fasta' -> 'sample1'."""
    base = os.path.basename(path)
    for suffix in (".fa.gz", ".fasta.gz", ".fa", ".fasta"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def find_evidence_files(directory):
    files = sorted(
        p
        for p in glob.glob(os.path.join(directory, "*"))
        if p.endswith(".fa") or p.endswith(".fa.gz")
    )
    if not files:
        logger.error("No *.fa / *.fa.gz evidence files found in %s", directory)
        sys.exit(1)
    return files


_TEMP_FILES = []


def is_bgzf(path):
    """Check the BGZF magic (a gzip 'extra field' subfield, BC) directly,
    rather than probing with pysam — pysam/htslib logs a scary-looking but
    harmless error straight to the C stderr stream if given a plain-gzip
    file, before Python ever sees the exception."""
    with open(path, "rb") as f:
        header = f.read(18)
    if len(header) < 18 or header[0:2] != b"\x1f\x8b":
        return False
    if not (header[3] & 0x04):  # FEXTRA flag
        return False
    xlen = int.from_bytes(header[10:12], "little")
    return xlen >= 6 and header[12:14] == b"BC"


def prepare_fasta(path):
    """Open `path` as a pysam.FastaFile, returning (FastaFile, resolved_path).

    pysam/htslib's faidx can only random-access a .gz fasta if it's
    bgzip-compressed, not plain gzip. If `path` is plain-gzip, transparently
    decompress it to a temporary plain fasta first and index/use that
    instead (resolved_path then points at the temp file, not `path`).
    Uncompressed and bgzip fasta are used directly, no copy made.
    """
    resolved = path
    if path.endswith(".gz") and not is_bgzf(path):
        logger.warning(
            "'%s' is plain-gzip, not bgzip, so it can't be indexed directly "
            "— decompressing to a temporary fasta for random access (this "
            "may take a while for a full genome; bgzip-compressing it ahead "
            "of time avoids this next run)",
            path,
        )
        tmp = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)
        tmp.close()
        with gzip.open(path, "rb") as src, open(tmp.name, "wb") as dst:
            shutil.copyfileobj(src, dst)
        resolved = tmp.name
        _TEMP_FILES.append(resolved)
        logger.debug("'%s' decompressed to temporary file '%s'", path, resolved)
    if not os.path.exists(resolved + ".fai"):
        logger.debug("Building fasta index for '%s'", resolved)
        pysam.faidx(resolved)
    return pysam.FastaFile(resolved), resolved


def resolve_seqname(fasta, seqname):
    """Look up `seqname` in `fasta`, tolerating a UCSC-style 'chr' prefix
    mismatch between the evidence file's coordinates and this fasta's own
    naming convention (e.g. evidence built against a UCSC-named genome, but
    checked against an Ensembl-named reference/query fasta, or vice versa).
    """
    references = fasta.references
    if seqname in references:
        return seqname
    alt = seqname[3:] if seqname.startswith("chr") else f"chr{seqname}"
    if alt in references:
        logger.debug("'%s' not in fasta, using '%s' instead", seqname, alt)
        return alt
    raise KeyError(
        f"'{seqname}' not found in fasta (also tried '{alt}'); this fasta's "
        f"sequence names look like: {list(references[:3])}"
    )


def clamp_window(fasta, seqname, center, window):
    seqlen = fasta.get_reference_length(seqname)
    start = max(0, center - window)
    end = min(seqlen, center + window)
    return start, end


def extract_window(fasta, seqname, center, window):
    seqname = resolve_seqname(fasta, seqname)
    start, end = clamp_window(fasta, seqname, center, window)
    return fasta.fetch(seqname, start, end), start, end


def build_minimap2_index(query_fasta_path, minimap2_bin, preset, threads):
    """Build a minimap2 index (.mmi) for the query genome once, so every
    breakpoint's orthology lookup can reuse it instead of each one making
    minimap2 re-index the whole genome from scratch (the dominant cost for a
    full genome — see locate_orthologous_windows_batch)."""
    index_path = tempfile.NamedTemporaryFile(suffix=".mmi", delete=False).name
    _TEMP_FILES.append(index_path)
    logger.info("Building minimap2 index (-x %s) for %s ...", preset, query_fasta_path)
    try:
        subprocess.run(
            [minimap2_bin, "-x", preset, "-t", str(threads), "-d", index_path, query_fasta_path],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        logger.error(
            "'%s' not found on PATH. Install minimap2, or pass --same-scaffold "
            "if the reference and query fastas already share coordinates.",
            minimap2_bin,
        )
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        logger.error("minimap2 indexing failed: %s", e.stderr)
        sys.exit(1)
    logger.info("minimap2 index built: %s", index_path)
    return index_path


def locate_orthologous_windows_batch(
    records, minimap2_index, query_fasta, minimap2_bin, preset, threads
):
    """Find the query-genome window orthologous to each record's reference
    window, in a single minimap2 call against the pre-built index — instead
    of one subprocess (and, previously, one from-scratch index build) per
    breakpoint.

    `records` is an iterable of objects with a unique `.key` and a
    `.ref_window_seq`. Returns {key: (seq, scaffold, start, end, strand,
    hit_span_diff)} for keys with a hit; keys with no hit are simply absent
    from the returned dict. hit_span_diff = (target span) - (query span)
    for that hit: a value close to the expected insertion size is itself a
    signal the insertion is present in the query genome (a real minimap2
    gap spanning the missing sequence), independent of the clip alignments.
    """
    records = list(records)
    if not records:
        return {}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        for rec in records:
            tmp.write(f">{rec.key}\n{rec.ref_window_seq}\n")
        tmp_path = tmp.name

    logger.info("Running one batched minimap2 query for %d breakpoint(s)...", len(records))
    try:
        try:
            result = subprocess.run(
                [
                    minimap2_bin,
                    "-x",
                    preset,
                    "-t",
                    str(threads),
                    "--secondary=no",
                    minimap2_index,
                    tmp_path,
                ],
                capture_output=True,
                text=True,
                check=True,
            )
        except FileNotFoundError:
            logger.error(
                "'%s' not found on PATH. Install minimap2, or pass --same-scaffold "
                "if the reference and query fastas already share coordinates.",
                minimap2_bin,
            )
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            logger.warning("minimap2 batch mapping failed: %s", e.stderr)
            return {}
    finally:
        os.unlink(tmp_path)

    best_by_key = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        key, nmatch = fields[0], int(fields[9])
        if key not in best_by_key or nmatch > int(best_by_key[key][9]):
            best_by_key[key] = fields

    hits = {}
    for key, fields in best_by_key.items():
        qstart, qend = int(fields[2]), int(fields[3])
        strand = fields[4]
        tname, tstart, tend = fields[5], int(fields[7]), int(fields[8])
        seq = query_fasta.fetch(tname, tstart, tend)
        if strand == "-":
            seq = str(Seq(seq).reverse_complement())
        hit_span_diff = (tend - tstart) - (qend - qstart)
        hits[key] = (seq, tname, tstart, tend, strand, hit_span_diff)
        logger.debug(
            "%s: orthologous window %s:%d-%d(%s), span diff %+d",
            key, tname, tstart, tend, strand, hit_span_diff,
        )
    logger.info(
        "%d/%d breakpoint(s) had a minimap2 hit in the query genome", len(hits), len(records)
    )
    return hits


@dataclass
class PendingRecord:
    key: str
    sample: str
    bp: Breakpoint
    ref_window_seq: str
    ref_start: int
    ref_end: int


@dataclass
class AlignResult:
    score: float = float("nan")
    identity_pct: float = float("nan")
    aligned_len: int = 0
    target_start: int = -1
    target_end: int = -1


def local_align(query_seq, target_seq):
    if not query_seq or not target_seq:
        return AlignResult()
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    alignments = aligner.align(query_seq, target_seq)
    best = alignments[0]
    counts = best.counts()
    identities = counts.identities
    mismatches = counts.mismatches
    denom = identities + mismatches
    identity_pct = 100.0 * identities / denom if denom else float("nan")
    target_blocks = best.aligned[1]
    target_start = int(target_blocks[0][0]) if len(target_blocks) else -1
    target_end = int(target_blocks[-1][1]) if len(target_blocks) else -1
    aligned_len = target_end - target_start if target_end >= 0 else 0
    return AlignResult(
        score=best.score,
        identity_pct=identity_pct,
        aligned_len=aligned_len,
        target_start=target_start,
        target_end=target_end,
    )


def classify(
    right_ref,
    right_query,
    left_ref,
    left_query,
    min_identity,
    min_coverage,
    right_len,
    left_len,
):
    def hits(res, clip_len):
        if clip_len == 0:
            return None
        return (
            res.identity_pct >= min_identity
            and (res.aligned_len / clip_len) >= min_coverage
        )

    right_in_query = hits(right_query, right_len)
    left_in_query = hits(left_query, left_len)
    right_in_ref = hits(right_ref, right_len)
    left_in_ref = hits(left_ref, left_len)

    if right_in_ref or left_in_ref:
        return "inconclusive (clip also matches reference window — check window/repeat content)"
    if right_in_query and left_in_query:
        return "likely present in query genome"
    if right_in_query or left_in_query:
        return "ambiguous — only one clip end matches query genome"
    return "likely absent from query genome (consistent with novel/somatic)"


def main():
    parser = argparse.ArgumentParser(
        prog="check_strain_insertion.py",
        description=(
            "Check whether candidate insertion breakpoints from a "
            "collect_evidence.py .fa/.fa.gz file already exist in a second "
            "(query) genome, e.g. checking a C57BL/6-called breakpoint "
            "against CD-1."
            "Useful when your samples came from a hybbrid mouse strain."
        ),
    )
    parser.add_argument(
        "--reference-fasta",
        required=True,
        help="Genome fasta the evidence coordinates are based on (e.g. C57BL/6)",
    )
    parser.add_argument(
        "--query-fasta",
        required=True,
        help="Second genome fasta to check for pre-existing insertions (e.g. CD-1)",
    )
    evidence_group = parser.add_mutually_exclusive_group(required=True)
    evidence_group.add_argument(
        "--evidence-dir",
        help="Directory to scan for *.fa/*.fa.gz evidence files (from collect_evidence.py)",
    )
    evidence_group.add_argument(
        "--evidence-fa",
        nargs="+",
        help="One or more explicit evidence .fa/.fa.gz file paths (e.g. one per sample)",
    )
    parser.add_argument(
        "--same-scaffold",
        action="store_true",
        help=(
            "Reference and query fastas share the same coordinate system "
            "(true for some reference-anchored inbred-strain assemblies; not "
            "expected for an independently-assembled or outbred genome like "
            "CD-1). If not set (default), the orthologous locus in "
            "--query-fasta is located with minimap2 before extracting its "
            "window."
        ),
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10000,
        help=(
            "Number of bp to extract upstream and downstream of each "
            "breakpoint position (per side; total extracted window is "
            "~2x this). Should comfortably exceed the expected insertion "
            "size — IAP LTRs are typically 2-8kb. Default: 10000"
        ),
    )
    parser.add_argument("--output", required=True, help="Path to output summary tsv")
    parser.add_argument(
        "--min-identity",
        type=float,
        default=85.0,
        help="Identity %% threshold used only for the summary 'note' column. Default: 85.0",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.6,
        help=(
            "Minimum fraction of a clip sequence that must be aligned, used "
            "only for the summary 'note' column. Default: 0.6"
        ),
    )
    parser.add_argument(
        "--minimap2-bin",
        default="minimap2",
        help="minimap2 executable to use for orthologous-locus lookup. Default: minimap2",
    )
    parser.add_argument(
        "--minimap2-preset",
        default="asm5",
        help="minimap2 -x preset for orthologous-locus lookup. Default: asm5",
    )
    parser.add_argument(
        "--minimap2-threads",
        type=int,
        default=4,
        help="Threads for minimap2 index building and mapping. Default: 4",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity. Default: INFO",
    )
    parser.add_argument(
        "--log-file",
        help="Optional path to also write log entries to (in addition to stderr)",
    )
    args = parser.parse_args()

    log_format = logging.Formatter(
        "%(asctime)s %(levelname)-8s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    handlers = [logging.StreamHandler(sys.stderr)]
    if args.log_file:
        handlers.append(logging.FileHandler(args.log_file))
    for handler in handlers:
        handler.setFormatter(log_format)
    logging.basicConfig(level=args.log_level, handlers=handlers)

    logger.info("Run configuration:")
    for key, value in sorted(vars(args).items()):
        logger.info("  --%s = %s", key.replace("_", "-"), value)

    for path in (args.reference_fasta, args.query_fasta):
        if not os.path.isfile(path):
            logger.error("file not found: %s", path)
            sys.exit(1)

    if args.evidence_dir:
        evidence_files = find_evidence_files(args.evidence_dir)
    else:
        evidence_files = args.evidence_fa
        for path in evidence_files:
            if not os.path.isfile(path):
                logger.error("file not found: %s", path)
                sys.exit(1)
    logger.info("Found %d evidence file(s) to process", len(evidence_files))

    ref_fasta, _ = prepare_fasta(args.reference_fasta)
    query_fasta, query_fasta_resolved_path = prepare_fasta(args.query_fasta)

    minimap2_index = None
    if not args.same_scaffold:
        minimap2_index = build_minimap2_index(
            query_fasta_resolved_path, args.minimap2_bin, args.minimap2_preset,
            args.minimap2_threads,
        )

    rows = []
    windows_out_path = os.path.splitext(args.output)[0] + ".windows.fa"
    try:
        # Pass 1: parse every evidence file and extract each breakpoint's
        # reference window (cheap — indexed random access on an already-open
        # fasta).
        pending = []
        for evidence_path in evidence_files:
            sample = sample_name(evidence_path)
            logger.info("Reading sample '%s' (%s)", sample, evidence_path)
            for bp in parse_evidence_fa(evidence_path):
                key = f"{sample}_{bp.name}"
                try:
                    ref_window_seq, ref_start, ref_end = extract_window(
                        ref_fasta, bp.seqname, bp.pos1, args.window
                    )
                except (KeyError, ValueError) as e:
                    logger.warning(
                        "%s: skipping, couldn't read reference window: %s", key, e
                    )
                    continue
                pending.append(
                    PendingRecord(key, sample, bp, ref_window_seq, ref_start, ref_end)
                )
        logger.info("%d breakpoint(s) to check across %d sample(s)", len(pending), len(evidence_files))

        # Pass 2: locate every orthologous window in one batched minimap2
        # call against the pre-built index (skipped entirely for
        # --same-scaffold, where coordinates are just reused directly).
        batch_hits = {}
        if not args.same_scaffold:
            batch_hits = locate_orthologous_windows_batch(
                pending, minimap2_index, query_fasta,
                args.minimap2_bin, args.minimap2_preset, args.minimap2_threads,
            )

        # Pass 3: align clips against both windows and classify each breakpoint.
        with open(windows_out_path, "w") as windows_out:
            for rec in pending:
                bp = rec.bp
                logger.info("Processing %s", rec.key)
                windows_out.write(
                    f">{rec.key}_reference:{bp.seqname}:{rec.ref_start}-{rec.ref_end}\n{rec.ref_window_seq}\n"
                )

                query_scaffold = query_start = query_end = query_strand = (
                    hit_span_diff
                ) = None
                if args.same_scaffold:
                    try:
                        query_window_seq, query_start, query_end = extract_window(
                            query_fasta, bp.seqname, bp.pos1, args.window
                        )
                        query_scaffold, query_strand = bp.seqname, "+"
                    except (KeyError, ValueError):
                        query_window_seq = None
                else:
                    hit = batch_hits.get(rec.key)
                    if hit is None:
                        query_window_seq = None
                    else:
                        (
                            query_window_seq,
                            query_scaffold,
                            query_start,
                            query_end,
                            query_strand,
                            hit_span_diff,
                        ) = hit

                if query_window_seq:
                    windows_out.write(
                        f">{rec.key}_query:{query_scaffold}:{query_start}-{query_end}({query_strand})\n"
                        f"{query_window_seq}\n"
                    )

                window_vs_window = (
                    local_align(rec.ref_window_seq, query_window_seq)
                    if query_window_seq
                    else AlignResult()
                )

                right_clip = bp.longest_right_clip
                left_clip = bp.longest_left_clip
                right_vs_ref = local_align(right_clip, rec.ref_window_seq)
                left_vs_ref = local_align(left_clip, rec.ref_window_seq)
                right_vs_query = (
                    local_align(right_clip, query_window_seq)
                    if query_window_seq
                    else AlignResult()
                )
                left_vs_query = (
                    local_align(left_clip, query_window_seq)
                    if query_window_seq
                    else AlignResult()
                )

                note = (
                    classify(
                        right_vs_ref,
                        right_vs_query,
                        left_vs_ref,
                        left_vs_query,
                        args.min_identity,
                        args.min_coverage,
                        len(right_clip),
                        len(left_clip),
                    )
                    if query_window_seq
                    else "no orthologous window found"
                )
                logger.info("%s: %s", rec.key, note)

                rows.append(
                    {
                        "sample": rec.sample,
                        "seqname": bp.seqname,
                        "pos": bp.pos1,
                        "strand": bp.strand,
                        "right_clip_len": len(right_clip),
                        "left_clip_len": len(left_clip),
                        "query_scaffold": query_scaffold,
                        "query_start": query_start,
                        "query_end": query_end,
                        "query_strand": query_strand,
                        "orthology_hit_span_diff": hit_span_diff,
                        "window_identity_pct": window_vs_window.identity_pct,
                        "window_aligned_len": window_vs_window.aligned_len,
                        "right_clip_vs_ref_identity_pct": right_vs_ref.identity_pct,
                        "right_clip_vs_ref_aligned_len": right_vs_ref.aligned_len,
                        "right_clip_vs_query_identity_pct": right_vs_query.identity_pct,
                        "right_clip_vs_query_aligned_len": right_vs_query.aligned_len,
                        "right_clip_vs_query_pos": (
                            f"{right_vs_query.target_start}-{right_vs_query.target_end}"
                            if right_vs_query.target_start >= 0
                            else None
                        ),
                        "left_clip_vs_ref_identity_pct": left_vs_ref.identity_pct,
                        "left_clip_vs_ref_aligned_len": left_vs_ref.aligned_len,
                        "left_clip_vs_query_identity_pct": left_vs_query.identity_pct,
                        "left_clip_vs_query_aligned_len": left_vs_query.aligned_len,
                        "left_clip_vs_query_pos": (
                            f"{left_vs_query.target_start}-{left_vs_query.target_end}"
                            if left_vs_query.target_start >= 0
                            else None
                        ),
                        "note": note,
                    }
                )
    finally:
        for tmp_path in _TEMP_FILES:
            os.unlink(tmp_path)

    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    logger.info(
        "Wrote %d breakpoint(s) from %d sample(s) to %s and %s",
        len(rows), len(evidence_files), args.output, windows_out_path,
    )


if __name__ == "__main__":
    main()
