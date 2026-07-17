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

Dependencies beyond this repo's usual env.txt: biopython>=1.80 (pairwise
local alignment), and the minimap2 binary on PATH (only needed unless
--same-scaffold is given).
"""

import argparse
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


def locate_orthologous_window(
    ref_window_seq, bp_name, minimap2_target_path, query_fasta, minimap2_bin, preset
):
    """Find the query-genome window orthologous to ref_window_seq via minimap2.

    Returns (seq, scaffold, start, end, strand, hit_span_diff) or None if no
    hit was found. hit_span_diff = (target span) - (query span) for the best
    hit: a value close to the expected insertion size is itself a signal the
    insertion is present in the query genome (a real minimap2 gap spanning
    the missing sequence), independent of the clip alignments below.
    """
    logger.debug("Running minimap2 (-x %s) for %s against %s", preset, bp_name, minimap2_target_path)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp.write(f">{bp_name}\n{ref_window_seq}\n")
        tmp_path = tmp.name
    try:
        try:
            result = subprocess.run(
                [
                    minimap2_bin,
                    "-x",
                    preset,
                    "--secondary=no",
                    minimap2_target_path,
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
            logger.warning("minimap2 failed for %s: %s", bp_name, e.stderr)
            return None
    finally:
        os.unlink(tmp_path)

    lines = [l for l in result.stdout.splitlines() if l.strip()]
    if not lines:
        logger.debug("%s: no minimap2 hit against the query genome", bp_name)
        return None
    best = max(lines, key=lambda l: int(l.split("\t")[9]))  # nmatch
    fields = best.split("\t")
    tname, tstart, tend, strand = fields[5], int(fields[7]), int(fields[8]), fields[4]
    qstart, qend = int(fields[2]), int(fields[3])

    seq = query_fasta.fetch(tname, tstart, tend)
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    hit_span_diff = (tend - tstart) - (qend - qstart)
    logger.debug(
        "%s: orthologous window %s:%d-%d(%s), span diff %+d",
        bp_name, tname, tstart, tend, strand, hit_span_diff,
    )
    return seq, tname, tstart, tend, strand, hit_span_diff


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
    parser.add_argument(
        "--evidence-fa",
        required=True,
        help="Path to the .fa/.fa.gz file produced by collect_evidence.py",
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
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity. Default: INFO",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )

    for path in (args.reference_fasta, args.query_fasta, args.evidence_fa):
        if not os.path.isfile(path):
            logger.error("file not found: %s", path)
            sys.exit(1)

    ref_fasta, _ = prepare_fasta(args.reference_fasta)
    query_fasta, query_fasta_resolved_path = prepare_fasta(args.query_fasta)

    rows = []
    windows_out_path = os.path.splitext(args.output)[0] + ".windows.fa"
    try:
        with open(windows_out_path, "w") as windows_out:
            for bp in parse_evidence_fa(args.evidence_fa):
                logger.info("Processing %s", bp.name)
                try:
                    ref_window_seq, ref_start, ref_end = extract_window(
                        ref_fasta, bp.seqname, bp.pos1, args.window
                    )
                except (KeyError, ValueError) as e:
                    logger.warning("%s: skipping, couldn't read reference window: %s", bp.name, e)
                    continue
                logger.debug(
                    "%s: reference window %s:%d-%d (%dbp)",
                    bp.name, bp.seqname, ref_start, ref_end, ref_end - ref_start,
                )
                windows_out.write(
                    f">{bp.name}_reference:{bp.seqname}:{ref_start}-{ref_end}\n{ref_window_seq}\n"
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
                    hit = locate_orthologous_window(
                        ref_window_seq,
                        bp.name,
                        query_fasta_resolved_path,
                        query_fasta,
                        args.minimap2_bin,
                        args.minimap2_preset,
                    )
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
                        f">{bp.name}_query:{query_scaffold}:{query_start}-{query_end}({query_strand})\n"
                        f"{query_window_seq}\n"
                    )

                window_vs_window = (
                    local_align(ref_window_seq, query_window_seq)
                    if query_window_seq
                    else AlignResult()
                )

                right_clip = bp.longest_right_clip
                left_clip = bp.longest_left_clip
                right_vs_ref = local_align(right_clip, ref_window_seq)
                left_vs_ref = local_align(left_clip, ref_window_seq)
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
                logger.info("%s: %s", bp.name, note)

                rows.append(
                    {
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
        "Wrote %d breakpoint(s) to %s and %s", len(rows), args.output, windows_out_path
    )


if __name__ == "__main__":
    main()
