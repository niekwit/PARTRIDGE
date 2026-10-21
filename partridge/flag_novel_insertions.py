"""
Automate the step-5 checks (README rules 1-5) for every candidate insertion
site in the audit table from collect_insertions.R (step 3), and flag sites
that pass all of them as possible novel insertions.

No site is removed: every row of the audit table is written back out, with
one set of columns per rule and a final "possible_novel_insertion" flag.
The flag is a pre-screen to prioritise manual validation, not a replacement
for it; the "classification" column used by 06_cluster_source_coverage.R is
left for you to fill in.

Input
-----
--audit-xlsx  audit table from collect_insertions.R ("name" = seqname:first_bp,
              "seqname", "strand", "first_bp", "mouse").
--bam/--sample  one or more coordinate-sorted, indexed BAMs with the Expt.ID
              of each (same order). Rows whose "mouse" matches none of the
              samples (e.g. "multi") are kept but not evaluated.
--genome      indexed reference FASTA (samtools faidx), used for the TSD
              sequence and, unless --blat-db is given, as the BLAT database.
--rmsk        optional RepeatMasker annotation, either UCSC rmsk.txt(.gz) or
              a BED file (name in column 4). Without it rule 2 is not
              evaluated.

Method, per site
----------------
Clipped reads are taken from the BAM rather than from the collect_evidence.py
evidence file, because the TSD (rule 4) needs the junction coordinates, which
that file does not keep. All mapped, non-secondary, non-QC-fail records are
used, including supplementary alignments (in practice most junction-spanning
reads are supplementary). A soft clip on either end of a read counts when it
is at least --min-clip-len bp long and its junction lies within --window bp of
first_bp:
  - right clip: junction = last aligned base (1-based), clip = bases after it;
  - left clip: junction = first aligned base (1-based), clip = bases before it.
For each side the most frequent junction is the breakpoint, and a majority-
vote consensus of the clips at that breakpoint, anchored at the junction, is
the LTR-end sequence. As in collect_evidence.py, all sequences are on the
reference forward strand.

Rule 1  LTR coverage: at least one concordant clip longer than
        --min-ltr-len (20) bp on both the right and the left side.
Rule 2  Not in a repeat: no RepeatMasker element overlaps the TSD extended by
        --repeat-flank bp on each side. Overlapping repeats are listed so you
        can judge whether they are "obviously repetitive". NA without --rmsk.
Rule 3  Concordant LTR start: every clip near the site starts the LTR at the
        same place. A clip whose junction is d bp from the breakpoint (the
        aligner clipped or extended a few bases) is concordant if |d| <=
        --max-shift and, after correcting for d, its first bases next to the
        junction match the consensus (<= 10% mismatches over up to 20 bp).
Rule 4  TSD: the right breakpoint lies inside or past the left one, so the
        duplicated stretch [left breakpoint, right breakpoint] is between
        --min-tsd and --max-tsd bp (defaults: more than 2, less than 30).
        Note: a base at the LTR terminus that happens to match the genome is
        absorbed into the alignment, so the TSD can be 1-2 bp longer than
        the "true" one (the test site gives 7 bp, AGGAGGC, where the README
        counts the 6 bp AGGAGG).
Rule 5  Source LTR pair: the right-clip consensus followed directly by the
        left-clip consensus (the query pasted into UCSC BLAT in the manual
        protocol) is aligned to the genome with BLAT. Hits with identity >
        --min-identity (93%, computed as on the UCSC BLAT web page via
        pslCalcMilliBad) that cover >= --min-query-cov of the query qualify.
        Among the best-scoring qualifying hits, each one with another
        qualifying hit on the same chromosome and strand starting
        --pair-min to --pair-max bp (2-8 kb) away is a candidate source:
        both LTRs of one full-length element. "source" is the first
        candidate (span of both LTRs, 1-based), "source_candidates" lists
        all, "source_unique" says whether there is exactly one.

possible_novel_insertion is TRUE when every evaluated rule passes;
"failed_rules" and "rules_not_evaluated" say why it is not.

Output
------
The audit table with extra columns, written as .xlsx or tab-separated text
depending on the --output extension.

Dependencies beyond this repo's usual environment: the BLAT binary
(bioconda::blat) on PATH or given with --blat. Standalone BLAT without an
.ooc file is used on purpose: an .ooc file masks over-represented 11-mers,
which include IAP LTR sequence.
"""

import argparse
import gzip
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass, field

import pandas as pd
import pysam

logger = logging.getLogger("flag_novel_insertions")

PSL_COLUMNS = [
    "matches",
    "misMatches",
    "repMatches",
    "nCount",
    "qNumInsert",
    "qBaseInsert",
    "tNumInsert",
    "tBaseInsert",
    "strand",
    "qName",
    "qSize",
    "qStart",
    "qEnd",
    "tName",
    "tSize",
    "tStart",
    "tEnd",
    "blockCount",
    "blockSizes",
    "qStarts",
    "tStarts",
]


@dataclass
class Clip:
    """
    One soft-clipped read end at a candidate insertion site.

    Attributes
    ----------
    junction : int
        1-based reference position of the last aligned base before the clip
        (right clips) or the first aligned base after it (left clips).
    seq : str
        The clipped bases on the reference forward strand, oriented so that
        seq[0] is the base next to the junction (left clips are reversed).
    """

    junction: int  # 1-based reference position of the aligned base at the junction
    seq: str  # clipped bases, oriented so that seq[0] is next to the junction


@dataclass
class Side:
    """
    All clips on one side (right or left) of a site and the summary that
    evaluate_side() derives from them.

    Attributes
    ----------
    clips : list of Clip
        Clips on this side, sequences oriented junction-first.
    breakpoint : int or None
        Most frequent junction position (1-based); None without clips.
    consensus : str
        Majority-vote consensus of the clips at the breakpoint, junction-first.
    n_concordant, n_discordant : int
        Clips that do / do not start the LTR at the breakpoint (rule 3).
    max_ltr_len : int
        Longest LTR stretch (bp) covered by one concordant clip (rule 1).
    """

    clips: list = field(default_factory=list)
    breakpoint: int = None
    consensus: str = ""
    n_concordant: int = 0
    n_discordant: int = 0
    max_ltr_len: int = 0


def parse_args():
    p = argparse.ArgumentParser(
        prog="flag_novel_insertions.py",
        description=(
            "Apply the step-5 validation rules to every candidate insertion "
            "site and flag possible novel insertions (no site is removed)."
        ),
    )
    p.add_argument("--audit-xlsx", required=True, help="isa.audit.xlsx from step 3")
    p.add_argument("--bam", required=True, nargs="+", help="indexed BAM(s)")
    p.add_argument(
        "--sample",
        required=True,
        nargs="+",
        help="Expt.ID of each BAM, same order, matched to the 'mouse' column",
    )
    p.add_argument("--genome", required=True, help="faidx-indexed reference FASTA")
    p.add_argument("--output", required=True, help="output .xlsx, .tsv or .txt")
    p.add_argument("--rmsk", help="UCSC rmsk.txt(.gz) or BED of repeats (rule 2)")
    p.add_argument(
        "--blat-db",
        help="BLAT database (.2bit recommended, much faster); default: --genome",
    )
    p.add_argument("--blat", default="blat", help="BLAT binary (default: on PATH)")
    p.add_argument(
        "--blat-args",
        default="-minIdentity=90",
        help="extra BLAT arguments, one string (default: %(default)s)",
    )
    p.add_argument("--window", type=int, default=30, help="bp around first_bp")
    p.add_argument("--min-clip-len", type=int, default=5)
    p.add_argument("--min-ltr-len", type=int, default=20, help="rule 1, exclusive")
    p.add_argument("--repeat-flank", type=int, default=0, help="rule 2")
    p.add_argument("--max-shift", type=int, default=5, help="rule 3")
    p.add_argument("--min-tsd", type=int, default=3, help="rule 4, inclusive")
    p.add_argument("--max-tsd", type=int, default=29, help="rule 4, inclusive")
    p.add_argument(
        "--min-identity", type=float, default=93.0, help="rule 5, %%, exclusive"
    )
    p.add_argument("--min-query-cov", type=float, default=0.9, help="rule 5")
    p.add_argument("--pair-min", type=int, default=2000, help="rule 5, bp")
    p.add_argument("--pair-max", type=int, default=8000, help="rule 5, bp")
    p.add_argument("--max-candidates", type=int, default=20, help="listed sources")
    args = p.parse_args()
    if len(args.bam) != len(args.sample):
        p.error("--bam and --sample need the same number of values")
    return args


def soft_clips(read):
    """
    Soft-clip lengths at both ends of an aligned read.

    Hard clips are skipped, so a soft clip inside a hard clip (e.g. 5H10S...)
    is still found.

    Parameters
    ----------
    read : pysam.AlignedSegment
        A mapped read.

    Returns
    -------
    tuple of int
        (left, right) soft-clip lengths in bp; 0 where the end is not soft
        clipped.
    """
    ops = [(op, n) for op, n in read.cigartuples if op != pysam.CHARD_CLIP]
    left = ops[0][1] if ops[0][0] == pysam.CSOFT_CLIP else 0
    right = ops[-1][1] if ops[-1][0] == pysam.CSOFT_CLIP else 0
    return left, right


def collect_clips(bam, seqname, first_bp, window, min_clip_len):
    """
    Collect the soft clips whose junction lies near a candidate site.

    Unmapped, secondary and QC-fail records are skipped; supplementary
    alignments are kept. Each end of a read is considered separately, so a
    read clipped on both ends can contribute a clip on either side.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Indexed BAM of the sample.
    seqname : str
        Chromosome of the site.
    first_bp : int
        1-based site position from the audit table.
    window : int
        Only clips with a junction within first_bp +/- window are kept.
    min_clip_len : int
        Minimum clip length in bp.

    Returns
    -------
    tuple of Side
        (right, left): right clips (bases after the alignment end) and left
        clips (bases before the alignment start), both junction-first.
    """
    right, left = Side(), Side()
    lo, hi = first_bp - window, first_bp + window
    for r in bam.fetch(seqname, max(0, lo - 1), hi):
        if r.is_unmapped or r.is_secondary or r.is_qcfail or r.query_sequence is None:
            continue
        lclip, rclip = soft_clips(r)
        if rclip >= min_clip_len and lo <= r.reference_end <= hi:
            right.clips.append(Clip(r.reference_end, r.query_sequence[-rclip:]))
        if lclip >= min_clip_len and lo <= r.reference_start + 1 <= hi:
            # reversed so that, as for right clips, seq[0] touches the junction
            left.clips.append(
                Clip(r.reference_start + 1, r.query_sequence[:lclip][::-1])
            )
    return right, left


def majority_consensus(seqs):
    """
    Column-wise majority-vote consensus of sequences anchored at position 0.

    Each column is called from the sequences long enough to reach it, so the
    consensus is as long as the longest sequence. Ties go to the base seen
    first.

    Parameters
    ----------
    seqs : list of str
        Non-empty list of sequences aligned at their first base.

    Returns
    -------
    str
        The consensus sequence.
    """
    length = max(len(s) for s in seqs)
    cons = []
    for i in range(length):
        counts = Counter(s[i] for s in seqs if i < len(s))
        cons.append(counts.most_common(1)[0][0])
    return "".join(cons)


def junction_matches(a, b, k=20):
    """
    Test whether two junction-first sequences agree next to the junction.

    Compares the first min(len(a), len(b), k) bases position by position
    (no gaps) and allows up to 10% mismatches, at least 1.

    Parameters
    ----------
    a, b : str
        Sequences whose first base is the base next to the junction.
    k : int
        Maximum number of bases compared.

    Returns
    -------
    bool
        True if the sequences match; also True if there is nothing to
        compare (an empty sequence).
    """
    n = min(len(a), len(b), k)
    if n == 0:
        return True
    mismatches = sum(x != y for x, y in zip(a[:n], b[:n]))
    return mismatches <= max(1, int(0.1 * n))


def evaluate_side(side, first_bp, max_shift, toward_genome):
    """
    Derive the breakpoint, LTR consensus and clip concordance of one side.

    The breakpoint is the most frequent junction (ties: closest to first_bp)
    and the consensus is built from the clips at the breakpoint only. Every
    clip is then compared with the consensus after correcting for its shift
    d from the breakpoint (rule 3); concordant clips also set max_ltr_len
    (rule 1). Results are written into `side`; nothing happens without clips.

    Parameters
    ----------
    side : Side
        Clips of one side, from collect_clips(); updated in place.
    first_bp : int
        1-based site position from the audit table, to break ties.
    max_shift : int
        Largest |d| (bp) a clip may be shifted and still be concordant.
    toward_genome : int
        Sign of a junction shift that puts extra genomic bases into the
        clip: -1 for right clips (read clipped before the breakpoint), +1
        for left clips (read clipped after it).
    """
    if not side.clips:
        return
    counts = Counter(c.junction for c in side.clips)
    top = max(counts.values())
    side.breakpoint = min(
        (j for j, n in counts.items() if n == top), key=lambda j: abs(j - first_bp)
    )
    side.consensus = majority_consensus(
        [c.seq for c in side.clips if c.junction == side.breakpoint]
    )
    for c in side.clips:
        # d > 0: d genomic bases precede the LTR in the clip;
        # d < 0: the first -d LTR bases were aligned to the genome.
        d = (c.junction - side.breakpoint) * toward_genome
        if abs(d) > max_shift or len(c.seq) <= d:
            side.n_discordant += 1
            continue
        a, b = (c.seq[d:], side.consensus) if d >= 0 else (c.seq, side.consensus[-d:])
        if junction_matches(a, b):
            side.n_concordant += 1
            side.max_ltr_len = max(side.max_ltr_len, len(c.seq) - d)
        else:
            side.n_discordant += 1


def read_repeats(path, seqnames):
    """
    Load a RepeatMasker annotation for the chromosomes of interest (rule 2).

    The format is detected from the first line: UCSC rmsk.txt (>= 13
    columns, numeric "bin" first) or BED (chrom, start, end, name).

    Parameters
    ----------
    path : str
        UCSC rmsk.txt or BED file, optionally gzipped.
    seqnames : set of str
        Chromosomes to keep.

    Returns
    -------
    dict of str to pandas.DataFrame
        Per chromosome, a table with 0-based half-open "start"/"end" and a
        "name" (repName/repClass for rmsk.txt, column 4 for BED).
    """
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
    if len(first) >= 13 and first[0].isdigit():  # UCSC rmsk.txt
        cols, names = [5, 6, 7, 10, 11], ["chrom", "start", "end", "name", "rclass"]
    else:  # BED
        cols, names = [0, 1, 2, 3], ["chrom", "start", "end", "name"]
    rep = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=cols,
        names=names,
        comment="#",
        dtype={"chrom": str},
    )
    rep = rep[rep["chrom"].isin(seqnames)]
    if "rclass" in rep:
        rep["name"] = rep["name"] + "/" + rep["rclass"]
    return {c: g for c, g in rep.groupby("chrom")}


def psl_milli_bad(h):
    """
    Per-mille "badness" of a PSL hit, as pslCalcMilliBad(psl, isMrna=TRUE) in
    the UCSC kent source.

    The UCSC BLAT web page reports identity as 100 - milliBad * 0.1; this
    reproduces it so the threshold in rule 5 means the same as in the manual
    protocol. Mismatches, query inserts and a log penalty for the query/target
    size difference count as bad; target inserts (introns) do not.

    Parameters
    ----------
    h : pandas.Series
        One PSL row (columns as in PSL_COLUMNS).

    Returns
    -------
    float
        Bad bases per 1000 aligned bases; 0 for an empty alignment.
    """
    q_ali = h.qEnd - h.qStart
    t_ali = h.tEnd - h.tStart
    if min(q_ali, t_ali) <= 0:
        return 0
    size_dif = max(q_ali - t_ali, 0)
    total = h.matches + h.repMatches + h.misMatches
    if total == 0:
        return 0
    return (
        1000 * (h.misMatches + h.qNumInsert + round(3 * math.log(1 + size_dif)))
    ) / total


def run_blat(blat, db, queries, extra_args):
    """
    Align all queries to the genome with standalone BLAT in a single run.

    One run loads the database once, which dominates the run time for a
    whole genome. Queries and output go to a temporary directory.

    Parameters
    ----------
    blat : str
        Path to the BLAT binary.
    db : str
        BLAT database (.2bit or FASTA).
    queries : dict of str to str
        Query name (audit-table index as str) to sequence.
    extra_args : str
        Additional BLAT arguments, whitespace-separated.

    Returns
    -------
    pandas.DataFrame
        One row per PSL hit (PSL_COLUMNS) plus "identity" (%, as on the UCSC
        BLAT web page), "score" (pslScore: matches + repMatches/2 - mismatches
        - query and target inserts) and "query_cov" (aligned fraction of the
        query).

    Raises
    ------
    subprocess.CalledProcessError
        If BLAT exits with an error.
    """
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "queries.fa")
        psl = os.path.join(tmp, "hits.psl")
        with open(fa, "w") as fh:
            for name, seq in queries.items():
                fh.write(f">{name}\n{seq}\n")
        cmd = [blat, db, fa, psl, "-noHead", *extra_args.split()]
        logger.info("running: %s", " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        hits = pd.read_csv(psl, sep="\t", header=None, names=PSL_COLUMNS)
    hits["qName"] = hits["qName"].astype(str)
    hits["identity"] = 100.0 - hits.apply(psl_milli_bad, axis=1) * 0.1
    hits["score"] = (
        hits.matches
        + hits.repMatches // 2
        - hits.misMatches
        - hits.qNumInsert
        - hits.tNumInsert
    )
    hits["query_cov"] = (hits.qEnd - hits.qStart) / hits.qSize
    return hits


def source_pairs(hits, args):
    """
    Find candidate source elements (LTR pairs) among one query's BLAT hits.

    Hits qualify with identity > args.min_identity and query coverage >=
    args.min_query_cov. Each best-scoring qualifying hit is paired with every
    qualifying hit on the same chromosome and strand whose start lies
    args.pair_min to args.pair_max bp away, i.e. the other LTR of the same
    full-length element (rule 5).

    Parameters
    ----------
    hits : pandas.DataFrame
        BLAT hits of one query, from run_blat().
    args : argparse.Namespace
        Parsed arguments (min_identity, min_query_cov, pair_min, pair_max).

    Returns
    -------
    tuple
        (qualifying hits as a DataFrame, sorted list of unique pairs as
        (chrom, start, end, strand) with 1-based start and end spanning both
        LTRs).
    """
    ok = hits[
        (hits.identity > args.min_identity) & (hits.query_cov >= args.min_query_cov)
    ]
    if ok.empty:
        return ok, []
    best = ok[ok.score == ok.score.max()]
    pairs = set()
    for _, h in best.iterrows():
        partners = ok[
            (ok.tName == h.tName)
            & (ok.strand == h.strand)
            & ((ok.tStart - h.tStart).abs().between(args.pair_min, args.pair_max))
        ]
        for _, m in partners.iterrows():
            pairs.add(
                (h.tName, min(h.tStart, m.tStart) + 1, max(h.tEnd, m.tEnd), h.strand)
            )
    return ok, sorted(pairs)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    blat = shutil.which(args.blat)
    if blat is None:
        sys.exit(f"ERROR: BLAT binary not found: {args.blat}")
    for path in args.bam + [args.genome]:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: file not found: {path}")

    audit = pd.read_excel(args.audit_xlsx)
    bam_of = dict(zip(args.sample, args.bam))
    genome = pysam.FastaFile(args.genome)
    repeats = (
        read_repeats(args.rmsk, set(audit["seqname"].astype(str)))
        if args.rmsk
        else None
    )

    results = {}  # audit index -> dict of new columns
    queries = {}  # audit index (as str) -> BLAT query
    for sample, bam_path in bam_of.items():
        rows = audit[audit["mouse"].astype(str) == sample]
        logger.info("%s: %d sites", sample, len(rows))
        with pysam.AlignmentFile(bam_path) as bam:
            for idx, row in rows.iterrows():
                seqname, first_bp = str(row["seqname"]), int(row["first_bp"])
                right, left = collect_clips(
                    bam, seqname, first_bp, args.window, args.min_clip_len
                )
                evaluate_side(right, first_bp, args.max_shift, toward_genome=-1)
                evaluate_side(left, first_bp, args.max_shift, toward_genome=+1)
                res = {
                    "evaluated": True,
                    "n_right_clips": len(right.clips),
                    "n_left_clips": len(left.clips),
                    "right_breakpoint": right.breakpoint,
                    "left_breakpoint": left.breakpoint,
                    "right_ltr_max_len": right.max_ltr_len,
                    "left_ltr_max_len": left.max_ltr_len,
                    "rule1_ltr_coverage": right.max_ltr_len > args.min_ltr_len
                    and left.max_ltr_len > args.min_ltr_len,
                }

                # rule 4 (TSD) first: rule 2 uses its interval
                tsd_start, tsd_end = left.breakpoint, right.breakpoint
                if tsd_start is not None and tsd_end is not None:
                    tsd_len = tsd_end - tsd_start + 1
                    res["tsd"] = f"{seqname}:{tsd_start}-{tsd_end}"
                    res["tsd_len"] = tsd_len
                    res["tsd_seq"] = (
                        genome.fetch(seqname, tsd_start - 1, tsd_end).upper()
                        if tsd_len > 0
                        else ""
                    )
                    res["rule4_tsd"] = args.min_tsd <= tsd_len <= args.max_tsd
                else:
                    res["rule4_tsd"] = False

                if repeats is not None:
                    lo = (tsd_start if tsd_start is not None else first_bp) - 1
                    hi = tsd_end if tsd_end is not None else first_bp
                    lo, hi = lo - args.repeat_flank, hi + args.repeat_flank
                    rep = repeats.get(seqname)
                    hit = (
                        rep[(rep.start < hi) & (rep.end > lo)]
                        if rep is not None
                        else []
                    )
                    res["overlapping_repeats"] = (
                        ";".join(hit["name"]) if len(hit) else ""
                    )
                    res["rule2_not_repetitive"] = len(hit) == 0

                res["n_discordant_clips"] = right.n_discordant + left.n_discordant
                res["rule3_concordant_ltr_start"] = (
                    bool(right.clips and left.clips) and res["n_discordant_clips"] == 0
                )

                if right.consensus and left.consensus:
                    # left consensus is stored junction-first; restore + strand
                    q = right.consensus + left.consensus[::-1]
                    res["blat_query"] = q
                    queries[str(idx)] = q
                res["rule5_source_ltr_pair"] = False
                results[idx] = res

    if queries:
        hits = run_blat(blat, args.blat_db or args.genome, queries, args.blat_args)
        for name, h in hits.groupby("qName"):
            ok, pairs = source_pairs(h, args)
            res = results[int(name)]
            if not ok.empty:
                top = ok.sort_values(["score", "identity"], ascending=False).iloc[0]
                res["best_hit"] = f"{top.tName}:{top.tStart + 1}-{top.tEnd}"
                res["best_hit_identity"] = round(top.identity, 1)
                res["best_hit_query_cov"] = round(top.query_cov, 3)
            res["n_source_candidates"] = len(pairs)
            if pairs:
                listed = [f"{c}:{s}-{e}({st})" for c, s, e, st in pairs]
                res["source"] = f"{pairs[0][0]}:{pairs[0][1]}-{pairs[0][2]}"
                res["source_candidates"] = ";".join(listed[: args.max_candidates])
                res["source_unique"] = len(pairs) == 1
                res["rule5_source_ltr_pair"] = True

    rules = {
        "1": "rule1_ltr_coverage",
        "2": "rule2_not_repetitive",
        "3": "rule3_concordant_ltr_start",
        "4": "rule4_tsd",
        "5": "rule5_source_ltr_pair",
    }
    for res in results.values():
        failed = [k for k, col in rules.items() if col in res and not res[col]]
        res["failed_rules"] = ",".join(failed)
        res["rules_not_evaluated"] = ",".join(
            k for k, c in rules.items() if c not in res
        )
        res["possible_novel_insertion"] = not failed

    new = pd.DataFrame.from_dict(results, orient="index").convert_dtypes()
    out = audit.join(new)
    # rows not evaluated (no BAM for their "mouse") are never flagged
    for col in ("evaluated", "possible_novel_insertion"):
        out[col] = out[col].fillna(False).astype(bool)
    if args.output.endswith(".xlsx"):
        out.to_excel(args.output, index=False)
    else:
        out.to_csv(args.output, sep="\t", index=False)
    logger.info(
        "%d of %d sites evaluated, %d flagged as possible novel insertions -> %s",
        out["evaluated"].sum(),
        len(out),
        out["possible_novel_insertion"].sum(),
        args.output,
    )


if __name__ == "__main__":
    main()
