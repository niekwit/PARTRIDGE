"""
Quality-control check: how much of a hybridisation-capture BAM's signal is
landing on known/annotated ERVK elements, versus everything sequenced?

Rationale: reads whose CIGAR fully matches the reference (no soft/hard
clip) within a small mismatch tolerance are, by construction, NOT the
clipped junction-spanning reads isa_hmmer2.py looks for — they represent
ordinary reads from loci already present in the reference genome. For this
assay that mostly means reads from known, reference-annotated ERVK/IAP
copies. Intersecting those reads with a BED of ERVK element coordinates and
expressing the result as a fraction of all reads in the BAM gives a simple
enrichment metric: what fraction of everything sequenced for this sample
is a clean hit on a known element we're trying to capture.

For each *.bam file in --bam-dir, reports:
  - total_reads:                     primary reads (excludes secondary/supplementary)
  - fully_matched_reads:             of those, reads with no soft/hard clip
                                      and <= --max-mismatches (NM tag)
  - ervk_overlapping_reads:          of the fully-matched reads, how many
                                      overlap >=1 interval in --bed
  - pct_total_reads_overlapping_ervk: ervk_overlapping_reads / total_reads * 100

Writes one combined --output-csv across all BAMs, and a --output-plot bar
chart (seaborn) of the percentage per sample.

Note: BAM and BED must already agree on chromosome-naming convention (e.g.
both Ensembl-style "1" or both UCSC-style "chr1") — this script does not
reconcile a mismatch between them (see partridge/check_strain_insertion.py's
resolve_seqname for that problem in a different context).

Dependencies beyond this repo's usual environment.yaml: pybedtools (and the
bedtools binary it wraps), seaborn.
"""

import argparse
import glob
import logging
import os
import sys
import tempfile

import pandas as pd
import pybedtools
import pysam
import seaborn as sns
from matplotlib import pyplot as plt
from pysam import CHARD_CLIP, CSOFT_CLIP

logger = logging.getLogger("ervk_enrichment_qc")


def find_bam_files(directory):
    files = sorted(glob.glob(os.path.join(directory, "*.bam")))
    if not files:
        logger.error("No *.bam files found in %s", directory)
        sys.exit(1)
    return files


def sample_name(bam_path):
    base = os.path.basename(bam_path)
    return base[: -len(".bam")] if base.endswith(".bam") else base


def fully_matches_reference(read, max_mismatches):
    """True if `read` has no soft/hard clip (every base is aligned) and its
    edit distance to the reference (NM tag) is within max_mismatches."""
    cigar = read.cigartuples
    if not cigar:
        return False
    if any(op in (CSOFT_CLIP, CHARD_CLIP) for op, _ in cigar):
        return False
    try:
        nm = read.get_tag("NM")
    except KeyError:
        return False
    return nm <= max_mismatches


def analyze_bam(bam_path, bed, max_mismatches, threads):
    total_reads = 0
    fully_matched = 0
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    try:
        with tmp, pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
            for read in bam:
                if read.is_secondary or read.is_supplementary:
                    continue
                total_reads += 1
                if read.is_unmapped:
                    continue
                if not fully_matches_reference(read, max_mismatches):
                    continue
                fully_matched += 1
                tmp.write(f"{read.reference_name}\t{read.reference_start}\t{read.reference_end}\n")

        overlapping = 0
        if fully_matched:
            qualifying = pybedtools.BedTool(tmp.name)
            overlapping = qualifying.intersect(bed, u=True).count()
    finally:
        os.unlink(tmp.name)

    pct = 100.0 * overlapping / total_reads if total_reads else float("nan")
    return total_reads, fully_matched, overlapping, pct


def main():
    parser = argparse.ArgumentParser(
        prog="ervk_enrichment_qc.py",
        description=(
            "For every BAM in --bam-dir, report what fraction of all reads "
            "are non-clipped, low-mismatch hits that overlap a known ERVK "
            "element in --bed — a quality-control check for enrichment of "
            "known/annotated insertions."
        ),
    )
    parser.add_argument(
        "--bam-dir",
        required=True,
        help="Directory of *.bam files (output of make_bams.sh)",
    )
    parser.add_argument("--bed", required=True, help="Path to the ERVK elements BED file")
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=2,
        help=(
            "Maximum edit distance (NM tag) allowed for a read to count as "
            "fully matching the reference, to allow for sequencing errors. "
            "Default: 2"
        ),
    )
    parser.add_argument("--output-csv", required=True, help="Path to output summary csv")
    parser.add_argument("--output-plot", required=True, help="Path to output bar chart (e.g. .png/.pdf)")
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for BAM decompression (per file). Default: 1",
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

    if not os.path.isdir(args.bam_dir):
        logger.error("directory not found: %s", args.bam_dir)
        sys.exit(1)
    if not os.path.isfile(args.bed):
        logger.error("file not found: %s", args.bed)
        sys.exit(1)

    bam_files = find_bam_files(args.bam_dir)
    logger.info("Found %d bam file(s) to process", len(bam_files))
    bed = pybedtools.BedTool(args.bed)

    rows = []
    for bam_path in bam_files:
        sample = sample_name(bam_path)
        logger.info("Processing %s (%s)", sample, bam_path)
        total, matched, overlapping, pct = analyze_bam(
            bam_path, bed, args.max_mismatches, args.threads
        )
        logger.info(
            "%s: %d total reads, %d fully-matched (<=%d mismatches, no clip), "
            "%d overlap an ERVK element (%.3f%% of total)",
            sample, total, matched, args.max_mismatches, overlapping, pct,
        )
        rows.append(
            {
                "sample": sample,
                "total_reads": total,
                "fully_matched_reads": matched,
                "ervk_overlapping_reads": overlapping,
                "pct_total_reads_overlapping_ervk": pct,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(args.output_csv, index=False)
    logger.info("Wrote %s", args.output_csv)

    plt.figure(figsize=(max(6, 0.6 * len(df)), 5))
    ax = sns.barplot(
        data=df, x="sample", y="pct_total_reads_overlapping_ervk", color="#2a78d6"
    )
    ax.set_ylabel("% of total reads overlapping a known ERVK element")
    ax.set_xlabel("")
    ax.set_title("ERVK enrichment QC")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(args.output_plot, dpi=150)
    logger.info("Wrote %s", args.output_plot)


if __name__ == "__main__":
    main()
