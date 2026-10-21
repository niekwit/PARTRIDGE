"""
Collect the read evidence for each candidate insertion site of one sample
(PARTRIDGE step 4), for manual validation in step 5 (README rules 1-5,
e.g. BLAT of the clipped sequences to find the source IAP).

Input
-----
--audit-xlsx  audit table from collect_insertions.R (step 3). Only rows whose
              "mouse" column equals --sample are processed; this column holds
              the sample's Expt.ID. Rows with mouse == "multi" (sites found in
              more than one sample) never match and are skipped.
--bam         coordinate-sorted, indexed BAM of that sample.
--sample      the sample's Expt.ID.
--output      gzipped text file to write (.fa.gz).

Method
------
For each site ("name" = "seqname:first_bp", and "strand"):

1. Clipped reads. Fetch all reads overlapping first_bp .. first_bp + 7
   (1-based). The site is assumed to span first_bp .. first_bp + 6, i.e. a
   6 bp target-site duplication; the real TSD length is not used. Keep reads
   that are soft-clipped on exactly one side (reads clipped on both sides or
   not at all are skipped):
   - left-clipped reads: the clipped bases (before the alignment start) go to
     LEFT_CLIP; for reverse-strand reads, the mate lies to the left, i.e.
     inside the inserted element, so the read name is kept to find the mate.
   - right-clipped reads: the clipped bases (after the alignment end) go to
     RIGHT_CLIP; for forward-strand reads the mate lies to the right, inside
     the element, and the read name is kept.
   Clips are not required to be at the breakpoint or to match an LTR: any
   one-sided clip overlapping the window is reported. Duplicate, secondary
   and supplementary alignments are not filtered out.
   Mates are dropped when the mate is itself one of the clipped reads at this
   site, so a pair is not reported twice.

2. Mates. The whole BAM is then read once (slow for large BAMs). The first
   forward-strand record with a kept left-clip read name goes to LEFT_MATE,
   the first reverse-strand record with a kept right-clip read name goes to
   RIGHT_MATE. Mates that map inside the inserted IAP, elsewhere in the
   genome (e.g. at the source copy) or not at all are all found this way.

Output format
-------------
A non-standard text file (not valid FASTA), one block per site, in the order
of the audit table:

    @chr15:98805467-98805473 (-)
    >RIGHT_CLIP
    TGTTATTCGACGCGTTCTCACGACC...
    TGTTATTCGACGCGTTCTCACGACCGGCCAGGAAGAACAC...
    >LEFT_CLIP
    >                    TATTTGTTTACCACTTAGAACACAGG
    >CCATGGCCGAGCTGACGTT...TATTTGTTTACCACTTAGAACACAGG
    >LEFT_MATE
    ...
    >RIGHT_MATE
    ...

- Header line "@seqname:start-end (strand)": start = first_bp from the audit
  table, end = start + 6 (assumed TSD, see above), strand = orientation of
  the inserted element.
- ">RIGHT_CLIP": one clipped sequence per line, shortest first. The first
  base of each line is the first base after the junction, so the lines are
  aligned at their left end.
- ">LEFT_CLIP": one clipped sequence per line, shortest first. Each line
  starts with ">" and is left-padded with spaces so all lines end in the
  same column, i.e. they are aligned at the junction (their right end).
- ">LEFT_MATE" / ">RIGHT_MATE": full mate sequences, one per line, may be
  empty.
All sequences are on the reference forward strand, as stored in the BAM
(not the orientation in which they were sequenced). There is no read name,
UMI or quality information, and duplicate molecules appear as repeated
lines.

Note: a site without any left-clipped read makes the script fail
(max() of an empty sequence) when writing the output.
"""

import argparse
import gzip
import os

import pandas as pd
from pysam import AlignmentFile, CSOFT_CLIP

parser = argparse.ArgumentParser(
    prog="collect_evidence.py",
    description=(
        "Gather evidence sequences and paired-mate reads flanking each "
        "candidate insertion breakpoint, for manual verification."
    ),
)
parser.add_argument(
    "--audit-xlsx",
    required=True,
    help="Path to isa.audit.xlsx (from collect_insertions.R)",
)
parser.add_argument("--bam", required=True, help="Path to the bam file for this sample")
parser.add_argument("--output", required=True, help="Path to output fa.gz file")
parser.add_argument(
    "--sample",
    required=True,
    help="Sample identifier to filter the audit table's 'mouse' column by",
)
args = parser.parse_args()

PATH_TO_AUDIT = args.audit_xlsx
PATH_TO_BAM = args.bam
PATH_TO_OUTPUT = args.output
SAMPLE = args.sample

if not os.path.isfile(PATH_TO_BAM):
    print(f"ERROR: BAM file not found: {PATH_TO_BAM}")
    exit(1)

d = pd.read_excel(PATH_TO_AUDIT)

reads = {}


class Breakpoint:
    def __init__(self, seqname, pos1, pos2, strand):
        self.seqname = seqname
        self.pos1 = pos1
        self.pos2 = pos2
        self.strand = strand
        self.left_clip = []
        self.right_clip = []
        self.left_mates = set()
        self.right_mates = set()
        self.covered_fwd = set()
        self.covered_rev = set()
        self.left_mate_seq = []
        self.right_mate_seq = []

    def add_left_clip(self, seq, name):
        self.left_clip.append(seq)

    def add_right_clip(self, seq, name):
        self.right_clip.append(seq)

    def cleanup(self):
        # make sure alrady used mates are not used again.
        self.left_mates = self.left_mates.difference(self.covered_fwd)
        self.right_mates = self.right_mates.difference(self.covered_rev)
        self.covered_fwd = None  # remove
        self.covered_rev = None  # remove
        self.left_clip = sorted(self.left_clip, key=lambda x: len(x))
        self.right_clip = sorted(self.right_clip, key=lambda x: len(x))

    def __str__(self):
        return f"{self.seqname}:{self.pos1}-{self.pos2} ({self.strand})"


with AlignmentFile(PATH_TO_BAM, mode="rb") as af:
    bps = []
    for i, (x, strand) in d[d["mouse"] == SAMPLE][["name", "strand"]].iterrows():
        seqname, pos = x.split(":")
        pos = int(pos)
        print(f"Processing {seqname}:{pos}")
        bp = Breakpoint(seqname, pos, pos + 6, strand)
        print(str(bp))
        for r in af.fetch(seqname, pos - 1, pos + 7):
            p = 0
            lclip = 0
            rclip = 0
            for m, l in r.cigartuples:
                if m == CSOFT_CLIP:
                    if p == 0:
                        lclip = l
                    else:
                        rclip = l
                else:
                    rclip = 0
                p += l
            if lclip and rclip:
                continue
            if not lclip and not rclip:
                continue
            ### reads are always mapped to + strand
            ### if read1 is not reversed, read 2 is to the right, read1 ^ is_forward = 1 ^ 1 = 0
            #    r1 ---->   <----- r2
            ### if read1 is reversed, read 2 is to the left read1 ^ is_forward = 1 ^ 0 = 1
            #     r2 ----->   <------- r1
            ### if read2 is reversed, read 1 is to the left read1 ^ is_forward = 0 ^ 0 = 0
            #    r1 ---->   <----- r2
            ### if read2 is not reversed, read 1 is to the right read1 ^ is_forward = 0 ^ 1 = 1
            #     r2 ----->   <------- r1
            ### for left clipped sequences, we are only interested in mates to the right, thus read1 ^ is_reverse = 1
            if lclip:
                bp.add_left_clip(r.query_sequence[:lclip], r.query_name)
                if r.is_reverse:
                    bp.left_mates.add(r.query_name)
            if rclip:
                bp.add_right_clip(
                    r.query_sequence[(len(r.query_sequence) - rclip) :], r.query_name
                )
                if r.is_forward:
                    bp.right_mates.add(r.query_name)

            # make sure a read is not used twice
            if r.is_forward:
                bp.covered_fwd.add(r.query_name)
            else:
                bp.covered_rev.add(r.query_name)
        bp.cleanup()
        bps.append(bp)
        print(
            f"missing {len(bp.left_mates)} left-mates and {len(bp.right_mates)} right-mates"
        )
    print(f"finding mates")
    af.reset()
    for r in af:
        if r.is_reverse:
            for bp in bps:
                if r.query_name in bp.right_mates:
                    bp.right_mate_seq.append(r.query_sequence)
                    bp.right_mates.remove(r.query_name)
        else:
            for bp in bps:
                if r.query_name in bp.left_mates:
                    bp.left_mate_seq.append(r.query_sequence)
                    bp.left_mates.remove(r.query_name)
    with gzip.open(PATH_TO_OUTPUT, "wt") as opt:
        for bp in bps:
            opt.write(f"@{str(bp)}\n")
            opt.write(f">RIGHT_CLIP\n")
            for seq in bp.right_clip:
                opt.write(f"{seq}\n")
            opt.write(f">LEFT_CLIP\n")
            rmaxlen = max(len(seq) for seq in bp.left_clip)
            for seq in bp.left_clip:
                opt.write(f">{' '*(rmaxlen-len(seq)) + seq}\n")
            opt.write(f">LEFT_MATE\n")
            for seq in bp.left_mate_seq:
                opt.write(f"{seq}\n")
            opt.write(f">RIGHT_MATE\n")
            for seq in bp.right_mate_seq:
                opt.write(f"{seq}\n")
