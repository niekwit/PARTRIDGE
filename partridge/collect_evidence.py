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
