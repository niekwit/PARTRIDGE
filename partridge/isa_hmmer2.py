import argparse
import concurrent.futures

import pysam
import pyhmmer
import cigar
from collections.abc import Generator
import sys
from concurrent.futures import ProcessPoolExecutor, Future
import multiprocessing
import time
import os
import math

# CONVENTION
# All coordiantes are 0-based
# the start coordinate is inclusive, thus includes position
# the end coordinate is exclusive, thus does not include this position
# example:
# 0 1 2 3 4 5 6 7 -> length = 8
#   =======   -> start = 1, end = 5, length = 5-1=4
# invert: start = length-end, end = length-start => start = 3, end = 7

# Flags indicating a model end
MIN_CLIPPED_BASES = 24 #must have MORE than this number of clipped bases to be submitted.
SIDE_NONE = 2
SIDE_LEFT = 0
SIDE_RIGHT = 1
MAPQ_CUTOFF = 40 # do not look at alignments with MAPQ smaller than this

DEBUG = [
    ("chr1",116518567),
    ("chr6",8573744),
    ("chr6",83190171),
    ("chr13",93389886),
    ("chr14",99124836),
]
DEBUG = None
class ReportedBreakpoint:
    """
    Class to store the final output and to format the output text file
    """
    record = {
        'rname': 'NA', #name of the read, may include the barcode
        'barcode': 'NA', #barcode. If none are present, takes the first 9 bases of the read.
        'seqname': 'NA', #seqname of the alignment
        'pair_read': 'NA', #number of the read
        'strand': '*', #strand of the alignment
        'breakpoint': 'NA', #breakpoint, thus last or first base of the alignment outside the LTR.
        'ltr.end': 'NA', #end of the ltr, 5' or 3', F or R
        'ltr.len': 'NA', # length of the LTR, max the length of the HMM, thus 25bp
        'ltr.score': 'NA', #bit score of the LTR
        'aln.len': 'NA', #length of the alignment, excluding soft clipped bases (which are the LTR)
        'bp_delta': 'NA', #number of basepaires between the border of the LTR and the border of the alignment
        'mapq': 'NA', #mapping quality of the alignment
        'ltr.model.border': 'NA', #how many bases of the model end were not matched? Should be 0 for a exactly adjacent match.
    }

    def __init__(self, hmmer: 'hmmerHit'):

        if not hmmer: return
        self.record['seqname'] = hmmer.seqname
        self.record['pair_read'] = "1" if hmmer.is_read1 else "2"
        self.record['barcode'] = hmmer.barcode
        self.record['strand'] = hmmer.strand
        self.record['rname'] = hmmer.query_name
        self.record['breakpoint'] = hmmer.reference_breakpoint
        self.record['bp_delta'] = hmmer.breakpoint_delta
        self.record['mapq'] = hmmer.mapq
        self.record['ltr.len'] = hmmer.len
        self.record['ltr.score'] = hmmer.bits
        self.record['aln.len'] = hmmer.aln_len
        self.record['ltr.end'] = hmmer.end
        self.record['ltr.model.border'] = hmmer.model_border


    @staticmethod
    def tsvtitle():
        return "\t".join([
            "rname",
            'barcode',
            'pair.read',
            "seqname",
            "strand",
            "breakpoint",
            "ltr.end",
            "ltr.len",
            "ltr.score",
            "aln.len",
            "bp.delta",
            "mapq",
            'ltr.model.border'
        ])

    def tsvline(self):
        return "\t".join([
            self.record['rname'],
            self.record['barcode'],
            self.record['pair_read'],
            str(self.record['seqname']),
            str(self.record['strand']),
            str(self.record['breakpoint']),
            str(self.record['ltr.end']),
            str(self.record['ltr.len']),
            str(self.record['ltr.score']),
            str(self.record['aln.len']),
            str(self.record['bp_delta']),
            str(self.record['mapq']),
            str(self.record['ltr.model.border']),
        ])


class Alignment:
    """
    Class storing a single bwa-mem2 alignment
    """

    def __init__(self, record: pysam.AlignedSegment):
        self.mapq = record.mapq
        self.seqname = record.reference_name
        self.query_name = record.query_name
        self.ref_start = record.reference_start
        self.ref_end = record.reference_end
        self.query_length = record.query_length
        self.forward = not record.is_reverse
        self.is_read1 = record.is_read1
        self._cigarstring = record.cigarstring
        self.query_sequence = record.query_sequence
        self._cigarobj = None
        self.rev = record.is_read1 ^ record.is_forward

        pre = 0
        post = 0
        self.aligned_len = 0
        left = True
        for c in self.cigar.cigar:
            if c.operator == "S":
                if left:
                    pre += c.length
                else:
                    post += c.length
            elif left and c.operator == "M":
                self.aligned_len += c.length
                left = False
        self.clipped = max(pre, post)
        if pre > post:
            self.side = SIDE_LEFT
        elif post > pre:
            self.side = SIDE_RIGHT
        else:
            self.side = SIDE_NONE

    def _makeCigarObj(self) -> None:
        self._cigarobj = cigar.Cigar(self._cigarstring)
        self._record_start, self._record_end = self._cigarobj.startend()

    def reference_position(self, query_position) -> int:
        """
        Get Reference position for query_position
        :param query_position: 0-based query position you would like to figure out the ref position for
        Important: as per SAM convention, the cigar-string is always relative to the + strand, even if the match is
        actually on the - strand.
        :return: integer (reference_position) or None if the ref position is not defined. 0-based
        """
        q = 0
        ref = self.ref_start
        # tuples storing the first and the last ref position of this match, query_pos then ref_pos
        firsttuple = None
        lasttuple = None
        for c in self.cigar.cigar:
            if c.operator in cigar.CIGAR_CONSUMES_BOTH:
                lasttuple = q, ref
                if firsttuple is None:
                    firsttuple = lasttuple
                if query_position >= q:
                    if query_position < q + c.length:
                        return ref + query_position - q
                q += c.length
                ref += c.length
            else:
                if c.operator in cigar.CIGAR_CONSUMES_QUERY:
                    q += c.length
                if c.operator in cigar.CIGAR_CONSUMES_REFERENCE:
                    ref += c.length
        if firsttuple is not None and query_position < firsttuple[0]:
            return firsttuple[1] + query_position - firsttuple[0]
        if lasttuple is not None and query_position > lasttuple[0]:
            return lasttuple[1] + query_position - lasttuple[0]
        return None

    @property
    def record_start(self) -> int:
        if self._cigarobj is None:
            self._makeCigarObj()
        return self._record_start

    @property
    def record_end(self) -> int:
        if self._cigarobj is None:
            self._makeCigarObj()
        return self._record_end

    @property
    def cigar(self) -> cigar.Cigar:
        if self._cigarobj is None:
            self._makeCigarObj()
        return self._cigarobj

    @property
    def strand(self):
        if self.forward:
            return "+"
        else:
            return "-"

    @staticmethod
    def revcomp(seq: str):
        seq = list(seq.upper())
        seq.reverse()
        for i in range(len(seq)):
            if seq[i] == "A":
                seq[i] = "T"
            elif seq[i] == "T":
                seq[i] = "A"
            elif seq[i] == "G":
                seq[i] = "C"
            elif seq[i] == "C":
                seq[i] = "G"
        return "".join(seq)

    def positionstr(self):
        return f"{self.seqname}:{self.ref_start}-{self.ref_end}{self.strand}"

    def __str__(self) -> str:
        if self.reference_breakpoint is not None:
            return f"{self.positionstr()} [query=({self.record_start}-{self.record_end}) BP={self.reference_breakpoint}, ∂={self.breakpoint_delta}]"
        else:
            return self.positionstr()

    def __len__(self) -> int:
        return self.aligned_len


class hmmerHit:
    """
    Stores a single line of HMMER result.
    """

    def __init__(self, alignment: Alignment, hit: pyhmmer.plan7.Domain):
        """
        HMMER Hits are 1-based
        thus "from" values and "spot" values have to be converted to 0-based by subtracting 1.
        Now, the problem is that env_from and env_to can be inverted if the model match is inverted.
        :param hit:
        """
        self.seqname = alignment.seqname
        self.strand = alignment.strand
        self.mapq = alignment.mapq
        self.query_name = alignment.query_name
        self.aln_len = len(alignment)
        self.is_read1 = alignment.is_read1
        self.cigar = str(alignment.cigar)
        self.accession = None
        self.accession = hit.alignment.hmm_accession
        self.model = str(hit.alignment.hmm_name)
        self.model_prime = self.model[len(self.model)-3]
        assert self.model_prime in ("3","5"), f"model {self.model} has invalid prime {self.model_prime}"
        self.model_end = self.model[len(self.model) - 2]
        assert self.model_end in ("R", "F"), f"model {self.model} has invalid directionality {self.model_prime}"
        self.end = self.model_prime + self.model_end
        if ( self.model_prime == "5" and self.model_end == "F" ) or ( self.model_prime == "3" and self.model_end == "R" ):
                self.side = SIDE_RIGHT
                self.model_border = hit.alignment.hmm_from -1
                self.query_breakpoint = hit.alignment.target_from-1
                self.breakpoint_delta = self.query_breakpoint - alignment.query_length + alignment.clipped
        if ( self.model_prime == "5" and self.model_end == "R" ) or ( self.model_prime == "3" and self.model_end == "F"):
                self.side = SIDE_LEFT
                self.model_border = hit.alignment.hmm_length - hit.alignment.hmm_to
                self.query_breakpoint = hit.alignment.target_to
                self.breakpoint_delta =  alignment.clipped - self.query_breakpoint
        self.reference_breakpoint = alignment.reference_position(self.query_breakpoint)
        self.bits = hit.score
        self.len = hit.alignment.target_length
        self.passes_filter = True
        self.alignment_side = alignment.side
        self.clipped = alignment.clipped
        if self.alignment_side != self.side:
            self.passes_filter = False
        if abs(self.breakpoint_delta) > 5:
            self.passes_filter = False
        if self.model_border != 0:
            self.passes_filter = False

        # breakpoint_delta: Distance between LTR and Alignment. Set by Record:borderAlignments
        # positive integer: Gap between both alignments in bp
        # negative integer: Overlap between both alignments in bp
        # 0 = both alignments are immediately adjacent
        # None= not set.
        #
        # Examples :
        # A = Alignment
        # L = LTR identified by HMMER
        #
        # AAAAAAAAAAAAA-----------
        # ----------LLLLLLLLLLLLLL
        # breakpoint_delta = -3
        #
        # AAAAAAAAAAAA---------------
        # --------------LLLLLLLLLLLLL
        # breakpoint_delta = +2
        #
        # AAAAAAAAAAA----------------
        # -----------LLLLLLLLLLLLLLLL
        # breakpoint_delta = 0
        # Example
        # Query:     ----AAAAAAAA----------
        #            ----------LLLLLLLLL---
        # a =4, l=3, breakpoint_delta = -2, A = 6, L=7
        # figure out padded ends

        self.extract_barcode(alignment)


    def extract_barcode(self, alignment: Alignment):
        if False and len(self.query_name)>45:
            self.barcode = alignment.query_name[-9:]
        else:
            if alignment.rev:
                self.barcode = Alignment.revcomp(alignment.query_sequence[-9:])
            else:
                self.barcode = alignment.query_sequence[0:9]

    def __len__(self) -> int:
        """
        Get length of match in MODEL
        :return:
        """
        if self.model is None:
            return 0
        return self.len

    def __float__(self) -> float:
        """
        Get ratio of length of match in ENV in relation to length of match in MODEL
        :return:
        """
        if len(self) == 0: return 0
        return self.bits / len(self)

    def endstring(self):
        return self.model_prime + ("L" if self.model_side == SIDE_LEFT else "R")

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.seqname}:{self.reference_breakpoint} {self.end} bp={self.query_breakpoint} {self.cigar} mb={self.model_border} bpd={self.breakpoint_delta} model_side={'L' if self.side==SIDE_LEFT else 'R'} alignment_side={'L' if self.alignment_side==SIDE_LEFT else ('R' if self.alignment_side==SIDE_RIGHT else 'ø')} alignment_clipped={self.clipped} pass={self.passes_filter}"

class RecordChunk:

    def __init__(self):
        self.records = []

    def add(self, record: Alignment = None) -> None:
        """
        THIS FUNCTION IS NOT THREADSAFE!!
        :param record:
        :return:
        """
        record.id = int(len(self.records)).to_bytes(2, byteorder=sys.byteorder, signed=False)
        self.records.append(record)

    def __getitem__(self, item: bytes) -> Alignment:
        return self.records[int.from_bytes(item, byteorder=sys.byteorder, signed=False)]

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Generator[Alignment]:
        for r in self.records:
            return r

    def digitalSequenceIterator(self, alphabet) -> Generator[pyhmmer.easel.DigitalSequence]:
        for record in self.records:
            yield pyhmmer.easel.DigitalSequence(alphabet, sequence=alphabet.encode(record.query_sequence),
                                                accession=record.id, name=record.query_name.encode('ascii'))

    def digitalSequenceBlock(self) -> pyhmmer.easel.DigitalSequenceBlock:
        dna = pyhmmer.easel.Alphabet.dna()
        return pyhmmer.easel.DigitalSequenceBlock(
            dna,
            iterable=self.digitalSequenceIterator(dna)
        )


class nHmmer:
    """
    Class to perform nhmmer search
    handles a list of Records (=chunk) for optimization reasons
    Has problems if the length of this list is > 50, I don't know why, this is a pyhmmer problem.
    """
    dna = pyhmmer.easel.Alphabet.dna()

    def __init__(self, records: RecordChunk, hmm: pyhmmer.plan7.OptimizedProfileBlock, threads: int):
        self.seq = records.digitalSequenceBlock()
        self.threads = threads
        self._records = records
        self.run(hmm)

    def run(self, hmm: pyhmmer.plan7.OptimizedProfileBlock):
        self._results = pyhmmer.hmmer.hmmsearch(queries=hmm, sequences=self.seq, cpus=self.threads, E=10e-3)

    def prepareResults(self) -> Generator[hmmerHit]:
        results = {}
        for topHitsForHmm in self._results:
            for hit in topHitsForHmm.reported:
                # exclude non-edge results
                if hit.accession is None: continue
                domain = hit.best_domain
                if hit.accession in results.keys():
                    # keep only best result
                    if results[hit.accession].c_evalue > domain.c_evalue:
                        results[hit.accession] = domain
                else:
                    results[hit.accession] = domain

        for key in results.keys():
            x = hmmerHit(self._records[key], results[key])
            if x.passes_filter:
                yield (x)
            if x.reference_breakpoint is not None and DEBUG is not None:
                for seqname, pos in DEBUG:
                    if x.seqname == seqname and x.reference_breakpoint > pos-250 and x.reference_breakpoint < pos+250:
                        print(x)


result_line = 0
result_block = 0


def initPPE(hmmfile, threads_per_process):
    global hmm, subthreads
    subthreads = threads_per_process
    hmm = pyhmmer.plan7.OptimizedProfileBlock(alphabet=pyhmmer.easel.Alphabet.dna(),
                                              iterable=pyhmmer.plan7.HMMPressedFile(hmmfile))

def workPPE(chunk: RecordChunk):
    global subthreads
    x = list(nHmmer(chunk, hmm, threads=subthreads).prepareResults())
    return (x)


def writeResult(future: concurrent.futures.Future):
    global result_line, result_block
    result_block += 1
    if result_block % 1000 == 0:
        print(f"> writing result line {result_line}")
    for r in future.result():
        result_line += 1
        output.write(ReportedBreakpoint(r).tsvline())
        output.write("\n")


class CommandLineManager:
    title = """Integration Site Analysis - HMMER Tool 2.0
--------------------------------------
by Jeremy Deuel <jeremy.deuel@usz.ch>, June 2024

"""

    def __init__(self):

        print(CommandLineManager.title)
        parser = argparse.ArgumentParser(prog="isa_hmmer2.py")
        parser.add_argument("--bam", required=True, help="Path to bam file. Mandatory")
        parser.add_argument("--output", required=True,
                             help="Path to output (tsv) file. Mandatory. If the file already exists, it will be overwritten")
        parser.add_argument("--threads", type=int, default=None,
                             help="Number of threads to use. Optional, will be set to the number of CPUs available if not set")
        parser.add_argument("--hmm", default="resources/iap.hmm",
                             help="Path to the hmm file. Optional, will be set to resources/iap.hmm if not set.")
        parser.add_argument("--chunksize", type=int, default=100,
                             help="Number of records to be merged together into a chunk for multiprocessing. Optional, will be set to 100 if not set.")
        args = parser.parse_args()

        self.bam = args.bam
        self.output = args.output
        self.threads = args.threads if args.threads is not None else multiprocessing.cpu_count()
        self.cores = self.threads
        # divide threads by number of subthreads
        self.subthreads = 1
        self.threads = max(1, math.floor(self.threads / self.subthreads) - 1)
        self.hmmfile = args.hmm
        self.chunksize = args.chunksize
        self.chunksleep = 10  # number of seconds to sleep if chunkbuffer is full
        self.chunkbuffer = math.ceil(
            500 / self.chunksize * self.chunksleep * self.subthreads * self.threads)  # max number of unprocessed chunks allowed

        print(f"""Parameters:
BAM-File:       {self.bam}
Output-File:    {self.output}
Processes       {self.threads}, will be used: {self.threads * self.subthreads + 1}, total cores requested {self.cores}, available {multiprocessing.cpu_count()}
HMM-File:       {self.hmmfile}
Chunk-Size:     {self.chunksize}
Chunk-Buffer:   {self.chunkbuffer} (fixed value)
Sleep-Time:     wait for {self.chunksleep}s before adding more input if chunk-buffer is full. (fixed value)
Subprocesses    {self.subthreads} (fixed value, number of threads per subprocess)
""")
        # if self.threads > multiprocessing.cpu_count():
        #    print("ERROR: Requested more threads than available, aborting.")
        #    exit(1)
        if self.threads < 1:
            print("ERROR: Requested less than one thread, aborting.")
            exit(1)
        if not os.path.isfile(self.bam):
            print("ERROR: BAM File does not exist, aborting.")
            exit(1)
        if not os.path.getsize(self.bam):
            print("ERROR: BAM File is empty, aborting.")
            exit(1)
        if os.path.isfile(self.output):
            print("WARNING: Output file exists, will be overwritten.")
        # if os.path.isfile(self.output) and not not os.access(self.output, os.W_OK):
        #    print("ERROR: Output file is not writable, aborting.")
        #    exit(1)
        if not os.path.isfile(self.hmmfile):
            print("ERROR: HMM-File does not exist, aborting.")
            exit(1)
        if self.chunksize < 1:
            print("ERROR: Chunk-Size is smaller than one, aborting.")
            exit(1)
        if self.chunksize > 999999:
            print("ERROR: Chunk-Size is unreasonably high, aborting.")


if __name__ == "__main__":
    t = time.time()
    clm = CommandLineManager()
    samfile = pysam.AlignmentFile(clm.bam, "rb")
    output = open(clm.output, "w")
    output.write(ReportedBreakpoint.tsvtitle())
    output.write("\n")
    ppe = ProcessPoolExecutor(max_workers=clm.threads,
                              mp_context=multiprocessing.get_context('spawn'),
                              initializer=initPPE, initargs=(clm.hmmfile, clm.subthreads))

    chunk = RecordChunk()
    record1 = None
    record2 = None
    record_name = None
    output_container = []
    chunk_number = 0
    record_number = 0
    for r in samfile:
        record_number += 1
        if r.is_unmapped: continue
        if r.is_duplicate: continue
        if r.is_qcfail: continue
        record = Alignment(r)
        if record.clipped > MIN_CLIPPED_BASES:
            chunk.add(record)
            if len(chunk) > clm.chunksize:
                while (chunk_number - result_block > clm.chunkbuffer):
                    print(f"   queue full, waiting for {clm.chunksleep}s before submitting further input")
                    time.sleep(clm.chunksleep)
                chunk_number += 1
                if chunk_number % 1000 == 0:
                    print(f"... submitted bam entry {record_number} for processing")
                ppe.submit(workPPE, chunk).add_done_callback(writeResult)
                chunk = RecordChunk()
    if len(chunk):
        chunk_number += 1
        ppe.submit(workPPE, chunk).add_done_callback(writeResult)
    print(f"...completed submission of whole file, total {record_number} records in {chunk_number} chunks.")
    ppe.shutdown()
    print(f"...completed writing of output file, total of {result_line} lines in {result_block} chunks.")
    output.close()
    samfile.close()
    print(f"done in {time.time() - t} seconds.")
