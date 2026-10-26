#!/usr/bin/env Rscript

# Standalone, non-interactive replacement for 06_cluster_source_coverage.R
# (PARTRIDGE step 6). For every manually validated insertion, counts the
# unique molecules (barcodes) spanning the edges of its source IAP in the
# sample(s) carrying the insertion. The source copy is fixed (homozygous) in
# inbred strains and is captured as efficiently as the sequence-identical
# novel copy, so this count is the denominator of the insertion's VAF
# (calculate_vaf.R).
#
# Input:
#   --audit-xlsx  audit table from collect_insertions.R, after manual
#                 validation (README step 5), with two added columns:
#                   classification  TRUE for validated insertions
#                   source          source IAP as "seqname:start-end"; any
#                                   sequence name is accepted (e.g. "chr6"
#                                   or Ensembl-style "6"), but it must match
#                                   the BAM's naming
#                 The "mouse" column holds the sample (Expt.ID) carrying the
#                 insertion, or "multi" when the insertion was found in more
#                 than one sample.
#   BAM files     one coordinate-sorted, indexed BAM per sample, named
#                 <Expt.ID>.bam (the Expt.ID must match the audit table).
#
# Method, per source region and sample:
#   1. samtools view all reads overlapping the source region.
#   2. Give each read a barcode: the UMI when the read name carries one as
#      its last ":"-separated field (9 nt, e.g. "...:1096:TCCACATCA");
#      otherwise a barcode synthesised per read pair from the first 4 bases
#      of the forward-oriented mate and the last 5 bases (reverse
#      complemented) of the reverse-oriented mate, i.e. from the two
#      fragment ends.
#   3. Keep reads with MAPQ >= 40 (after building the barcodes, so a pair's
#      synthesised barcode still uses both mates).
#   4. Take a 51 bp window centred on each outer edge of the source region
#      (start - 25 .. start + 25 and end - 25 .. end + 25) and keep reads that
#      span a whole window, i.e. that cross the junction between the IAP and
#      its flanking genome with at least 25 bp on each side.
#   5. Coverage = number of unique barcodes among those reads, pooled over
#      both windows.
#
# Output (--output-xlsx): the audit table with a source_coverage column:
#   - validated single-sample rows: coverage in that sample
#   - validated "multi" rows: Inf (xlsx cannot store Inf, so the cell reads
#     back as NA), and one extra row per sample (appended at the end, "mouse"
#     set to the sample) with the coverage in that sample
#   - non-validated rows: NA
# This is the layout calculate_vaf.R expects.
# Optional (--source-table-xlsx): coverage of every validated source in every
# sample (columns source, Expt.ID, coverage), e.g. to compare capture
# efficiency between samples.
#
# Differences from 06_cluster_source_coverage.R:
#   - the window at the source's start is centred on the 5' edge. In the old
#     script, `virus.gr - 25` shrinks the range instead of shifting it, so its
#     left window was start + 25 .. start + 75, entirely inside the 5' LTR,
#     and never tested the 5' junction (the right window was correct). On
#     test_data this changes the source coverage from 287 to 191, and the
#     README example VAF from 6/287 = 2.1% to 6/191 = 3.1%.
#   - all paths are arguments; sources are taken from the audit table instead
#     of a hard-coded list
#   - BAMs are found by "<Expt.ID>.bam" (the old pattern "^\w+\.bam$" missed
#     sample IDs containing "." or "-", e.g. "10.3C")
#   - source regions without a "chr" prefix are accepted
#   - a UMI is detected by a ":<9 nt>" suffix of the read name instead of a
#     read-name length >= 46, the same rule collect_insertions.R uses in
#     spirit (it checks the last character)
#   - sources wider than --max-source-width get NA (with a warning) instead
#     of 0, and a region without reads gets 0 instead of an error
#   - errors (missing BAM, failing samtools, malformed source) stop the run
#     instead of being turned into NA by try()
#   - an empty (NA) classification counts as not validated instead of
#     raising an error
#   - coverage of each source/sample pair is computed once
#
# Requires samtools on PATH (or --samtools).

suppressMessages({
  library(GenomicRanges)
  library(dplyr)
  library(openxlsx)
  library(Biostrings)
})

# Reads must have MAPQ >= MIN_MAPQ
MIN_MAPQ <- 40
# Reads must extend at least FLANK bp on both sides of a source edge
FLANK <- 25
# UMI appended to the read name as its last ":"-separated field
UMI_REGEX <- ":[ACGTN]{9}$"
UMI_LENGTH <- 9

usage <- paste(
  "Usage: Rscript cluster_source_coverage.R [options]",
  "",
  "Required, exactly one of:",
  "  --bam-dir DIR           Directory with <Expt.ID>.bam files",
  "  --bam FILE [FILE ...]   One or more explicit BAM paths (<Expt.ID>.bam)",
  "",
  "Required:",
  "  --audit-xlsx FILE   Validated audit table (collect_insertions.R output",
  "                      plus classification and source columns)",
  "  --output-xlsx FILE  Output path for the audit table with source_coverage",
  "",
  "Optional:",
  "  --source-table-xlsx FILE  Also write the coverage of every validated",
  "                            source in every sample",
  "  --samtools PATH           samtools executable (default: samtools)",
  "  --max-source-width N      Sources wider than N bp get NA coverage",
  "                            (default: 10000)",
  "",
  "  -h, --help          Show this help message",
  sep = "\n"
)

parse_args <- function(argv) {
  opts <- list(
    bam = character(0),
    samtools = "samtools",
    max_source_width = "10000"
  )
  i <- 1
  require_value <- function(key) {
    if (i == length(argv) || startsWith(argv[i + 1], "--")) {
      stop(sprintf("Missing value for argument %s", key), call. = FALSE)
    }
  }
  single <- c(
    "--bam-dir" = "bam_dir",
    "--audit-xlsx" = "audit_xlsx",
    "--output-xlsx" = "output_xlsx",
    "--source-table-xlsx" = "source_table_xlsx",
    "--samtools" = "samtools",
    "--max-source-width" = "max_source_width"
  )
  while (i <= length(argv)) {
    key <- argv[i]
    if (key %in% c("-h", "--help")) {
      cat(usage, "\n")
      quit(save = "no", status = 0)
    } else if (key == "--bam") {
      j <- i + 1
      vals <- character(0)
      while (j <= length(argv) && !startsWith(argv[j], "--")) {
        vals <- c(vals, argv[j])
        j <- j + 1
      }
      if (length(vals) == 0) {
        stop("--bam requires at least one file path", call. = FALSE)
      }
      opts$bam <- vals
      i <- j
    } else if (key %in% names(single)) {
      require_value(key)
      opts[[single[[key]]]] <- argv[i + 1]
      i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s\n\n%s", key, usage), call. = FALSE)
    }
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c("audit_xlsx", "output_xlsx")
missing <- required[
  !vapply(required, function(x) !is.null(opts[[x]]), logical(1))
]
if (length(missing) > 0) {
  stop(
    sprintf(
      "Missing required argument(s): %s\n\n%s",
      paste0("--", gsub("_", "-", missing), collapse = ", "),
      usage
    ),
    call. = FALSE
  )
}
if (is.null(opts$bam_dir) && length(opts$bam) == 0) {
  stop(sprintf("Provide either --bam-dir or --bam\n\n%s", usage), call. = FALSE)
}
if (!is.null(opts$bam_dir) && length(opts$bam) > 0) {
  stop("Provide either --bam-dir or --bam, not both", call. = FALSE)
}
max_source_width <- suppressWarnings(as.numeric(opts$max_source_width))
if (is.na(max_source_width) || max_source_width <= 0) {
  stop("--max-source-width must be a positive number", call. = FALSE)
}
if (!nzchar(Sys.which(opts$samtools))) {
  stop(sprintf("samtools not found: %s", opts$samtools), call. = FALSE)
}

## Locate BAM files -----------------------------------------------------------
# The sample (Expt.ID) of each BAM is its file name without ".bam", the same
# convention as the "<Expt.ID>.isa.tsv" files of isa_hmmer2.py.

if (!is.null(opts$bam_dir)) {
  bam_files <- list.files(opts$bam_dir, pattern = "\\.bam$", full.names = TRUE)
  if (length(bam_files) == 0) {
    stop(sprintf("No *.bam files found in %s", opts$bam_dir), call. = FALSE)
  }
} else {
  bam_files <- opts$bam
  missing_files <- bam_files[!file.exists(bam_files)]
  if (length(missing_files) > 0) {
    stop(
      sprintf(
        "BAM file(s) not found: %s",
        paste(missing_files, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
names(bam_files) <- sub("\\.bam$", "", basename(bam_files))
dup_samples <- unique(names(bam_files)[duplicated(names(bam_files))])
if (length(dup_samples) > 0) {
  stop(
    sprintf(
      "Several BAM files for sample(s): %s",
      paste(dup_samples, collapse = ", ")
    ),
    call. = FALSE
  )
}
no_index <- bam_files[
  !file.exists(paste0(bam_files, ".bai")) &
    !file.exists(sub("\\.bam$", ".bai", bam_files)) &
    !file.exists(paste0(bam_files, ".csi"))
]
if (length(no_index) > 0) {
  stop(
    sprintf(
      "BAM file(s) without index (run samtools index): %s",
      paste(no_index, collapse = ", ")
    ),
    call. = FALSE
  )
}

## Read and check the audit table ---------------------------------------------

audit <- read.xlsx(opts$audit_xlsx)

missing_cols <- setdiff(
  c("name", "mouse", "classification", "source"),
  names(audit)
)
if (length(missing_cols) > 0) {
  stop(
    sprintf(
      paste(
        "Audit table %s lacks column(s) %s; add classification and source",
        "during manual validation (README step 5)"
      ),
      opts$audit_xlsx,
      paste(missing_cols, collapse = ", ")
    ),
    call. = FALSE
  )
}

# classification may be read as logical or as text ("TRUE"/"FALSE"); empty
# cells count as not validated
audit$mouse <- trimws(as.character(audit$mouse))
audit$source <- trimws(as.character(audit$source))
validated <- as.logical(audit$classification) %in% TRUE

no_source <- validated & (is.na(audit$source) | audit$source == "")
if (any(no_source)) {
  warning(
    sprintf(
      "Validated insertion(s) without a source, coverage set to NA: %s",
      paste(audit$name[no_source], collapse = ", ")
    ),
    call. = FALSE
  )
}
todo <- validated & !no_source

missing_bams <- setdiff(
  audit$mouse[todo & audit$mouse != "multi"],
  names(bam_files)
)
if (length(missing_bams) > 0) {
  stop(
    sprintf(
      "No BAM file for sample(s) in the audit table: %s",
      paste(missing_bams, collapse = ", ")
    ),
    call. = FALSE
  )
}

## Source coverage --------------------------------------------------------------

# Parse "seqname:start-end"; any sequence name without ":" is accepted
parse_region <- function(region) {
  m <- regmatches(region, regexec("^([^:]+):(\\d+)-(\\d+)$", region))[[1]]
  if (length(m) == 0) {
    stop(
      sprintf(
        "Malformed source region '%s', expected seqname:start-end",
        region
      ),
      call. = FALSE
    )
  }
  GRanges(m[2], IRanges(start = as.numeric(m[3]), end = as.numeric(m[4])))
}

# Number of reference bases covered by each CIGAR (M, D, N, = and X)
cigar_ref_length <- function(cigar) {
  vapply(
    regmatches(cigar, gregexpr("[0-9]+[MDN=X]", cigar)),
    function(ops) sum(as.numeric(sub("[MDN=X]$", "", ops))),
    numeric(1)
  )
}

# Barcode per read: the UMI from the read name if present, otherwise a
# barcode synthesised from both fragment ends of the read pair
read_barcodes <- function(rname, flag, seq) {
  barcode <- rep(NA_character_, length(rname))

  has_umi <- grepl(UMI_REGEX, rname)
  barcode[has_umi] <- substr(
    rname[has_umi],
    nchar(rname[has_umi]) - UMI_LENGTH + 1,
    nchar(rname[has_umi])
  )

  if (any(!has_umi)) {
    synth <- tibble(
      rname = rname[!has_umi],
      seq = seq[!has_umi],
      first = bitwAnd(flag[!has_umi], 0x40) > 0,
      reverse = bitwAnd(flag[!has_umi], 0x10) > 0
    )
    # rev: read is the reverse-oriented end of the fragment (read 1 on the
    # reverse strand, or read 2 on the forward strand). Take the first 4
    # bases of the forward end and the last 5 bases (reverse complemented)
    # of the reverse end, i.e. the bases at the two fragment ends.
    # Records without a stored sequence ("*") contribute nothing.
    synth$seq[synth$seq == "*"] <- ""
    synth$rev <- !xor(synth$first, synth$reverse)
    synth$end_seq <- substr(synth$seq, 1, 4)
    synth$end_seq[synth$rev] <- as.character(reverseComplement(DNAStringSet(
      substr(
        synth$seq[synth$rev],
        nchar(synth$seq[synth$rev]) - 4,
        nchar(synth$seq[synth$rev])
      )
    )))
    # Concatenate read 1's and read 2's end sequences per read name
    synth <- synth |>
      mutate(
        barcode_fwd = ifelse(first, end_seq, ""),
        barcode_rev = ifelse(first, "", end_seq)
      ) |>
      group_by(rname) |>
      summarise(barcode = paste0(barcode_fwd, barcode_rev, collapse = ""))
    barcode[!has_umi] <- synth$barcode[match(rname[!has_umi], synth$rname)]
  }
  barcode
}

get_source_coverage <- function(region, bam) {
  region_gr <- parse_region(region)
  if (width(region_gr) > max_source_width) {
    warning(
      sprintf(
        "Source %s is wider than %d bp, coverage set to NA",
        region,
        max_source_width
      ),
      call. = FALSE
    )
    return(NA_real_)
  }

  sam <- suppressWarnings(system2(
    opts$samtools,
    c("view", shQuote(bam), shQuote(region)),
    stdout = TRUE
  ))
  status <- attr(sam, "status")
  if (!is.null(status) && status != 0) {
    stop(
      sprintf("samtools view failed for %s in %s", region, bam),
      call. = FALSE
    )
  }
  if (length(sam) == 0) {
    return(0)
  }

  fields <- strsplit(sam, "\t", fixed = TRUE)
  field <- function(k) vapply(fields, `[[`, character(1), k)
  rname <- field(1)
  flag <- as.integer(field(2))
  pos <- as.numeric(field(4))
  mapq <- as.numeric(field(5))
  cigar <- field(6)
  seq <- field(10)

  reads <- GRanges(
    field(3),
    IRanges(start = pos, end = pos + cigar_ref_length(cigar)),
    barcode = read_barcodes(rname, flag, seq)
  )
  reads <- reads[mapq >= MIN_MAPQ]

  # 51 bp windows centred on the two outer edges of the source
  edges <- c(
    GRanges(
      seqnames(region_gr),
      IRanges(start(region_gr) - FLANK, start(region_gr) + FLANK)
    ),
    GRanges(
      seqnames(region_gr),
      IRanges(end(region_gr) - FLANK, end(region_gr) + FLANK)
    )
  )

  # Reads that contain a whole window (window "within" read)
  hits <- findOverlaps(edges, reads, type = "within")
  length(unique(reads$barcode[subjectHits(hits)]))
}

# Source/sample pairs to compute: the sample of each single-sample insertion,
# every sample for "multi" insertions, and every sample for the optional
# source table
pairs <- bind_rows(
  tibble(
    source = audit$source[todo & audit$mouse != "multi"],
    Expt.ID = audit$mouse[todo & audit$mouse != "multi"]
  ),
  expand.grid(
    source = unique(audit$source[todo & audit$mouse == "multi"]),
    Expt.ID = names(bam_files),
    stringsAsFactors = FALSE
  ),
  if (!is.null(opts$source_table_xlsx)) {
    expand.grid(
      source = unique(audit$source[todo]),
      Expt.ID = names(bam_files),
      stringsAsFactors = FALSE
    )
  }
) |>
  distinct()

pairs$coverage <- NA_real_
for (k in seq_len(nrow(pairs))) {
  message(sprintf("Source %s in %s", pairs$source[k], pairs$Expt.ID[k]))
  pairs$coverage[k] <- get_source_coverage(
    pairs$source[k],
    bam_files[[pairs$Expt.ID[k]]]
  )
}

coverage_of <- function(source, sample) {
  pairs$coverage[match(
    paste(source, sample),
    paste(pairs$source, pairs$Expt.ID)
  )]
}

## Write outputs ----------------------------------------------------------------
# Same layout as 06_cluster_source_coverage.R: "multi" rows get Inf and are
# expanded into one appended row per sample.

audit$source_coverage <- NA_real_
single <- todo & audit$mouse != "multi"
audit$source_coverage[single] <- coverage_of(
  audit$source[single],
  audit$mouse[single]
)

multi <- which(todo & audit$mouse == "multi")
audit$source_coverage[multi] <- Inf
expanded <- bind_rows(lapply(multi, function(i) {
  rows <- audit[rep(i, length(bam_files)), , drop = FALSE]
  rows$mouse <- names(bam_files)
  rows$source_coverage <- coverage_of(rows$source, rows$mouse)
  rows
}))
out <- bind_rows(audit, expanded)
rownames(out) <- NULL

write.xlsx(out, opts$output_xlsx)
message(sprintf("Wrote %d row(s) to %s", nrow(out), opts$output_xlsx))

if (!is.null(opts$source_table_xlsx)) {
  source_table <- pairs |>
    filter(source %in% audit$source[todo]) |>
    arrange(source, Expt.ID)
  write.xlsx(as.data.frame(source_table), opts$source_table_xlsx)
  message(sprintf("Wrote source table to %s", opts$source_table_xlsx))
}

message("Done.")
