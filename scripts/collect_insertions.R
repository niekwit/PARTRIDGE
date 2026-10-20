#!/usr/bin/env Rscript

# Standalone, non-interactive replacement for 03_collect_insertions.Rmd
# (PARTRIDGE step 3). Aggregates the tsv files produced by isa_hmmer2.py into
# a set of candidate novel insertion sites, joins them against a metadata
# table, and writes out summary plots and Excel audit files.
#
# Input: one <Expt.ID>.isa.tsv per sample (isa_hmmer2.py, step 2). Each row is
# one read whose soft-clipped part matched an IAP LTR end (HMMER), i.e. a read
# spanning the junction between genome and an inserted LTR. Columns used here:
#   rname       read name; carries the UMI as its last ":" field when the
#               FASTQs were demultiplexed with UMIs (e.g. "...:1096:TCCACATCA")
#   barcode     9-mer taken from the read sequence by isa_hmmer2.py, used as a
#               substitute UMI when the read name has none
#   seqname     chromosome of the read's genomic alignment
#   breakpoint  genomic coordinate where the alignment meets the LTR
#   ltr.end     two characters: which LTR terminus matched ("5" or "3") and
#               the orientation of that match on the read ("F" or "R"),
#               which reflects the orientation of the inserted element
#   ltr.score   HMMER bit score of the LTR match
#   bp.delta    gap (>0) or overlap (<0) in bp between the genomic alignment
#               and the LTR match on the read; 0 = exactly adjacent
#   mapq        mapping quality of the genomic alignment
#
# Overview:
#   1. Per sample (novel_sites()): collapse reads to unique molecules per
#      junction, cluster nearby junctions, and keep clusters that look like
#      both ends of one inserted element.
#   2. Across samples: one row per site and sample with read support, and the
#      number of samples the site occurs in (somatic vs shared/germline).
#   3. Join sample metadata, plot, and write the insertion table (input for
#      calculate_vaf.R) and the audit table (input for manual validation,
#      collect_evidence.py and 06_cluster_source_coverage.R).

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
})

usage <- paste(
  "Usage: Rscript collect_insertions.R [options]",
  "",
  "Required, exactly one of:",
  "  --tsv-dir DIR           Directory to scan for *.isa.tsv files",
  "  --tsv FILE [FILE ...]   One or more explicit tsv file paths",
  "",
  "Required:",
  "  --metadata FILE     Path to metadata xlsx (must contain an Expt.ID column)",
  "  --somatic-pdf FILE  Output path for the somatic (single-mouse) plot",
  "  --germline-pdf FILE Output path for the germline (multi-mouse) plot",
  "                      (only written if germline sites are found)",
  "  --output-xlsx FILE  Output path for the aggregated insertions table",
  "  --audit-xlsx FILE   Output path for the manual-review audit table",
  "",
  "Optional:",
  "  --cache-rds FILE    Path to the aggregated-breakpoints rds cache. If it",
  "                      already exists, it is loaded instead of re-parsing",
  "                      the tsv files; otherwise it is created. If not given,",
  "                      no cache is read or written.",
  "  --facet-column NAME Metadata column to facet the plots by (e.g. a",
  "                      genotype column). If not given, plots are not",
  "                      faceted.",
  "",
  "  -h, --help          Show this help message",
  sep = "\n"
)

parse_args <- function(argv) {
  opts <- list(tsv = character(0))
  i <- 1
  require_value <- function(key) {
    if (i == length(argv) || startsWith(argv[i + 1], "--")) {
      stop(sprintf("Missing value for argument %s", key), call. = FALSE)
    }
  }
  while (i <= length(argv)) {
    key <- argv[i]
    if (key %in% c("-h", "--help")) {
      cat(usage, "\n")
      quit(save = "no", status = 0)
    } else if (key == "--tsv-dir") {
      require_value(key)
      opts$tsv_dir <- argv[i + 1]
      i <- i + 2
    } else if (key == "--tsv") {
      j <- i + 1
      vals <- character(0)
      while (j <= length(argv) && !startsWith(argv[j], "--")) {
        vals <- c(vals, argv[j])
        j <- j + 1
      }
      if (length(vals) == 0) {
        stop("--tsv requires at least one file path", call. = FALSE)
      }
      opts$tsv <- vals
      i <- j
    } else if (key == "--metadata") {
      require_value(key)
      opts$metadata <- argv[i + 1]
      i <- i + 2
    } else if (key == "--cache-rds") {
      require_value(key)
      opts$cache_rds <- argv[i + 1]
      i <- i + 2
    } else if (key == "--somatic-pdf") {
      require_value(key)
      opts$somatic_pdf <- argv[i + 1]
      i <- i + 2
    } else if (key == "--germline-pdf") {
      require_value(key)
      opts$germline_pdf <- argv[i + 1]
      i <- i + 2
    } else if (key == "--output-xlsx") {
      require_value(key)
      opts$output_xlsx <- argv[i + 1]
      i <- i + 2
    } else if (key == "--audit-xlsx") {
      require_value(key)
      opts$audit_xlsx <- argv[i + 1]
      i <- i + 2
    } else if (key == "--facet-column") {
      require_value(key)
      opts$facet_column <- argv[i + 1]
      i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s\n\n%s", key, usage), call. = FALSE)
    }
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c(
  "metadata",
  "somatic_pdf",
  "germline_pdf",
  "output_xlsx",
  "audit_xlsx"
)
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
if (is.null(opts$tsv_dir) && length(opts$tsv) == 0) {
  stop(sprintf("Provide either --tsv-dir or --tsv\n\n%s", usage), call. = FALSE)
}
if (!is.null(opts$tsv_dir) && length(opts$tsv) > 0) {
  stop("Provide either --tsv-dir or --tsv, not both", call. = FALSE)
}

## Identify novel sites from a single tsv file -------------------------------

novel_sites <- function(file) {
  d_file <- read.table(file, header = TRUE)

  # Replace the sequence-derived barcode with the real UMI where the read name
  # carries one. Standard Illumina read names end in a digit (y-coordinate);
  # names ending in A/C/G/T have a UMI appended, of which the last 9
  # characters are taken (9-nt xGen UMI). Both mates of a pair share the UMI,
  # so they collapse into one molecule below. UMIs ending in N fall back to
  # the sequence-derived barcode.
  nc_filter <- substr(d_file$rname, nchar(d_file$rname), nchar(d_file$rname))
  nc_filter <- nc_filter %in% c("A", "T", "C", "G")
  nc <- nchar(d_file$rname)[nc_filter]
  d_file$barcode[nc_filter] <- substr(d_file$rname[nc_filter], nc - 8, nc)

  # Count unique molecules per junction:
  # - keep reads whose LTR match is adjacent to or overlaps the genomic
  #   alignment (bp.delta <= 0); reads with unaligned bases between the
  #   genome and the LTR are dropped
  # - collapse reads sharing a barcode at the same junction (ltr.end +
  #   position) into one molecule (PCR duplicates), keeping its best LTR
  #   score and best MAPQ
  # - per junction: score = summed LTR score over molecules, n = number of
  #   unique molecules (barcodes), mapq = best MAPQ of any molecule
  # Note that barcodes are only compared within one junction, so the same
  # barcode at different positions counts as different molecules.
  d_agg <- d_file |>
    filter(bp.delta < 1) |>
    group_by(ltr.end, seqname, breakpoint, barcode) |>
    summarise(score = max(ltr.score), mapq = max(mapq)) |>
    group_by(ltr.end, seqname, breakpoint) |>
    summarise(score = sum(score), n = n(), mapq = max(mapq)) |>
    arrange(seqname, breakpoint)

  # Cluster junctions into candidate insertion sites. Rows are sorted by
  # chromosome and position; a row joins the previous row's group when it is
  # on the same chromosome and less than 150 bp downstream of it. This is
  # chaining (single linkage), so a group can span more than 150 bp in total.
  # The two junctions of one insertion (5' LTR end and 3' LTR end) are
  # separated only by the target-site duplication, so they fall in one group.
  group <- 1
  d_agg$group <- 0
  for (i in seq_len(nrow(d_agg) - 1)) {
    d_agg$group[i] <- group
    if (
      (d_agg$seqname[i] == d_agg$seqname[i + 1]) &
        (d_agg$breakpoint[i] + 150 > d_agg$breakpoint[i + 1])
    ) {
      d_agg$group[i + 1] <- group
    } else {
      group <- group + 1
    }
  }

  # Split ltr.end into the LTR terminus ("5"/"3") and the orientation ("F"/"R")
  d_agg$ltr.end.prime <- substr(d_agg$ltr.end, 1, 1)
  d_agg$ltr.end.dir <- substr(d_agg$ltr.end, 2, 2)

  # Summarise each group and keep those that look like a complete insertion:
  # - nprimes == 2: both the 5' and the 3' LTR terminus are seen, i.e. both
  #   flanks of the inserted element are supported
  # - ndir == 1: all junctions agree on the element's orientation
  # - max_mapq > 39: at least one molecule maps uniquely (MAPQ >= 40)
  # Per group this also records:
  # - dir: orientation of the element (F/R; becomes strand +/- later)
  # - tsd: distance between the outermost junctions, an estimate of the
  #   target-site duplication length (not filtered here; README step 5)
  # - first_bp: leftmost junction, used as the site's coordinate/ID
  # The per-junction rows of passing groups are joined back on (one row per
  # junction), so n and score are still split by LTR end at this point.
  d_final <- d_agg |>
    group_by(group) |>
    summarise(
      nprimes = length(unique(ltr.end.prime)),
      ndir = length(unique(ltr.end.dir)),
      dir = dplyr::first(ltr.end.dir),
      tsd = max(breakpoint) - min(breakpoint),
      max_mapq = max(mapq),
      first_bp = min(breakpoint)
    ) |>
    filter(nprimes == 2 & ndir == 1) |>
    filter(max_mapq > 39) |>
    left_join(d_agg)

  return(d_final)
}

## Aggregate breakpoints ------------------------------------------------------
# Run novel_sites() on every sample's tsv and stack the results into one
# table, d. The sample (Expt.ID) is extracted from each tsv filename by
# stripping the ".isa.tsv" suffix. isa_hmmer2.py always names its output
# "<bam-basename>.isa.tsv" (only --outdir is configurable, not the
# filename/extension), so this is guaranteed to work on any isa_hmmer2.py
# output, not just this repo's example.
#
# Parsing many large tsv files is slow, so the stacked table can be cached in
# an rds file (--cache-rds). If the cache exists it is used as is: delete it
# after adding/changing tsv files, otherwise they are silently ignored.

if (!is.null(opts$cache_rds) && file.exists(opts$cache_rds)) {
  message(sprintf("Loading cached breakpoints from %s", opts$cache_rds))
  d <- readRDS(opts$cache_rds)
} else {
  if (!is.null(opts$tsv_dir)) {
    tsv_files <- list.files(
      opts$tsv_dir,
      pattern = "\\.isa\\.tsv$",
      full.names = TRUE
    )
    if (length(tsv_files) == 0) {
      stop(
        sprintf("No *.isa.tsv files found in %s", opts$tsv_dir),
        call. = FALSE
      )
    }
  } else {
    tsv_files <- opts$tsv
    missing_files <- tsv_files[!file.exists(tsv_files)]
    if (length(missing_files) > 0) {
      stop(
        sprintf(
          "tsv file(s) not found: %s",
          paste(missing_files, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  d <- NULL
  for (f in tsv_files) {
    message(sprintf("Processing %s", f))
    d_file <- novel_sites(f)
    d_file$Expt.ID <- gsub("\\.isa\\.tsv$", "", basename(f))
    d <- rbind(d, d_file)
  }
  if (!is.null(opts$cache_rds)) {
    saveRDS(d, opts$cache_rds)
  }
}

## Summarise -------------------------------------------------------------
# One row per site and sample. A site is identified by chromosome + leftmost
# junction (name = "seqname:first_bp"). Per site and sample:
# - score: summed LTR score over all its junctions
# - n: unique molecules at the 5' junction plus unique molecules at the 3'
#   junction; this is the insertion read support used as the numerator of
#   the VAF (calculate_vaf.R)
# - strand: orientation of the inserted element ("F" -> "+", "R" -> "-")
# Because the site ID is the exact leftmost coordinate, the same insertion
# called at a slightly different first_bp in another sample gets a separate
# name and is not merged.

d_all <- d |>
  group_by(seqname, first_bp, Expt.ID) |>
  summarise(
    score = sum(score),
    n = sum(n),
    strand = ifelse(dplyr::first(dir) == "F", "+", "-")
  ) |>
  mutate(name = paste0(seqname, ":", first_bp))

# n_mice: number of samples in which a site is found. Sites in one sample
# only are candidate somatic insertions; sites in several samples are
# inherited/strain-specific (germline) insertions or recurrent artefacts.
d_all <- d_all |>
  group_by(name) |>
  summarise(n_mice = length(unique(Expt.ID))) |>
  full_join(d_all)

## Read metadata ------------------------------------------------------------
# Metadata xlsx must have an "Expt.ID" column linking to the tsv-derived
# experimental ID, used to join in other metadata (e.g. mouse, genotype).
# It must also have a "mouse" column, which is used as the x axis of the
# plots. The metadata is joined to the sites by Expt.ID only, and is checked
# first so the join cannot silently drop or duplicate sites:
# - Expt.ID must be present, non-empty and unique; it is read as text with
#   surrounding whitespace removed, so numeric IDs (e.g. 76) match the
#   tsv-derived "76"
# - every sample with sites must have a metadata row; otherwise its sites
#   would be missing from all outputs and n_mice would disagree with the
#   audit table
# - metadata columns may not share a name with the site table's columns
#   (name, seqname, first_bp, score, n, strand, n_mice), as these would be
#   duplicated with .x/.y suffixes in the outputs

edf <- read.xlsx(opts$metadata)

if (!("Expt.ID" %in% names(edf))) {
  stop(
    sprintf(
      "Metadata %s has no Expt.ID column; columns are: %s",
      opts$metadata,
      paste(names(edf), collapse = ", ")
    ),
    call. = FALSE
  )
}

edf$Expt.ID <- trimws(as.character(edf$Expt.ID))
if (any(is.na(edf$Expt.ID) | edf$Expt.ID == "")) {
  stop(
    sprintf("Metadata %s has rows with an empty Expt.ID", opts$metadata),
    call. = FALSE
  )
}

dup_ids <- unique(edf$Expt.ID[duplicated(edf$Expt.ID)])
if (length(dup_ids) > 0) {
  stop(
    sprintf(
      "Expt.ID must be unique in the metadata; duplicated: %s",
      paste(dup_ids, collapse = ", ")
    ),
    call. = FALSE
  )
}

missing_ids <- setdiff(unique(d_all$Expt.ID), edf$Expt.ID)
if (length(missing_ids) > 0) {
  stop(
    sprintf(
      paste(
        "Sample(s) with insertion sites but no metadata row: %s",
        "Add them to %s, or leave their tsv files out of the input.",
        sep = "\n"
      ),
      paste(missing_ids, collapse = ", "),
      opts$metadata
    ),
    call. = FALSE
  )
}

clashing <- setdiff(intersect(names(edf), names(d_all)), "Expt.ID")
if (length(clashing) > 0) {
  stop(
    sprintf(
      "Metadata column(s) clash with site table columns, rename them: %s",
      paste(clashing, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (!is.null(opts$facet_column) && !(opts$facet_column %in% names(edf))) {
  stop(
    sprintf(
      "--facet-column '%s' not found in metadata columns: %s",
      opts$facet_column,
      paste(names(edf), collapse = ", ")
    ),
    call. = FALSE
  )
}

## Draw raw output and write the aggregated / audit tables -------------------
# Dot plots of site (y) by mouse (x), coloured by log10 read support (n);
# one plot for single-sample (somatic) sites and one for sites shared by
# several samples (germline), optionally faceted by a metadata column.

somatic_plot <- ggplot(
  d_all |> filter(n_mice == 1) |> inner_join(edf, by = "Expt.ID"),
  aes(x = mouse, y = name, color = log10(n))
) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
if (!is.null(opts$facet_column)) {
  somatic_plot <- somatic_plot +
    facet_grid(
      as.formula(paste0("~", opts$facet_column)),
      scales = "free_x",
      space = "free_x"
    )
}
ggsave(opts$somatic_pdf, plot = somatic_plot, width = 12, height = 24)

if (sum(d_all$n_mice > 1)) {
  germline_plot <- ggplot(
    d_all |> filter(n_mice > 1) |> inner_join(edf, by = "Expt.ID"),
    aes(x = mouse, y = name, color = log10(n))
  ) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!is.null(opts$facet_column)) {
    germline_plot <- germline_plot +
      facet_grid(
        as.formula(paste0("~", opts$facet_column)),
        scales = "free_x",
        space = "free_x"
      )
  }
  ggsave(opts$germline_pdf, plot = germline_plot, width = 12, height = 8)
} else {
  message(sprintf(
    "No germline (multi-mouse) sites found, skipping %s",
    opts$germline_pdf
  ))
}

# Insertion table (isa_hmmer2_out.xlsx in the README): one row per site and
# sample with read support (n, score) plus the sample's metadata.
openxlsx::write.xlsx(d_all |> inner_join(edf, by = "Expt.ID"), opts$output_xlsx)

# Audit table (isa.audit.xlsx in the README): one row per site, to be
# validated manually (README step 5) by adding "classification" and "source"
# columns. Its "mouse" column holds the Expt.ID of the sample with the site,
# not the metadata "mouse" column, or "multi" when the site is found in more
# than one sample; downstream scripts (collect_evidence.py,
# 06_cluster_source_coverage.R) use it to find the sample's BAM.
openxlsx::write.xlsx(
  d_all |>
    inner_join(edf, by = "Expt.ID") |>
    group_by(name, seqname, first_bp, strand) |>
    summarise(
      nmice = length(unique(Expt.ID)),
      firstmouse = dplyr::first(Expt.ID)
    ) |>
    mutate(mouse = ifelse(nmice > 1, "multi", as.character(firstmouse))) |>
    select(name, seqname, strand, first_bp, mouse),
  opts$audit_xlsx
)

message("Done.")
