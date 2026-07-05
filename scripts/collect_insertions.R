#!/usr/bin/env Rscript

# Standalone, non-interactive replacement for 03_collect_insertions.Rmd.
# Aggregates the tsv files produced by isa_hmmer2.py into a set of novel
# insertion sites, joins them against a metadata table, and writes out
# summary plots and Excel audit files.

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
  "  --cache-rds FILE    Path to the aggregated-breakpoints rds cache. If it",
  "                      already exists, it is loaded instead of re-parsing",
  "                      the tsv files; otherwise it is created.",
  "  --somatic-pdf FILE  Output path for the somatic (single-mouse) plot",
  "  --germline-pdf FILE Output path for the germline (multi-mouse) plot",
  "                      (only written if germline sites are found)",
  "  --output-xlsx FILE  Output path for the aggregated insertions table",
  "  --audit-xlsx FILE   Output path for the manual-review audit table",
  "",
  "Optional:",
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
  "cache_rds",
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
  # replace synthetised barcodes with true ones, for reads that have one
  nc_filter <- substr(d_file$rname, nchar(d_file$rname), nchar(d_file$rname))
  nc_filter <- nc_filter %in% c("A", "T", "C", "G")
  nc <- nchar(d_file$rname)[nc_filter]
  d_file$barcode[nc_filter] <- substr(d_file$rname[nc_filter], nc - 8, nc)

  # aggregate
  d_agg <- d_file |>
    filter(bp.delta < 1) |>
    group_by(ltr.end, seqname, breakpoint, barcode) |>
    summarise(score = max(ltr.score), mapq = max(mapq)) |>
    group_by(ltr.end, seqname, breakpoint) |>
    summarise(score = sum(score), n = n(), mapq = max(mapq)) |>
    arrange(seqname, breakpoint)

  # group in chunks of 150bp
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

  # split end setting
  d_agg$ltr.end.prime <- substr(d_agg$ltr.end, 1, 1)
  d_agg$ltr.end.dir <- substr(d_agg$ltr.end, 2, 2)

  # final filtering
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
# The sample (Expt.ID) is extracted from each tsv filename by stripping the
# ".isa.tsv" suffix. isa_hmmer2.py always names its output "<bam-basename>.isa.tsv"
# (only --outdir is configurable, not the filename/extension), so this is
# guaranteed to work on any isa_hmmer2.py output, not just this repo's example.

if (file.exists(opts$cache_rds)) {
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
  saveRDS(d, opts$cache_rds)
}

## Summarise -------------------------------------------------------------

d_all <- d |>
  group_by(seqname, first_bp, Expt.ID) |>
  summarise(
    score = sum(score),
    n = sum(n),
    strand = ifelse(dplyr::first(dir) == "F", "+", "-")
  ) |>
  mutate(name = paste0(seqname, ":", first_bp))

d_all <- d_all |>
  group_by(name) |>
  summarise(n_mice = length(unique(Expt.ID))) |>
  full_join(d_all)

## Read metadata ------------------------------------------------------------
# Metadata xlsx must have an "Expt.ID" column linking to the tsv-derived
# experimental ID, used to join in other metadata (e.g. mouse, genotype).

edf <- read.xlsx(opts$metadata)

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

somatic_plot <- ggplot(
  d_all |> filter(n_mice == 1) |> inner_join(edf),
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
    d_all |> filter(n_mice > 1) |> inner_join(edf),
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

openxlsx::write.xlsx(d_all |> inner_join(edf), opts$output_xlsx)

openxlsx::write.xlsx(
  d_all |>
    inner_join(edf) |>
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
