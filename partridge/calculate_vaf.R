# Estimate the variant allele frequency (VAF) of validated novel IAP
# insertions by normalising insertion support to source-locus coverage
# (README step 6).
#
# VAF = n / source_coverage, where
#   n               = unique UMIs supporting the insertion junctions
#                     (isa_hmmer2_out.xlsx, step 3)
#   source_coverage = unique UMIs spanning the source IAP's junctions
#                     (isa.audit.source_coverage.xlsx, step 6)
# The source copy is fixed (homozygous) in inbred strains and is captured with
# the same efficiency as the sequence-identical novel copy, so the ratio
# estimates the fraction of alleles carrying the insertion.
#
# The confidence interval is the exact conditional CI for the ratio of two
# Poisson counts from stats::poisson.test (Clopper-Pearson interval on
# n / (n + source_coverage), transformed to the ratio scale).
#
# Usage:
#   Rscript partridge/calculate_vaf.R \
#     --insertions test_data/isa_hmmer2_out.xlsx \
#     --source-coverage test_data/isa.audit.source_coverage.xlsx \
#     --output test_data/isa.vaf.xlsx [--conf-level 0.95]

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

parse_args <- function(args) {
  defaults <- list(conf_level = "0.95")
  required <- c("insertions", "source_coverage", "output")
  if (length(args) %% 2 != 0 || !all(startsWith(args[c(TRUE, FALSE)], "--"))) {
    stop("arguments must be given as --key value pairs")
  }
  keys <- gsub("-", "_", sub("^--", "", args[c(TRUE, FALSE)]))
  opts <- modifyList(defaults, as.list(setNames(args[c(FALSE, TRUE)], keys)))
  missing <- setdiff(required, names(opts))
  if (length(missing) > 0) {
    stop(
      "missing required argument(s): ",
      paste0("--", gsub("_", "-", missing), collapse = ", ")
    )
  }
  opts$conf_level <- as.numeric(opts$conf_level)
  opts
}

# Exact CI for the ratio n / source_coverage; NA when there is no usable
# source coverage
vaf_ci <- function(n, source_coverage, conf_level) {
  if (is.na(n) || is.na(source_coverage) || !is.finite(source_coverage)) {
    return(c(NA_real_, NA_real_))
  }
  if (source_coverage == 0) {
    return(c(NA_real_, NA_real_))
  }
  poisson.test(
    c(n, source_coverage),
    T = c(1, 1),
    conf.level = conf_level
  )$conf.int[1:2]
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

insertions <- read.xlsx(opts$insertions) |>
  select(name, Expt.ID, n, score)

# In the audit file "mouse" holds the Expt.ID; "multi" rows are placeholders
# (source_coverage = Inf) that 06_cluster_source_coverage.R expands into one
# row per sample
audit <- read.xlsx(opts$source_coverage) |>
  mutate(classification = as.logical(classification)) |>
  filter(classification %in% TRUE, mouse != "multi") |>
  rename(Expt.ID = mouse)

# Samples without a call at an insertion (possible for expanded "multi"
# insertions) have no supporting reads, so n = 0
vaf <- audit |>
  left_join(insertions, by = c("name", "Expt.ID")) |>
  mutate(
    n = coalesce(n, 0),
    source_coverage = as.numeric(source_coverage),
    vaf = ifelse(
      is.finite(source_coverage) & source_coverage > 0,
      n / source_coverage,
      NA_real_
    )
  )

ci <- t(mapply(
  vaf_ci,
  vaf$n,
  vaf$source_coverage,
  MoreArgs = list(conf_level = opts$conf_level)
))
vaf$vaf_ci_lower <- ci[, 1]
vaf$vaf_ci_upper <- ci[, 2]

no_coverage <- sum(is.na(vaf$vaf))
if (no_coverage > 0) {
  warning(
    no_coverage,
    " insertion/sample row(s) have missing or zero source coverage; ",
    "VAF set to NA"
  )
}

write.xlsx(vaf, opts$output)
message("Wrote ", nrow(vaf), " row(s) to ", opts$output)
