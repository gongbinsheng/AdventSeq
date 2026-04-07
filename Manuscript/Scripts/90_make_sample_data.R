#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[[1]]), mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

args <- parse_cli_args(list(
  results_dir = "Results",
  sample_info = "",
  output_results_dir = file.path("sample_data", "Results"),
  featurecounts_rows = "1000"
))
args$results_dir <- normalizePath(args$results_dir, mustWork = FALSE)
if (is.null(args$sample_info) || !nzchar(args$sample_info)) {
  args$sample_info <- file.path(args$results_dir, "sample_info.tsv")
}
args$sample_info <- normalizePath(args$sample_info, mustWork = FALSE)
args$output_results_dir <- normalizePath(args$output_results_dir, mustWork = FALSE)
args$featurecounts_rows <- as.integer(args$featurecounts_rows)

selected_ids <- c("1_S1", "8_S8", "19_S3", "28_S6", "31_S9", "A6_S6", "A7_S7", "4_S4", "11_S11", "13_S13")

copy_featurecounts_matrix <- function(src, dest, n_data_rows = 1000) {
  ensure_dir(dirname(dest))
  in_con <- gzfile(src, open = "rt")
  on.exit(close(in_con), add = TRUE)
  header <- readLines(in_con, n = 2, warn = FALSE)
  body <- readLines(in_con, n = n_data_rows, warn = FALSE)

  out_con <- gzfile(dest, open = "wt")
  on.exit(close(out_con), add = TRUE)
  writeLines(c(header, body), con = out_con, useBytes = TRUE)
}

copy_matching_files <- function(src_dir, dest_dir, truncate_featurecounts = FALSE) {
  ensure_dir(dest_dir)
  files <- list.files(src_dir, full.names = TRUE)
  keep <- vapply(files, function(path) {
    name <- basename(path)
    any(startsWith(name, paste0(selected_ids, ".")))
  }, logical(1))
  files <- files[keep]
  copied <- character(0)

  for (src in files) {
    dest <- file.path(dest_dir, basename(src))
    if (truncate_featurecounts) {
      copy_featurecounts_matrix(src, dest, n_data_rows = args$featurecounts_rows)
    } else {
      ensure_dir(dirname(dest))
      file.copy(src, dest, overwrite = TRUE)
    }
    copied <- c(copied, dest)
  }
  copied
}

metadata <- read_sample_metadata(args$sample_info)
subset_meta <- metadata[metadata$sample_id %in% selected_ids, c("library_ID", "title"), drop = FALSE]
subset_meta <- subset_meta[match(selected_ids, subset_meta$library_ID), , drop = FALSE]

ensure_dir(args$output_results_dir)
write_tsv(subset_meta, file.path(args$output_results_dir, "sample_info.tsv"))

copied_fastp <- copy_matching_files(
  file.path(args$results_dir, "fastp"),
  file.path(args$output_results_dir, "fastp")
)
copied_stats <- copy_matching_files(
  file.path(args$results_dir, "stats"),
  file.path(args$output_results_dir, "stats")
)
copied_featurecounts <- copy_matching_files(
  file.path(args$results_dir, "featureCounts"),
  file.path(args$output_results_dir, "featureCounts"),
  truncate_featurecounts = TRUE
)
copied_featurecounts_summary <- copy_matching_files(
  file.path(args$results_dir, "featureCounts", "summary"),
  file.path(args$output_results_dir, "featureCounts", "summary")
)
copied_read_count <- copy_matching_files(
  file.path(args$results_dir, "read_count"),
  file.path(args$output_results_dir, "read_count")
)
copied_viraquant <- copy_matching_files(
  file.path(args$results_dir, "ViraQuant"),
  file.path(args$output_results_dir, "ViraQuant")
)

fixture_manifest <- data.frame(
  source = c(
    "sample_info",
    rep("fastp", length(copied_fastp)),
    rep("stats", length(copied_stats)),
    rep("featureCounts_matrix", length(copied_featurecounts)),
    rep("featureCounts_summary", length(copied_featurecounts_summary)),
    rep("read_count", length(copied_read_count)),
    rep("ViraQuant", length(copied_viraquant))
  ),
  path = c(
    file.path(args$output_results_dir, "sample_info.tsv"),
    copied_fastp,
    copied_stats,
    copied_featurecounts,
    copied_featurecounts_summary,
    copied_read_count,
    copied_viraquant
  ),
  stringsAsFactors = FALSE
)

write_tsv(fixture_manifest, file.path(dirname(args$output_results_dir), "sample_data_manifest.tsv"))
write_log_lines(
  c(
    sprintf("Selected sample IDs: %s", paste(selected_ids, collapse = ", ")),
    sprintf("sample_info rows written: %d", nrow(subset_meta)),
    sprintf("fastp files copied: %d", length(copied_fastp)),
    sprintf("stats files copied: %d", length(copied_stats)),
    sprintf("featureCounts matrices truncated/copied: %d", length(copied_featurecounts)),
    sprintf("featureCounts summaries copied: %d", length(copied_featurecounts_summary)),
    sprintf("read_count files copied: %d", length(copied_read_count)),
    sprintf("ViraQuant files copied: %d", length(copied_viraquant))
  ),
  file.path(dirname(args$output_results_dir), "sample_data_creation.log")
)
