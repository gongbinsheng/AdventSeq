#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[[1]]), mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

args <- resolve_common_args(list(
  results_dir = "Results",
  sample_info = "",
  output_dir = "Analyses"
))
out_dirs <- prepare_output_dirs(args$output_dir)

metadata <- read_sample_metadata(args$sample_info)
stats_df <- read_stats_qc(args$results_dir)
stats_df <- add_metadata(stats_df, metadata)

included_stats <- stats_df[stats_df$include_analysis, , drop = FALSE]

summary_df <- summarize_group_metrics(
  included_stats,
  group_cols = c("host_mapper", "experiment_protocol"),
  value_cols = c(
    "mapped_rate",
    "primary_mapped_rate",
    "properly_paired_rate",
    "singleton_rate",
    "view_primary_mapped_pairs"
  )
)

write_tsv(stats_df, file.path(out_dirs$tables, "03_host_mapping_metrics.tsv"))
write_tsv(summary_df, file.path(out_dirs$tables, "03_host_mapping_summary.tsv"))

plot_grouped_boxplot(
  included_stats$mapped_rate,
  paste(included_stats$host_mapper, included_stats$experiment_protocol, sep = "__"),
  file.path(out_dirs$plots, "03_mapped_rate_by_mapper_protocol.pdf"),
  main = "Mapped read rate by host mapper and protocol",
  ylab = "Mapped rate"
)

plot_grouped_boxplot(
  included_stats$properly_paired_rate,
  paste(included_stats$host_mapper, included_stats$experiment_protocol, sep = "__"),
  file.path(out_dirs$plots, "03_properly_paired_rate_by_mapper_protocol.pdf"),
  main = "Properly paired rate by host mapper and protocol",
  ylab = "Properly paired rate"
)

log_lines <- c(
  sprintf("Host mapping rows parsed: %d", nrow(stats_df)),
  sprintf("Included host mapping rows: %d", nrow(included_stats)),
  sprintf("Rows missing with_itself_and_mate_mapped: %d", sum(is.na(stats_df$with_itself_and_mate_mapped))),
  sprintf("Unique host mappers represented: %d", length(unique(stats_df$host_mapper)))
)
write_log_lines(log_lines, file.path(out_dirs$logs, "03_host_mapping_qc.log"))
