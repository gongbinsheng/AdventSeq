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
fastp <- read_fastp_qc(args$results_dir)
norm_df <- build_normalization_denominators(args$results_dir)

fastp <- add_metadata(fastp, metadata)
norm_df <- add_metadata(norm_df, metadata)

included_fastp <- fastp[fastp$include_analysis, , drop = FALSE]
included_norm <- norm_df[norm_df$include_analysis, , drop = FALSE]

sample_norm_summary <- stats::aggregate(
  cbind(after_filtering_total_reads, before_filtering_total_reads) ~ sample_id + library_type +
    experiment_protocol + truth_group,
  data = included_fastp,
  FUN = function(x) x[[1]]
)
sample_norm_summary$after_vs_before_fraction <- with(
  sample_norm_summary,
  after_filtering_total_reads / before_filtering_total_reads
)

non_host_summary <- summarize_group_metrics(
  included_norm,
  group_cols = c("host_mapper", "experiment_protocol"),
  value_cols = c("denom_total_reads", "denom_non_host_reads")
)

write_tsv(fastp, file.path(out_dirs$tables, "02_fastp_summary.tsv"))
write_tsv(norm_df, file.path(out_dirs$tables, "02_normalization_denominators.tsv"))
write_tsv(sample_norm_summary, file.path(out_dirs$tables, "02_fastp_sample_summary.tsv"))
write_tsv(non_host_summary, file.path(out_dirs$tables, "02_non_host_denominator_summary.tsv"))

plot_grouped_boxplot(
  included_fastp$after_filtering_total_reads,
  included_fastp$experiment_protocol,
  file.path(out_dirs$plots, "02_total_reads_by_protocol.pdf"),
  main = "fastp after-filtering total reads by protocol",
  ylab = "Reads"
)

plot_grouped_boxplot(
  included_norm$denom_non_host_reads,
  paste(included_norm$host_mapper, included_norm$experiment_protocol, sep = "__"),
  file.path(out_dirs$plots, "02_non_host_reads_by_mapper_protocol.pdf"),
  main = "Secondary non-host denominator by mapper/protocol",
  ylab = "Reads"
)

log_lines <- c(
  sprintf("fastp rows parsed: %d", nrow(fastp)),
  sprintf("Normalization rows emitted: %d", nrow(norm_df)),
  sprintf("Rows missing total-read denominator: %d", sum(is.na(norm_df$denom_total_reads))),
  sprintf("Rows missing non-host denominator: %d", sum(is.na(norm_df$denom_non_host_reads))),
  sprintf("Rows with non-positive non-host denominator: %d", sum(!is.na(norm_df$denom_non_host_reads) & norm_df$denom_non_host_reads <= 0))
)
write_log_lines(log_lines, file.path(out_dirs$logs, "02_fastp_qc.log"))
