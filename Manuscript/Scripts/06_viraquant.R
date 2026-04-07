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
  output_dir = "Analyses",
  abundance_field = "mean_depth_ge_1",
  normalization_mode = "both",
  scale_factor = "1000000"
))
args$scale_factor <- as.numeric(args$scale_factor)

if (!(args$abundance_field %in% VIRAQUANT_FIELDS)) {
  stop("--abundance-field must be one of: ", paste(VIRAQUANT_FIELDS, collapse = ", "))
}
if (!(args$normalization_mode %in% c("both", "total_reads", "non_host_reads"))) {
  stop("--normalization-mode must be one of: both, total_reads, non_host_reads")
}

metadata <- read_sample_metadata(args$sample_info)
out_dirs <- prepare_output_dirs(args$output_dir)
norm_df <- build_normalization_denominators(args$results_dir)
vq <- read_viraquant_data(args$results_dir)

if (!nrow(vq)) {
  stop("No in-scope ViraQuant files were found.")
}

abundance_raw <- suppressWarnings(as.numeric(vq[[args$abundance_field]]))
abundance_status <- rep("ok", nrow(vq))

if (args$abundance_field %in% COVERAGE_ABUNDANCE_FIELDS) {
  zero_fill <- is.na(abundance_raw) & suppressWarnings(as.numeric(vq$mapped_reads)) == 0
  abundance_raw[zero_fill] <- 0
  abundance_status[zero_fill] <- "coerced_zero_from_na"
}
abundance_status[is.na(abundance_raw) & abundance_status == "ok"] <- "missing_selected_abundance"

vq$selected_abundance_field <- args$abundance_field
vq$abundance_raw <- abundance_raw
vq$abundance_field_status <- abundance_status

vq <- add_normalization_columns(vq, norm_df, scale_factor = args$scale_factor, abundance_col = "abundance_raw")
vq$record_status <- mapply(combine_flags, vq$abundance_field_status, vq$normalization_status, USE.NAMES = FALSE)
vq <- add_metadata(vq, metadata)

primary_norm_col <- switch(
  args$normalization_mode,
  total_reads = "norm_total_reads",
  non_host_reads = "norm_non_host_reads",
  both = "norm_total_reads"
)
vq$primary_normalization_mode <- args$normalization_mode
vq$primary_normalized_abundance <- vq[[primary_norm_col]]

included_vq <- vq[vq$include_analysis, , drop = FALSE]

sample_summary <- stats::aggregate(
  cbind(
    abundance_raw,
    norm_total_reads,
    norm_non_host_reads,
    primary_normalized_abundance,
    mapped_reads,
    pass_k_90pct_ge_1
  ) ~ sample_id + host_mapper + virus_mapper + selected_abundance_field +
    primary_normalization_mode + library_type + experiment_protocol + truth_group,
  data = included_vq,
  FUN = sum
)
names(sample_summary)[names(sample_summary) == "pass_k_90pct_ge_1"] <- "n_pass_k_90pct_ge_1"

quality_summary <- summarize_group_metrics(
  included_vq,
  group_cols = c("host_mapper", "virus_mapper", "selected_abundance_field", "experiment_protocol"),
  value_cols = c("frac_ge_1", "mean_depth_ge_1", "mapped_per_bp")
)

write_tsv(vq, file.path(out_dirs$tables, "06_viraquant_records.tsv"))
write_tsv(sample_summary, file.path(out_dirs$tables, "06_viraquant_sample_summary.tsv"))
write_tsv(quality_summary, file.path(out_dirs$tables, "06_viraquant_quality_summary.tsv"))

plot_grouped_boxplot(
  sample_summary$primary_normalized_abundance,
  paste(
    sample_summary$truth_group,
    sample_summary$experiment_protocol,
    sample_summary$host_mapper,
    sample_summary$virus_mapper,
    sep = "__"
  ),
  file.path(out_dirs$plots, "06_viraquant_primary_abundance_boxplot.pdf"),
  main = sprintf("ViraQuant %s (%s)", args$abundance_field, args$normalization_mode),
  ylab = "Primary normalized abundance"
)

plot_grouped_boxplot(
  included_vq$pass_k_90pct_ge_1,
  paste(included_vq$host_mapper, included_vq$virus_mapper, sep = "__"),
  file.path(out_dirs$plots, "06_viraquant_pass_k90_boxplot.pdf"),
  main = "ViraQuant pass_k_90pct_ge_1 by mapper combination",
  ylab = "pass_k_90pct_ge_1"
)

log_lines <- c(
  sprintf("ViraQuant rows written: %d", nrow(vq)),
  sprintf("Selected abundance field: %s", args$abundance_field),
  sprintf("Normalization mode: %s", args$normalization_mode),
  sprintf("Rows coerced from NA to zero: %d", sum(vq$abundance_field_status == 'coerced_zero_from_na')),
  sprintf("Rows with missing selected abundance after coercion: %d", sum(vq$abundance_field_status == 'missing_selected_abundance'))
)
write_log_lines(log_lines, file.path(out_dirs$logs, "06_viraquant.log"))
