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
included <- metadata[metadata$include_analysis, , drop = FALSE]
excluded <- metadata[!metadata$include_analysis, , drop = FALSE]

expected_design <- expand.grid(
  library_type = sort(unique(included$library_type)),
  technical_replicate = sort(unique(included$technical_replicate)),
  experiment_protocol = sort(unique(included$experiment_protocol)),
  stringsAsFactors = FALSE
)
observed_design <- stats::aggregate(
  sample_id ~ library_type + technical_replicate + experiment_protocol,
  data = included,
  FUN = length
)
names(observed_design)[names(observed_design) == "sample_id"] <- "observed_samples"
design_table <- merge(
  expected_design,
  observed_design,
  by = c("library_type", "technical_replicate", "experiment_protocol"),
  all.x = TRUE,
  sort = FALSE
)
design_table$observed_samples[is.na(design_table$observed_samples)] <- 0
design_table$expected_samples <- 1L
design_table$design_status <- ifelse(design_table$observed_samples == 1L, "ok", "unexpected_count")

inventory <- bind_rows_fill(list(
  transform(list_fastp_files(args$results_dir), in_scope = TRUE),
  transform(list_stats_files(args$results_dir), in_scope = TRUE),
  list_featurecounts_matrices(args$results_dir),
  list_featurecounts_summaries(args$results_dir),
  list_read_count_files(args$results_dir, in_scope_only = FALSE),
  list_viraquant_files(args$results_dir, in_scope_only = FALSE)
))
inventory <- inventory[order(inventory$source, inventory$sample_id), , drop = FALSE]

fastp_expected <- data.frame(
  source = "fastp",
  sample_id = included$sample_id,
  stringsAsFactors = FALSE
)
fastp_actual <- list_fastp_files(args$results_dir)[, c("sample_id", "path"), drop = FALSE]
fastp_complete <- merge(fastp_expected, fastp_actual, by = "sample_id", all.x = TRUE, sort = FALSE)
fastp_complete$exists <- !is.na(fastp_complete$path)

stats_expected <- expand.grid(
  source = "stats",
  sample_id = included$sample_id,
  host_mapper = HOST_MAPPERS,
  stringsAsFactors = FALSE
)
stats_actual <- list_stats_files(args$results_dir)[, c("sample_id", "host_mapper", "path"), drop = FALSE]
stats_complete <- merge(
  stats_expected,
  stats_actual,
  by = c("sample_id", "host_mapper"),
  all.x = TRUE,
  sort = FALSE
)
stats_complete$exists <- !is.na(stats_complete$path)

feature_expected <- expand.grid(
  source = "featureCounts_matrix",
  sample_id = included$sample_id,
  host_mapper = HOST_MAPPERS,
  count_kind = c("count_pairs", "count_reads"),
  stringsAsFactors = FALSE
)
feature_actual <- list_featurecounts_matrices(args$results_dir)[
  , c("sample_id", "host_mapper", "count_kind", "path"),
  drop = FALSE
]
feature_complete <- merge(
  feature_expected,
  feature_actual,
  by = c("sample_id", "host_mapper", "count_kind"),
  all.x = TRUE,
  sort = FALSE
)
feature_complete$exists <- !is.na(feature_complete$path)

read_count_expected <- expand.grid(
  source = "read_count",
  sample_id = included$sample_id,
  host_mapper = HOST_MAPPERS,
  virus_mapper = VIRUS_MAPPERS,
  virus_reference = VIRUS_REFERENCES,
  stringsAsFactors = FALSE
)
read_count_actual <- list_read_count_files(args$results_dir, in_scope_only = TRUE)[
  , c("sample_id", "host_mapper", "virus_mapper", "virus_reference", "path"),
  drop = FALSE
]
read_count_complete <- merge(
  read_count_expected,
  read_count_actual,
  by = c("sample_id", "host_mapper", "virus_mapper", "virus_reference"),
  all.x = TRUE,
  sort = FALSE
)
read_count_complete$exists <- !is.na(read_count_complete$path)

viraquant_expected <- expand.grid(
  source = "ViraQuant",
  sample_id = included$sample_id,
  host_mapper = HOST_MAPPERS,
  virus_mapper = VIRUS_MAPPERS,
  stringsAsFactors = FALSE
)
viraquant_actual <- list_viraquant_files(args$results_dir, in_scope_only = TRUE)[
  , c("sample_id", "host_mapper", "virus_mapper", "path"),
  drop = FALSE
]
viraquant_complete <- merge(
  viraquant_expected,
  viraquant_actual,
  by = c("sample_id", "host_mapper", "virus_mapper"),
  all.x = TRUE,
  sort = FALSE
)
viraquant_complete$exists <- !is.na(viraquant_complete$path)

completeness <- bind_rows_fill(list(
  fastp_complete,
  stats_complete,
  feature_complete,
  read_count_complete,
  viraquant_complete
))

completeness_found <- stats::aggregate(
  exists ~ sample_id + source,
  data = completeness,
  FUN = sum
)
names(completeness_found)[names(completeness_found) == "exists"] <- "found"
completeness_expected <- stats::aggregate(
  exists ~ sample_id + source,
  data = completeness,
  FUN = length
)
names(completeness_expected)[names(completeness_expected) == "exists"] <- "expected"
completeness_fraction <- stats::aggregate(
  exists ~ sample_id + source,
  data = completeness,
  FUN = mean
)
names(completeness_fraction)[names(completeness_fraction) == "exists"] <- "fraction"
completeness_summary <- merge(
  merge(completeness_found, completeness_expected, by = c("sample_id", "source"), sort = FALSE),
  completeness_fraction,
  by = c("sample_id", "source"),
  sort = FALSE
)
completeness_summary$missing <- completeness_summary$expected - completeness_summary$found
completeness_summary <- add_metadata(completeness_summary, metadata)

write_tsv(metadata, file.path(out_dirs$tables, "01_sample_metadata.tsv"))
write_tsv(excluded, file.path(out_dirs$tables, "01_excluded_samples.tsv"))
write_tsv(design_table, file.path(out_dirs$tables, "01_design_table.tsv"))
write_tsv(inventory, file.path(out_dirs$tables, "01_file_inventory.tsv"))
write_tsv(completeness, file.path(out_dirs$tables, "01_file_completeness.tsv"))
write_tsv(completeness_summary, file.path(out_dirs$tables, "01_file_completeness_summary.tsv"))

protocol_counts <- table(included$experiment_protocol)
plot_simple_bar(
  as.numeric(protocol_counts),
  names(protocol_counts),
  file.path(out_dirs$plots, "01_included_samples_by_protocol.pdf"),
  main = "Included samples by protocol",
  ylab = "Samples"
)

heat_df <- completeness_summary[, c("sample_id", "source", "fraction")]
heat_mat <- xtabs(fraction ~ sample_id + source, data = heat_df)
with_pdf(file.path(out_dirs$plots, "01_missingness_heatmap.pdf"), {
  ord_rows <- order(rownames(heat_mat))
  ord_cols <- order(colnames(heat_mat))
  plot_mat <- heat_mat[ord_rows, ord_cols, drop = FALSE]
  graphics::image(
    x = seq_len(ncol(plot_mat)),
    y = seq_len(nrow(plot_mat)),
    z = t(plot_mat[nrow(plot_mat):1, , drop = FALSE]),
    axes = FALSE,
    col = grDevices::colorRampPalette(c("#b2182b", "#f7f7f7", "#2166ac"))(100),
    xlab = "Source",
    ylab = "Sample",
    main = "Prompt-scope completeness fraction"
  )
  graphics::axis(1, at = seq_len(ncol(plot_mat)), labels = colnames(plot_mat), las = 2)
  graphics::axis(2, at = seq_len(nrow(plot_mat)), labels = rev(rownames(plot_mat)), las = 2)
  graphics::box()
})

log_lines <- c(
  sprintf("Sample info rows: %d", nrow(metadata)),
  sprintf("Included samples: %d", nrow(included)),
  sprintf("Excluded U937 samples: %d", nrow(excluded)),
  sprintf("Protocols represented: %d", length(unique(included$experiment_protocol))),
  sprintf("Design rows with unexpected counts: %d", sum(design_table$design_status != 'ok')),
  sprintf("Completeness rows with missing files: %d", sum(!completeness$exists))
)
write_log_lines(log_lines, file.path(out_dirs$logs, "01_metadata_inventory.log"))
