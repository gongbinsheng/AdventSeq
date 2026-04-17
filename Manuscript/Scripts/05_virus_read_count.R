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
  count_method = "reads",
  scale_factor = "1000000"
))
args$scale_factor <- as.numeric(args$scale_factor)

if (!(args$count_method %in% c("reads", "pairs"))) {
  stop("--count-method must be one of: reads, pairs")
}

metadata <- read_sample_metadata(args$sample_info)
out_dirs <- prepare_output_dirs(args$output_dir)
norm_df <- build_normalization_denominators(args$results_dir)
read_count <- read_read_count_data(args$results_dir, count_method = args$count_method)

detail <- read_count$detail
summary_df <- read_count$summary

if (!nrow(detail)) {
  stop("No in-scope read_count files were found.")
}

detail$taxon_name <- detail$label
detail$abundance_raw <- detail$count_number
detail <- detail[, c(
  "sample_id",
  "host_mapper",
  "virus_mapper",
  "virus_reference",
  "count_method",
  "taxon_name",
  "canonical_virus",
  "abundance_raw",
  "path"
)]

targeted_keys <- unique(summary_df[summary_df$virus_reference == "7viruses", c(
  "sample_id",
  "host_mapper",
  "virus_mapper",
  "virus_reference",
  "count_method"
)])
if (nrow(targeted_keys)) {
  targeted_expected <- merge(
    targeted_keys,
    data.frame(canonical_virus = TARGETED_VIRUS_LEVELS, stringsAsFactors = FALSE),
    all = TRUE
  )
  targeted_observed <- detail[detail$virus_reference == "7viruses", , drop = FALSE]
  targeted_full <- merge(
    targeted_expected,
    targeted_observed,
    by = c("sample_id", "host_mapper", "virus_mapper", "virus_reference", "count_method", "canonical_virus"),
    all.x = TRUE,
    sort = FALSE
  )
  targeted_full$taxon_name[is.na(targeted_full$taxon_name)] <- targeted_full$canonical_virus[is.na(targeted_full$taxon_name)]
  targeted_full$abundance_raw[is.na(targeted_full$abundance_raw)] <- 0
  untargeted_detail <- detail[detail$virus_reference == "RVDBv31", , drop = FALSE]
  detail <- rbind(targeted_full, untargeted_detail)
}

detail <- add_normalization_columns(detail, norm_df, scale_factor = args$scale_factor, abundance_col = "abundance_raw")
detail <- add_metadata(detail, metadata)

summary_df$abundance_raw <- summary_df$total_virus
summary_df <- add_normalization_columns(summary_df, norm_df, scale_factor = args$scale_factor, abundance_col = "abundance_raw")
summary_df <- add_metadata(summary_df, metadata)

targeted_detail <- detail[detail$virus_reference == "7viruses" & detail$include_analysis, , drop = FALSE]
untargeted_detail <- detail[detail$virus_reference == "RVDBv31" & detail$include_analysis, , drop = FALSE]

targeted_sample_summary <- stats::aggregate(
  cbind(abundance_raw, norm_total_reads, norm_non_host_reads) ~ sample_id + host_mapper + virus_mapper +
    count_method + library_type + experiment_protocol + truth_group,
  data = targeted_detail,
  FUN = sum
)
targeted_detected <- stats::aggregate(
  abundance_raw ~ sample_id + host_mapper + virus_mapper + count_method,
  data = transform(targeted_detail, abundance_raw = abundance_raw > 0),
  FUN = sum
)
names(targeted_detected)[names(targeted_detected) == "abundance_raw"] <- "viruses_detected"
targeted_sample_summary <- merge(
  targeted_sample_summary,
  targeted_detected,
  by = c("sample_id", "host_mapper", "virus_mapper", "count_method"),
  all.x = TRUE,
  sort = FALSE
)

targeted_detection_summary <- stats::aggregate(
  cbind(sample_detected = abundance_raw > 0, viruses_detected) ~ truth_group + host_mapper + virus_mapper +
    experiment_protocol,
  data = targeted_sample_summary,
  FUN = mean
)
names(targeted_detection_summary)[names(targeted_detection_summary) == "sample_detected"] <- "sample_detection_rate"

untargeted_sample_summary <- stats::aggregate(
  cbind(abundance_raw, norm_total_reads, norm_non_host_reads) ~ sample_id + host_mapper + virus_mapper +
    count_method + library_type + experiment_protocol + truth_group,
  data = untargeted_detail,
  FUN = sum
)
taxa_detected <- stats::aggregate(
  abundance_raw ~ sample_id + host_mapper + virus_mapper + count_method,
  data = transform(untargeted_detail, abundance_raw = abundance_raw > 0),
  FUN = sum
)
names(taxa_detected)[names(taxa_detected) == "abundance_raw"] <- "taxa_detected"
untargeted_sample_summary <- merge(
  untargeted_sample_summary,
  taxa_detected,
  by = c("sample_id", "host_mapper", "virus_mapper", "count_method"),
  all.x = TRUE,
  sort = FALSE
)

untargeted_top_hits <- top_n_by_group(
  untargeted_detail,
  group_cols = c("sample_id", "host_mapper", "virus_mapper"),
  value_col = "norm_total_reads",
  n = 5
)

write_tsv(detail, file.path(out_dirs$tables, "05_read_count_records.tsv"))
write_tsv(summary_df, file.path(out_dirs$tables, "05_read_count_sample_totals.tsv"))
write_tsv(targeted_sample_summary, file.path(out_dirs$tables, "05_targeted_read_count_summary.tsv"))
write_tsv(targeted_detection_summary, file.path(out_dirs$tables, "05_targeted_detection_summary.tsv"))
write_tsv(untargeted_sample_summary, file.path(out_dirs$tables, "05_untargeted_read_count_summary.tsv"))
write_tsv(untargeted_top_hits, file.path(out_dirs$tables, "05_untargeted_top_hits.tsv"))

plot_grouped_boxplot(
  targeted_sample_summary$norm_total_reads,
  paste(
    targeted_sample_summary$truth_group,
    targeted_sample_summary$experiment_protocol,
    targeted_sample_summary$host_mapper,
    targeted_sample_summary$virus_mapper,
    sep = "__"
  ),
  file.path(out_dirs$plots, "05_targeted_norm_total_reads_boxplot.pdf"),
  main = sprintf("Targeted read_count abundance (%s)", args$count_method),
  ylab = "Normalized abundance per 1e6 reads"
)

plot_grouped_boxplot(
  untargeted_sample_summary$taxa_detected,
  paste(
    untargeted_sample_summary$experiment_protocol,
    untargeted_sample_summary$host_mapper,
    untargeted_sample_summary$virus_mapper,
    sep = "__"
  ),
  file.path(out_dirs$plots, "05_untargeted_taxa_detected_boxplot.pdf"),
  main = sprintf("Untargeted taxa detected (%s)", args$count_method),
  ylab = "Detected taxa"
)

log_lines <- c(
  sprintf("read_count detail rows written: %d", nrow(detail)),
  sprintf("Targeted detail rows written: %d", nrow(targeted_detail)),
  sprintf("Untargeted detail rows written: %d", nrow(untargeted_detail)),
  sprintf("Rows with missing total-read denominator: %d", sum(is.na(detail$denom_total_reads))),
  sprintf("Rows with missing/non-positive non-host denominator: %d", sum(is.na(detail$norm_non_host_reads))),
  sprintf("Count method used: %s", args$count_method)
)
write_log_lines(log_lines, file.path(out_dirs$logs, "05_virus_read_count.log"))
