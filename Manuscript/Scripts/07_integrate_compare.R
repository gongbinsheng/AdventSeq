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

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", check.names = FALSE)
}

read_count_records <- read_tsv(file.path(out_dirs$tables, "05_read_count_records.tsv"))
targeted_read_count <- read_tsv(file.path(out_dirs$tables, "05_targeted_read_count_summary.tsv"))
viraquant_records <- read_tsv(file.path(out_dirs$tables, "06_viraquant_records.tsv"))
viraquant_sample <- read_tsv(file.path(out_dirs$tables, "06_viraquant_sample_summary.tsv"))
host_feature_summary <- read_tsv(file.path(out_dirs$tables, "04_featurecounts_summary.tsv"))

read_count_records <- read_count_records[read_count_records$virus_reference == "7viruses" & read_count_records$include_analysis, , drop = FALSE]
viraquant_records <- viraquant_records[viraquant_records$include_analysis, , drop = FALSE]
targeted_read_count <- targeted_read_count[targeted_read_count$truth_group %in% c("negative_control", "positive_titration"), , drop = FALSE]
viraquant_sample <- viraquant_sample[viraquant_sample$truth_group %in% c("negative_control", "positive_titration"), , drop = FALSE]
host_primary <- host_feature_summary[
  host_feature_summary$include_analysis &
    host_feature_summary$count_kind == "count_pairs",
  c("sample_id", "host_mapper", "assigned_fraction", "experiment_protocol", "library_type", "truth_group"),
  drop = FALSE
]

read_count_records$virus <- read_count_records$canonical_virus
rc_vs_vq <- merge(
  read_count_records[, c(
    "sample_id",
    "host_mapper",
    "virus_mapper",
    "virus",
    "norm_total_reads",
    "norm_non_host_reads",
    "abundance_raw",
    "experiment_protocol",
    "truth_group"
  )],
  viraquant_records[, c(
    "sample_id",
    "host_mapper",
    "virus_mapper",
    "virus",
    "norm_total_reads",
    "norm_non_host_reads",
    "abundance_raw",
    "selected_abundance_field"
  )],
  by = c("sample_id", "host_mapper", "virus_mapper", "virus"),
  all = FALSE,
  suffixes = c("_read_count", "_viraquant"),
  sort = FALSE
)

if (nrow(rc_vs_vq)) {
  concordance_parts <- split(
    rc_vs_vq,
    interaction(rc_vs_vq$host_mapper, rc_vs_vq$virus_mapper, rc_vs_vq$selected_abundance_field, drop = TRUE)
  )
  virus_level_concordance <- do.call(rbind, lapply(concordance_parts, function(df) {
    data.frame(
      host_mapper = df$host_mapper[[1]],
      virus_mapper = df$virus_mapper[[1]],
      selected_abundance_field = df$selected_abundance_field[[1]],
      spearman_total = safe_cor(log1p(df$norm_total_reads_read_count), log1p(df$norm_total_reads_viraquant)),
      spearman_non_host = safe_cor(log1p(df$norm_non_host_reads_read_count), log1p(df$norm_non_host_reads_viraquant)),
      spearman_raw = safe_cor(log1p(df$abundance_raw_read_count), log1p(df$abundance_raw_viraquant)),
      matched_rows = nrow(df),
      stringsAsFactors = FALSE
    )
  }))
} else {
  virus_level_concordance <- data.frame()
}

targeted_read_count$method <- "read_count"
targeted_read_count$primary_signal <- targeted_read_count$norm_total_reads
viraquant_sample$method <- "ViraQuant"
viraquant_sample$primary_signal <- viraquant_sample$primary_normalized_abundance

method_detection <- rbind(
  targeted_read_count[, c(
    "sample_id",
    "host_mapper",
    "virus_mapper",
    "experiment_protocol",
    "truth_group",
    "method",
    "primary_signal"
  )],
  viraquant_sample[, c(
    "sample_id",
    "host_mapper",
    "virus_mapper",
    "experiment_protocol",
    "truth_group",
    "method",
    "primary_signal"
  )]
)
method_detection$detected <- method_detection$primary_signal > 0

method_detection_summary <- stats::aggregate(
  cbind(detected, primary_signal) ~ method + truth_group + host_mapper + virus_mapper + experiment_protocol,
  data = method_detection,
  FUN = mean
)
names(method_detection_summary)[names(method_detection_summary) == "detected"] <- "detection_rate"

host_tradeoff_read_count <- merge(
  host_primary,
  targeted_read_count[, c("sample_id", "host_mapper", "virus_mapper", "norm_total_reads")],
  by = c("sample_id", "host_mapper"),
  all = FALSE,
  sort = FALSE
)
host_tradeoff_read_count$method <- "read_count"
host_tradeoff_read_count$viral_signal <- host_tradeoff_read_count$norm_total_reads
host_tradeoff_read_count$source_signal_column <- "norm_total_reads"

host_tradeoff_viraquant <- merge(
  host_primary,
  viraquant_sample[, c("sample_id", "host_mapper", "virus_mapper", "primary_normalized_abundance")],
  by = c("sample_id", "host_mapper"),
  all = FALSE,
  sort = FALSE
)
host_tradeoff_viraquant$method <- "ViraQuant"
host_tradeoff_viraquant$viral_signal <- host_tradeoff_viraquant$primary_normalized_abundance
host_tradeoff_viraquant$source_signal_column <- "primary_normalized_abundance"

host_tradeoff <- rbind(
  host_tradeoff_read_count[, c(
    "sample_id",
    "host_mapper",
    "assigned_fraction",
    "experiment_protocol",
    "library_type",
    "truth_group",
    "virus_mapper",
    "method",
    "viral_signal",
    "source_signal_column"
  )],
  host_tradeoff_viraquant[, c(
    "sample_id",
    "host_mapper",
    "assigned_fraction",
    "experiment_protocol",
    "library_type",
    "truth_group",
    "virus_mapper",
    "method",
    "viral_signal",
    "source_signal_column"
  )]
)

host_tradeoff_cor <- split(host_tradeoff, interaction(host_tradeoff$method, host_tradeoff$host_mapper, host_tradeoff$virus_mapper, drop = TRUE))
host_tradeoff_cor <- lapply(host_tradeoff_cor, function(df) {
  data.frame(
    method = df$method[[1]],
    host_mapper = df$host_mapper[[1]],
    virus_mapper = df$virus_mapper[[1]],
    spearman_cor = safe_cor(df$assigned_fraction, log1p(df$viral_signal)),
    n_rows = nrow(df),
    stringsAsFactors = FALSE
  )
})
host_tradeoff_cor <- do.call(rbind, host_tradeoff_cor)

write_tsv(rc_vs_vq, file.path(out_dirs$tables, "07_read_count_viraquant_matched.tsv"))
write_tsv(virus_level_concordance, file.path(out_dirs$tables, "07_virus_level_concordance.tsv"))
write_tsv(method_detection_summary, file.path(out_dirs$tables, "07_targeted_method_detection_summary.tsv"))
write_tsv(host_tradeoff, file.path(out_dirs$tables, "07_host_viral_tradeoff.tsv"))
write_tsv(host_tradeoff_cor, file.path(out_dirs$tables, "07_host_viral_tradeoff_correlations.tsv"))

if (nrow(rc_vs_vq)) {
  with_pdf(file.path(out_dirs$plots, "07_read_count_vs_viraquant_scatter.pdf"), {
    graphics::plot(
      log1p(rc_vs_vq$norm_total_reads_read_count),
      log1p(rc_vs_vq$norm_total_reads_viraquant),
      pch = 19,
      col = as.integer(as.factor(rc_vs_vq$host_mapper)),
      xlab = "log1p(read_count norm_total_reads)",
      ylab = "log1p(ViraQuant norm_total_reads)",
      main = "Targeted virus concordance: read_count vs ViraQuant"
    )
    graphics::legend(
      "topleft",
      legend = levels(as.factor(rc_vs_vq$host_mapper)),
      col = seq_along(levels(as.factor(rc_vs_vq$host_mapper))),
      pch = 19,
      cex = 0.8
    )
    graphics::grid()
  })
}

if (nrow(host_tradeoff)) {
  with_pdf(file.path(out_dirs$plots, "07_host_vs_viral_tradeoff_scatter.pdf"), {
    graphics::plot(
      host_tradeoff$assigned_fraction,
      log1p(host_tradeoff$viral_signal),
      pch = 19,
      col = as.integer(as.factor(host_tradeoff$method)),
      xlab = "Host assigned fraction",
      ylab = "log1p(viral signal)",
      main = "Host assignment versus viral recovery"
    )
    graphics::legend(
      "bottomleft",
      legend = levels(as.factor(host_tradeoff$method)),
      col = seq_along(levels(as.factor(host_tradeoff$method))),
      pch = 19,
      cex = 0.8
    )
    graphics::grid()
  })
}

log_lines <- c(
  sprintf("Matched read_count/ViraQuant virus rows: %d", nrow(rc_vs_vq)),
  sprintf("Method detection summary rows: %d", nrow(method_detection_summary)),
  sprintf("Host tradeoff rows: %d", nrow(host_tradeoff))
)
write_log_lines(log_lines, file.path(out_dirs$logs, "07_integrate_compare.log"))
