#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(scales)
})

fig_dir <- "Figures"
tbl_dir <- "Tables"

dir.create(fig_dir, showWarnings = FALSE)
dir.create(tbl_dir, showWarnings = FALSE)

read_tsv_q <- function(path) {
  read_tsv(path, show_col_types = FALSE)
}

pick_existing_path <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (!length(existing)) {
    stop("None of the expected input files exist: ", paste(paths, collapse = ", "))
  }
  existing[[1]]
}

save_plot <- function(plot_obj, filename, width = 8, height = 5) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in"
  )
}

write_tsv_q <- function(df, filename) {
  write_tsv(df, file.path(tbl_dir, filename), na = "NA")
}

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    plot.title.position = "plot",
    legend.position = "bottom"
  )

si_labels <- label_number(scale_cut = cut_short_scale())

mapper_palette <- c(
  "bwa_mem" = "#1b3a4b",
  "HISAT2" = "#2a9d8f",
  "STAR" = "#e76f51"
)

truth_palette <- c(
  "negative_control" = "#6c757d",
  "positive_titration" = "#c1121f"
)

truth_labels <- c(
  "negative_control" = "Control",
  "positive_titration" = "Positive titration"
)

method_levels <- c(
  "bwa_mem / bowtie2",
  "bwa_mem / bwa_mem",
  "HISAT2 / bowtie2",
  "HISAT2 / bwa_mem",
  "STAR / bowtie2",
  "STAR / bwa_mem"
)

sample_meta <- read_tsv_q("Analyses/full_run/tables/01_sample_metadata.tsv")
fastp <- read_tsv_q("Analyses/full_run/tables/02_fastp_sample_summary.tsv")
denominators <- read_tsv_q("Analyses/full_run/tables/02_normalization_denominators.tsv") %>%
  filter(include_analysis)
gene <- read_tsv_q("Analyses/full_run/tables/04_gene_detection_summary.tsv") %>%
  filter(include_analysis, count_kind == "count_pairs")
mapcor <- read_tsv_q("Analyses/full_run/tables/04_mapper_concordance.tsv") %>%
  filter(include_analysis, count_kind == "count_pairs")
target <- read_tsv_q("Analyses/full_run/tables/05_targeted_read_count_summary.tsv")
untarget <- read_tsv_q("Analyses/full_run/tables/05_untargeted_read_count_summary.tsv")
vq <- read_tsv_q("Analyses/full_run/tables/06_viraquant_sample_summary.tsv")
vq_records <- read_tsv_q("Analyses/full_run/tables/06_viraquant_records.tsv") %>%
  filter(level == "virus")
matched <- read_tsv_q("Analyses/full_run/tables/07_read_count_viraquant_matched.tsv")
concordance <- read_tsv_q("Analyses/full_run/tables/07_virus_level_concordance.tsv")
kraken <- read_tsv_q(pick_existing_path(c(
  "Analyses/tables/08_kraken2_rvdbv31_viraquant_scan_taxid_summary.tsv",
  "Analyses/full_run/tables/08_kraken2_rvdbv31_viraquant_scan_taxid_summary.tsv"
))) %>%
  filter(include_analysis)
pca_scores <- read_tsv_q("Analyses/full_run/tables/04_host_pca_scores.tsv") %>%
  filter(include_analysis, count_kind == "count_pairs")
fastp_json_files <- list.files("Results/fastp", pattern = "\\.fastp\\.json$", full.names = TRUE)

parse_fastp_json <- function(path) {
  obj <- fromJSON(path)
  sample_id <- sub("\\.fastp\\.json$", "", basename(path))
  tibble(
    sample_id = sample_id,
    before_total_reads = obj$summary$before_filtering$total_reads,
    adapter_trimmed_reads = obj$adapter_cutting$adapter_trimmed_reads,
    polyx_trimmed_reads = obj$polyx_trimming$total_polyx_trimmed_reads,
    polyA_reads = obj$polyx_trimming$polyx_trimmed_reads[["A"]],
    polyT_reads = obj$polyx_trimming$polyx_trimmed_reads[["T"]],
    polyC_reads = obj$polyx_trimming$polyx_trimmed_reads[["C"]],
    polyG_reads = obj$polyx_trimming$polyx_trimmed_reads[["G"]],
    q20_before = obj$summary$before_filtering$q20_rate,
    q20_after = obj$summary$after_filtering$q20_rate,
    q30_before = obj$summary$before_filtering$q30_rate,
    q30_after = obj$summary$after_filtering$q30_rate
  )
}

fastp_extended <- map_dfr(fastp_json_files, parse_fastp_json) %>%
  inner_join(
    sample_meta %>%
      filter(include_analysis) %>%
      select(sample_id, experiment_protocol),
    by = "sample_id"
  ) %>%
  mutate(
    adapter_trimmed_fraction = adapter_trimmed_reads / before_total_reads,
    polyx_trimmed_fraction = polyx_trimmed_reads / before_total_reads,
    q20_gain = q20_after - q20_before,
    q30_gain = q30_after - q30_before,
    polyA_fraction = polyA_reads / polyx_trimmed_reads,
    polyT_fraction = polyT_reads / polyx_trimmed_reads,
    polyC_fraction = polyC_reads / polyx_trimmed_reads,
    polyG_fraction = polyG_reads / polyx_trimmed_reads
  )

protocol_summary <- sample_meta %>%
  filter(include_analysis) %>%
  group_by(experiment_protocol) %>%
  summarise(
    included_samples = n(),
    control_samples = sum(library_type == "Control"),
    inoculum_1to1_samples = sum(library_type == "1:1"),
    inoculum_10to1_samples = sum(library_type == "10:1"),
    .groups = "drop"
  ) %>%
  left_join(
    fastp %>%
      group_by(experiment_protocol) %>%
      summarise(
        median_after_filtering_reads = round(median(after_filtering_total_reads)),
        median_retained_fraction = round(median(after_vs_before_fraction), 3),
        .groups = "drop"
      ),
    by = "experiment_protocol"
  ) %>%
  left_join(
    denominators %>%
      group_by(experiment_protocol, host_mapper) %>%
      summarise(
        median_non_host_reads = round(median(denom_non_host_reads)),
        .groups = "drop"
      ) %>%
      mutate(host_mapper = paste0("median_non_host_reads_", host_mapper)) %>%
      pivot_wider(
        names_from = host_mapper,
        values_from = median_non_host_reads
      ),
    by = "experiment_protocol"
  ) %>%
  arrange(desc(median_after_filtering_reads))

write_tsv_q(protocol_summary, "Table1.tsv")

assigned_summary <- denominators %>%
  mutate(assigned_fraction = with_itself_and_mate_mapped / after_filtering_total_reads) %>%
  group_by(host_mapper) %>%
  summarise(
    median_host_assigned_fraction = round(median(assigned_fraction), 3),
    .groups = "drop"
  )

target_summary <- target %>%
  group_by(host_mapper, virus_mapper, truth_group) %>%
  summarise(
    median_target_norm_total_reads = round(median(norm_total_reads), 3),
    median_target_viruses_detected = round(median(viruses_detected), 2),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = truth_group,
    values_from = c(
      median_target_norm_total_reads,
      median_target_viruses_detected
    )
  )

untarget_summary <- untarget %>%
  group_by(host_mapper, virus_mapper, truth_group) %>%
  summarise(
    median_untarget_norm_total_reads = round(median(norm_total_reads), 3),
    median_untarget_taxa_detected = round(median(taxa_detected), 1),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = truth_group,
    values_from = c(
      median_untarget_norm_total_reads,
      median_untarget_taxa_detected
    )
  )

vq_summary <- vq %>%
  group_by(host_mapper, virus_mapper, truth_group) %>%
  summarise(
    median_viraquant_primary = round(median(primary_normalized_abundance), 3),
    median_viraquant_pass90 = round(median(n_pass_k_90pct_ge_1), 2),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = truth_group,
    values_from = c(median_viraquant_primary, median_viraquant_pass90)
  )

kraken_summary <- kraken %>%
  group_by(host_mapper, virus_mapper) %>%
  summarise(
    median_kraken_taxids = round(median(kraken_taxids_total), 1),
    median_rvdb_read_count_taxids = round(median(rvdb_taxids_total), 1),
    median_viraquant_scan_taxids = round(median(scan_taxids_total), 1),
    median_shared_all_three_taxids = round(median(shared_all_three_taxids), 1),
    median_kraken_only_taxids = round(median(kraken_only_taxids), 1),
    median_shared_all_three_rvdb_fraction = round(median(shared_all_three_rvdb_fraction), 3),
    .groups = "drop"
  )

method_summary <- target_summary %>%
  left_join(untarget_summary, by = c("host_mapper", "virus_mapper")) %>%
  left_join(vq_summary, by = c("host_mapper", "virus_mapper")) %>%
  left_join(kraken_summary, by = c("host_mapper", "virus_mapper")) %>%
  left_join(assigned_summary, by = "host_mapper") %>%
  arrange(match(host_mapper, c("bwa_mem", "HISAT2", "STAR")), virus_mapper)

write_tsv_q(method_summary, "Table2.tsv")

coverage_case <- matched %>%
  left_join(
    vq_records %>%
      select(
        sample_id,
        host_mapper,
        virus_mapper,
        virus,
        frac_ge_1,
        mean_depth_ge_1,
        pass_k_90pct_ge_1
      ),
    by = c("sample_id", "host_mapper", "virus_mapper", "virus")
  ) %>%
  filter(sample_id == "1_S1", virus %in% c("HEV", "Zika")) %>%
  mutate(
    method = factor(
      paste(host_mapper, virus_mapper, sep = " / "),
      levels = method_levels
    )
  ) %>%
  arrange(virus, method) %>%
  select(
    sample_id,
    experiment_protocol,
    host_mapper,
    virus_mapper,
    virus,
    abundance_raw_read_count,
    norm_total_reads_read_count,
    frac_ge_1,
    mean_depth_ge_1,
    pass_k_90pct_ge_1
  )

write_tsv_q(coverage_case, "Table3.tsv")

# Per-host-mapper summary (one row per host mapper). Pairwise mapper
# correlations are a property of mapper PAIRS, not of a single mapper, so they
# are reported separately (see mapper_pair_correlations below) rather than being
# cross-joined onto every row, which previously repeated the same value in each
# row and was confusing/redundant.
host_support_summary <- assigned_summary %>%
  left_join(
    gene %>%
      group_by(host_mapper) %>%
      summarise(
        median_genes_detected = round(median(genes_detected)),
        mean_genes_detected = round(mean(genes_detected)),
        .groups = "drop"
      ),
    by = "host_mapper"
  )

# Pairwise count_pairs mapper concordance (one row per mapper pair).
mapper_pair_correlations <- mapcor %>%
  group_by(mapper_a, mapper_b) %>%
  summarise(median_spearman = round(median(spearman_cor), 3), .groups = "drop") %>%
  mutate(mapper_pair = paste(mapper_a, mapper_b, sep = " vs ")) %>%
  select(mapper_pair, median_spearman)

write_tsv_q(host_support_summary, "Suppl_Table1.tsv")
write_tsv_q(mapper_pair_correlations, "Suppl_Table1_mapper_correlations.tsv")
write_tsv_q(concordance, "Suppl_Table2.tsv")
write_tsv_q(
  fastp_extended %>%
    group_by(experiment_protocol) %>%
    summarise(
      median_adapter_trimmed_reads_pct = round(median(adapter_trimmed_fraction) * 100, 2),
      median_polyx_trimmed_reads_pct = round(median(polyx_trimmed_fraction) * 100, 2),
      median_q20_gain_pct_points = round(median(q20_gain) * 100, 2),
      median_q30_gain_pct_points = round(median(q30_gain) * 100, 2),
      median_polyA_share_pct = round(median(polyA_fraction) * 100, 1),
      median_polyT_share_pct = round(median(polyT_fraction) * 100, 1),
      median_polyC_share_pct = round(median(polyC_fraction) * 100, 1),
      median_polyG_share_pct = round(median(polyG_fraction) * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(median_adapter_trimmed_reads_pct)),
  "Suppl_Table3.tsv"
)

plot_nonhost <- denominators %>%
  mutate(experiment_protocol = factor(experiment_protocol, levels = unique(protocol_summary$experiment_protocol))) %>%
  ggplot(aes(x = host_mapper, y = denom_non_host_reads, color = host_mapper)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.2) +
  geom_jitter(width = 0.15, alpha = 0.75, size = 1.7) +
  scale_color_manual(values = mapper_palette) +
  scale_y_log10(labels = si_labels) +
  facet_wrap(~experiment_protocol, scales = "free_y") +
  labs(
    title = "Host subtraction leaves sharply different non-host read burdens",
    subtitle = "bwa_mem consistently passed far fewer reads into the viral stage than HISAT2 or STAR",
    x = "Host mapper",
    y = "QC-passing reads not assigned to host"
  ) +
  base_theme +
  theme(legend.position = "none")

save_plot(plot_nonhost, "Figure2a.pdf", width = 11, height = 6.5)

plot_untarget <- kraken %>%
  mutate(
    host_mapper = factor(host_mapper, levels = c("bwa_mem", "HISAT2", "STAR")),
    truth_group = factor(truth_group, levels = c("negative_control", "positive_titration"), labels = truth_labels)
  ) %>%
  select(host_mapper, virus_mapper, truth_group, kraken_taxids_total, rvdb_taxids_total, scan_taxids_total) %>%
  pivot_longer(
    cols = c(kraken_taxids_total, rvdb_taxids_total, scan_taxids_total),
    names_to = "method",
    values_to = "taxids_total"
  ) %>%
  mutate(
    method = recode(
      method,
      kraken_taxids_total = "Kraken2",
      rvdb_taxids_total = "RVDB read_count",
      scan_taxids_total = "ViraQuant_scan"
    ),
    method = factor(method, levels = c("Kraken2", "RVDB read_count", "ViraQuant_scan"))
  ) %>%
  ggplot(aes(x = host_mapper, y = taxids_total, fill = method, color = method)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7,
    alpha = 0.35,
    position = position_dodge(width = 0.78)
  ) +
  geom_point(
    alpha = 0.4,
    size = 1.2,
    position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.78)
  ) +
  scale_fill_manual(values = c("Kraken2" = "#6c757d", "RVDB read_count" = "#d1495b", "ViraQuant_scan" = "#00798c")) +
  scale_color_manual(values = c("Kraken2" = "#6c757d", "RVDB read_count" = "#d1495b", "ViraQuant_scan" = "#00798c")) +
  scale_y_log10(labels = si_labels) +
  facet_grid(truth_group ~ virus_mapper) +
  labs(
    title = "Host subtraction still dominates untargeted taxid inflation in scan mode",
    subtitle = "At the ncbitaxon level, ViraQuant_scan matched RVDB read_count in this benchmark while Kraken2 diverged by host mapper",
    x = "Host mapper",
    y = "Untargeted taxids detected",
    fill = NULL,
    color = NULL
  ) +
  base_theme

save_plot(plot_untarget, "Figure2b.pdf", width = 8.5, height = 6.5)

plot_target <- target %>%
  mutate(
    host_mapper = factor(host_mapper, levels = c("bwa_mem", "HISAT2", "STAR")),
    truth_group = factor(truth_group, levels = c("negative_control", "positive_titration")),
    norm_total_plot = norm_total_reads + 0.01
  ) %>%
  ggplot(aes(x = host_mapper, y = norm_total_plot, fill = truth_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.85, position = position_dodge(width = 0.75)) +
  geom_point(
    aes(color = truth_group),
    alpha = 0.55,
    size = 1.4,
    position = position_jitterdodge(dodge.width = 0.75)
  ) +
  scale_fill_manual(values = truth_palette, labels = truth_labels) +
  scale_color_manual(values = truth_palette, labels = truth_labels) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~virus_mapper) +
  labs(
    title = "Targeted viral burden still depends strongly on host subtraction strategy",
    subtitle = "Control-sample signal is lowest after bwa_mem host subtraction, especially with bwa_mem viral alignment",
    x = "Host mapper",
    y = "Targeted panel abundance (norm_total_reads + 0.01)"
  ) +
  base_theme +
  theme(legend.title = element_blank())

save_plot(plot_target, "Figure3.pdf", width = 8.5, height = 5.8)

coverage_plot_data <- coverage_case %>%
  mutate(
    method = factor(paste(host_mapper, virus_mapper, sep = " / "), levels = method_levels),
    virus = factor(virus, levels = c("HEV", "Zika")),
    abundance_plot = abundance_raw_read_count + 1
  )

plot_case_reads <- ggplot(coverage_plot_data, aes(x = method, y = abundance_plot, fill = virus)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_manual(values = c("HEV" = "#d1495b", "Zika" = "#00798c")) +
  scale_y_log10(labels = si_labels) +
  labs(
    title = "Representative sample 1_S1: mapped-read totals alone can mislead",
    subtitle = "HEV often accumulated more mapped reads than Zika under permissive subtraction",
    x = "Host mapper / virus mapper",
    y = "Mapped reads in targeted panel (+1)"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title = element_blank())

save_plot(plot_case_reads, "Figure4a.pdf", width = 8.8, height = 5)

plot_case_breadth <- ggplot(coverage_plot_data, aes(x = method, y = frac_ge_1 * 100, fill = virus)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(
    aes(
      label = ifelse(pass_k_90pct_ge_1 > 0, "90% pass", "90% fail"),
      group = virus
    ),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 2.8
  ) +
  scale_fill_manual(values = c("HEV" = "#d1495b", "Zika" = "#00798c")) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 115)) +
  labs(
    title = "Coverage breadth separates local pileups from genome-wide support",
    subtitle = "The same sample shows HEV breadth near 4% versus Zika breadth near 95-99%",
    x = "Host mapper / virus mapper",
    y = "Genome positions with depth >= 1 (%)"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title = element_blank())

save_plot(plot_case_breadth, "Figure4b.pdf", width = 8.8, height = 5)

plot_pca <- pca_scores %>%
  mutate(
    host_mapper = factor(host_mapper, levels = c("bwa_mem", "HISAT2", "STAR")),
    library_type = factor(library_type, levels = c("Control", "1:1", "10:1"))
  ) %>%
  ggplot(aes(x = PC1, y = PC2, color = experiment_protocol, shape = library_type)) +
  geom_point(size = 2.4, alpha = 0.9) +
  facet_wrap(~host_mapper, scales = "free") +
  labs(
    title = "Host-expression sample structure is preserved across host mappers",
    subtitle = "Representative count-pair PCA shows protocol-driven clustering with dose-related separation within protocols",
    x = "PC1",
    y = "PC2",
    color = "Protocol",
    shape = "Condition"
  ) +
  base_theme

save_plot(plot_pca, "Suppl_Figure1.pdf", width = 10.5, height = 5.6)
