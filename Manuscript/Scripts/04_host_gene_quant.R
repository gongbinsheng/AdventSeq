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
matrix_inventory <- list_featurecounts_matrices(args$results_dir)
summary_df <- read_featurecounts_summaries(args$results_dir, in_scope_only = TRUE)

matrix_inventory <- add_metadata(matrix_inventory, metadata)
summary_df <- add_metadata(summary_df, metadata)

included_summary <- summary_df[summary_df$include_analysis, , drop = FALSE]

gene_detected_rows <- list()
replicate_cor_rows <- list()
mapper_cor_rows <- list()
pca_rows <- list()
count_matrices <- list()

build_pca_scores <- function(counts, host_mapper, count_kind) {
  library_sizes <- colSums(counts, na.rm = TRUE)
  valid <- library_sizes > 0
  if (sum(valid) < 2) {
    return(data.frame())
  }
  counts <- counts[, valid, drop = FALSE]
  library_sizes <- library_sizes[valid]
  log_cpm <- log1p(t(t(counts) / library_sizes) * 1e6)
  pca <- stats::prcomp(t(log_cpm), center = TRUE, scale. = FALSE)
  scores <- as.data.frame(pca$x[, seq_len(min(3, ncol(pca$x))), drop = FALSE])
  scores$sample_id <- rownames(scores)
  scores$host_mapper <- host_mapper
  scores$count_kind <- count_kind
  scores
}

for (host_mapper in HOST_MAPPERS) {
  for (count_kind in c("count_pairs", "count_reads")) {
    file_info <- matrix_inventory[
      matrix_inventory$host_mapper == host_mapper &
        matrix_inventory$count_kind == count_kind,
      ,
      drop = FALSE
    ]
    if (!nrow(file_info)) {
      next
    }

    merged_counts <- merge_featurecounts_group(file_info)
    count_cols <- setdiff(names(merged_counts), c("Geneid", "Length"))
    counts <- as.matrix(merged_counts[, count_cols, drop = FALSE])
    storage.mode(counts) <- "numeric"
    rownames(counts) <- merged_counts$Geneid
    count_matrices[[paste(host_mapper, count_kind, sep = "__")]] <- counts

    write_gz_tsv(
      merged_counts,
      file.path(out_dirs$tables, sprintf("04_host_counts_%s_%s.tsv.gz", host_mapper, count_kind))
    )

    genes_detected <- data.frame(
      sample_id = colnames(counts),
      host_mapper = host_mapper,
      count_kind = count_kind,
      genes_detected = colSums(counts > 0, na.rm = TRUE),
      library_size = colSums(counts, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    gene_detected_rows[[paste(host_mapper, count_kind, sep = "__")]] <- add_metadata(genes_detected, metadata)

    scores <- build_pca_scores(counts, host_mapper, count_kind)
    if (nrow(scores)) {
      scores <- add_metadata(scores, metadata)
      pca_rows[[paste(host_mapper, count_kind, sep = "__")]] <- scores
    }

    paired_groups <- stats::na.omit(unique(metadata$sample_group))
    rep_rows <- lapply(paired_groups, function(sample_group) {
      group_meta <- metadata[metadata$sample_group == sample_group, , drop = FALSE]
      reps <- group_meta$sample_id[group_meta$technical_replicate == 1]
      reps2 <- group_meta$sample_id[group_meta$technical_replicate == 2]
      if (length(reps) != 1 || length(reps2) != 1) {
        return(NULL)
      }
      if (!(reps %in% colnames(counts) && reps2 %in% colnames(counts))) {
        return(NULL)
      }
      data.frame(
        sample_group = sample_group,
        sample_id_rep1 = reps,
        sample_id_rep2 = reps2,
        host_mapper = host_mapper,
        count_kind = count_kind,
        spearman_cor = safe_cor(log1p(counts[, reps]), log1p(counts[, reps2])),
        stringsAsFactors = FALSE
      )
    })
    rep_rows <- Filter(Negate(is.null), rep_rows)
    if (length(rep_rows)) {
      replicate_cor_rows[[paste(host_mapper, count_kind, sep = "__")]] <- do.call(rbind, rep_rows)
    }
  }
}

for (count_kind in c("count_pairs", "count_reads")) {
  mapper_pairs <- utils::combn(HOST_MAPPERS, 2, simplify = FALSE)
  for (pair in mapper_pairs) {
    counts_a <- count_matrices[[paste(pair[[1]], count_kind, sep = "__")]]
    counts_b <- count_matrices[[paste(pair[[2]], count_kind, sep = "__")]]
    if (is.null(counts_a) || is.null(counts_b)) {
      next
    }
    shared_samples <- intersect(colnames(counts_a), colnames(counts_b))
    if (!length(shared_samples)) {
      next
    }
    shared_genes <- intersect(rownames(counts_a), rownames(counts_b))
    rows <- lapply(shared_samples, function(sample_id) {
      data.frame(
        sample_id = sample_id,
        mapper_a = pair[[1]],
        mapper_b = pair[[2]],
        count_kind = count_kind,
        spearman_cor = safe_cor(
          log1p(counts_a[shared_genes, sample_id]),
          log1p(counts_b[shared_genes, sample_id])
        ),
        stringsAsFactors = FALSE
      )
    })
    mapper_cor_rows[[paste(pair[[1]], pair[[2]], count_kind, sep = "__")]] <- add_metadata(do.call(rbind, rows), metadata)
  }
}

gene_detected_df <- if (length(gene_detected_rows)) do.call(rbind, gene_detected_rows) else data.frame()
replicate_cor_df <- if (length(replicate_cor_rows)) do.call(rbind, replicate_cor_rows) else data.frame()
mapper_cor_df <- if (length(mapper_cor_rows)) do.call(rbind, mapper_cor_rows) else data.frame()
pca_df <- if (length(pca_rows)) do.call(rbind, pca_rows) else data.frame()

write_tsv(matrix_inventory, file.path(out_dirs$tables, "04_featurecounts_matrix_inventory.tsv"))
write_tsv(summary_df, file.path(out_dirs$tables, "04_featurecounts_summary.tsv"))
write_tsv(gene_detected_df, file.path(out_dirs$tables, "04_gene_detection_summary.tsv"))
write_tsv(replicate_cor_df, file.path(out_dirs$tables, "04_technical_replicate_correlations.tsv"))
write_tsv(mapper_cor_df, file.path(out_dirs$tables, "04_mapper_concordance.tsv"))
write_tsv(pca_df, file.path(out_dirs$tables, "04_host_pca_scores.tsv"))

plot_grouped_boxplot(
  included_summary$assigned_fraction,
  paste(included_summary$host_mapper, included_summary$experiment_protocol, sep = "__"),
  file.path(out_dirs$plots, "04_assigned_fraction_by_mapper_protocol.pdf"),
  main = "Assigned feature fraction by mapper/protocol",
  ylab = "Assigned fraction"
)

if (nrow(gene_detected_df)) {
  plot_grouped_boxplot(
    gene_detected_df$genes_detected[gene_detected_df$count_kind == "count_pairs"],
    gene_detected_df$host_mapper[gene_detected_df$count_kind == "count_pairs"],
    file.path(out_dirs$plots, "04_genes_detected_count_pairs_by_mapper.pdf"),
    main = "Genes detected from count-pairs matrices",
    ylab = "Detected genes"
  )
}

if (nrow(pca_df)) {
  pair_scores <- pca_df[pca_df$count_kind == "count_pairs", , drop = FALSE]
  if (nrow(pair_scores) && all(c("PC1", "PC2") %in% names(pair_scores))) {
    with_pdf(file.path(out_dirs$plots, "04_host_pca_count_pairs.pdf"), {
      graphics::plot(
        pair_scores$PC1,
        pair_scores$PC2,
        pch = 19,
        col = as.integer(as.factor(pair_scores$experiment_protocol)),
        xlab = "PC1",
        ylab = "PC2",
        main = "Host count-pair PCA across mapper/sample combinations"
      )
      graphics::legend(
        "topright",
        legend = levels(as.factor(pair_scores$experiment_protocol)),
        col = seq_along(levels(as.factor(pair_scores$experiment_protocol))),
        pch = 19,
        cex = 0.7
      )
      graphics::grid()
    })
  }
}

log_lines <- c(
  sprintf("featureCounts matrix files in scope: %d", nrow(matrix_inventory)),
  sprintf("featureCounts summary rows in scope: %d", nrow(summary_df)),
  sprintf("Gene-detection summary rows: %d", nrow(gene_detected_df)),
  sprintf("Technical replicate correlations: %d", nrow(replicate_cor_df)),
  sprintf("Mapper concordance rows: %d", nrow(mapper_cor_df))
)
write_log_lines(log_lines, file.path(out_dirs$logs, "04_host_gene_quant.log"))
