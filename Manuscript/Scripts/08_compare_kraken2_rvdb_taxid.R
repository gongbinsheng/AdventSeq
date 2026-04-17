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
  sample_id = "",
  host_mapper = "",
  virus_mapper = ""
))

if (!(args$count_method %in% c("reads", "pairs"))) {
  stop("--count-method must be one of: reads, pairs")
}

metadata <- read_sample_metadata(args$sample_info)
out_dirs <- prepare_output_dirs(args$output_dir)

list_kraken2_files <- function(results_dir) {
  files <- list.files(file.path(results_dir, "Kraken2"), pattern = "\\.Kraken2_viral\\.report$", full.names = TRUE)
  parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.unmapped2host\\.Kraken2_viral\\.report$",
    c("sample_id", "host_mapper")
  )
}

list_rvdb_files <- function(results_dir) {
  files <- list.files(file.path(results_dir, "read_count"), pattern = "RVDBv31\\.unsorted\\.read_count\\.txt$", full.names = TRUE)
  parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.unmapped2host\\.(bwa_mem|bowtie2)\\.RVDBv31\\.unsorted\\.read_count\\.txt$",
    c("sample_id", "host_mapper", "virus_mapper")
  )
}

list_viraquant_scan_files <- function(results_dir) {
  files <- list.files(file.path(results_dir, "ViraQuant_scan", "ncbitaxon"), pattern = "\\.RVDBv31\\.ViraQuant_scan\\.tsv\\.gz$", full.names = TRUE)
  parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.unmapped2host\\.(bwa_mem|bowtie2)\\.RVDBv31\\.ViraQuant_scan\\.tsv\\.gz$",
    c("sample_id", "host_mapper", "virus_mapper")
  )
}

parse_kraken2_report <- function(path, sample_id, host_mapper) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  keep <- lengths(parts) >= 6
  parts <- parts[keep]
  if (!length(parts)) {
    return(data.frame())
  }
  data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    kraken_percent = vapply(parts, function(x) suppressWarnings(as.numeric(trimws(x[[1]]))), numeric(1)),
    kraken_clade_count = vapply(parts, function(x) suppressWarnings(as.numeric(trimws(x[[2]]))), numeric(1)),
    kraken_direct_count = vapply(parts, function(x) suppressWarnings(as.numeric(trimws(x[[3]]))), numeric(1)),
    kraken_rank_code = vapply(parts, function(x) trimws(x[[4]]), character(1)),
    taxid = vapply(parts, function(x) trimws(x[[5]]), character(1)),
    kraken_name = vapply(parts, function(x) trimws(x[[6]]), character(1)),
    stringsAsFactors = FALSE
  )
}

parse_rvdb_taxid_file <- function(path, sample_id, host_mapper, virus_mapper, count_method) {
  tab <- utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, fill = TRUE)
  names(tab) <- c("label", "count_method", "count_number")
  tab$count_number <- suppressWarnings(as.numeric(tab$count_number))

  section <- rep("taxid", nrow(tab))
  fallback_idx <- which(tab$label == "#Fallback to organism")
  if (length(fallback_idx)) {
    section[seq.int(from = fallback_idx[[1]], to = nrow(tab))] <- "fallback"
    section[fallback_idx[[1]]] <- "fallback_header"
  }
  tab$section <- section

  summary_rows <- startsWith(tab$label, "#") & tab$section == "taxid"
  main_rows <- !summary_rows & tab$section == "taxid" & tab$count_method == count_method
  fallback_rows <- tab$section == "fallback" & tab$count_method == count_method

  main_df <- data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    virus_mapper = virus_mapper,
    taxid = tab$label[main_rows],
    rvdb_count = tab$count_number[main_rows],
    stringsAsFactors = FALSE
  )
  main_df$rvdb_label <- main_df$taxid

  fallback_df <- data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    virus_mapper = virus_mapper,
    fallback_name = tab$label[fallback_rows],
    rvdb_count = tab$count_number[fallback_rows],
    stringsAsFactors = FALSE
  )
  fallback_df$fallback_name <- trimws(fallback_df$fallback_name)

  summary_df <- data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    virus_mapper = virus_mapper,
    label = tab$label[summary_rows & tab$count_method == count_method],
    count_number = tab$count_number[summary_rows & tab$count_method == count_method],
    stringsAsFactors = FALSE
  )

  list(main = main_df, fallback = fallback_df, summary = summary_df)
}

parse_viraquant_scan_file <- function(path, sample_id, host_mapper, virus_mapper) {
  tab <- utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
  if (!nrow(tab)) {
    return(data.frame())
  }

  keep <- tab$level == "ncbitaxon"
  if ("mean_depth_ge_1" %in% names(tab)) {
    tab$mean_depth_ge_1 <- suppressWarnings(as.numeric(tab$mean_depth_ge_1))
    keep <- keep & !is.na(tab$mean_depth_ge_1) & tab$mean_depth_ge_1 >= 1
  }

  tab <- tab[keep, , drop = FALSE]
  if (!nrow(tab)) {
    return(data.frame())
  }

  numeric_cols <- intersect(
    c("n_contigs", "length", "mapped_reads", "mapped_per_bp", "breadth_cov_gt0", "mean_depth_all",
      "frac_ge_1", "mean_depth_ge_1", "median_depth_ge_1", "best_contig_mapped_per_bp",
      "best_contig_depth_at_90pct", "best_contig_frac_ge_1"),
    names(tab)
  )
  for (col in numeric_cols) {
    tab[[col]] <- suppressWarnings(as.numeric(tab[[col]]))
  }

  data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    virus_mapper = virus_mapper,
    taxid = as.character(tab$virus),
    scan_taxid_name = if ("ncbitaxonname" %in% names(tab)) as.character(tab$ncbitaxonname) else NA_character_,
    scan_n_contigs = if ("n_contigs" %in% names(tab)) tab$n_contigs else NA_real_,
    scan_length = if ("length" %in% names(tab)) tab$length else NA_real_,
    scan_mapped_reads = if ("mapped_reads" %in% names(tab)) tab$mapped_reads else NA_real_,
    scan_mapped_per_bp = if ("mapped_per_bp" %in% names(tab)) tab$mapped_per_bp else NA_real_,
    scan_breadth_cov_gt0 = if ("breadth_cov_gt0" %in% names(tab)) tab$breadth_cov_gt0 else NA_real_,
    scan_mean_depth_all = if ("mean_depth_all" %in% names(tab)) tab$mean_depth_all else NA_real_,
    scan_frac_ge_1 = if ("frac_ge_1" %in% names(tab)) tab$frac_ge_1 else NA_real_,
    scan_mean_depth_ge_1 = if ("mean_depth_ge_1" %in% names(tab)) tab$mean_depth_ge_1 else NA_real_,
    scan_median_depth_ge_1 = if ("median_depth_ge_1" %in% names(tab)) tab$median_depth_ge_1 else NA_real_,
    scan_best_contig_id = if ("best_contig_id" %in% names(tab)) as.character(tab$best_contig_id) else NA_character_,
    stringsAsFactors = FALSE
  )
}

build_detection_pattern <- function(kraken_present, rvdb_present, scan_present) {
  pieces <- character()
  if (isTRUE(kraken_present)) {
    pieces <- c(pieces, "kraken2")
  }
  if (isTRUE(rvdb_present)) {
    pieces <- c(pieces, "rvdb_read_count")
  }
  if (isTRUE(scan_present)) {
    pieces <- c(pieces, "viraquant_scan")
  }
  if (!length(pieces)) {
    return("none")
  }
  paste(pieces, collapse = ";")
}

kraken_files <- list_kraken2_files(args$results_dir)
rvdb_files <- list_rvdb_files(args$results_dir)
scan_files <- list_viraquant_scan_files(args$results_dir)

if (nzchar(args$sample_id)) {
  kraken_files <- kraken_files[kraken_files$sample_id == args$sample_id, , drop = FALSE]
  rvdb_files <- rvdb_files[rvdb_files$sample_id == args$sample_id, , drop = FALSE]
  scan_files <- scan_files[scan_files$sample_id == args$sample_id, , drop = FALSE]
}
if (nzchar(args$host_mapper)) {
  kraken_files <- kraken_files[kraken_files$host_mapper == args$host_mapper, , drop = FALSE]
  rvdb_files <- rvdb_files[rvdb_files$host_mapper == args$host_mapper, , drop = FALSE]
  scan_files <- scan_files[scan_files$host_mapper == args$host_mapper, , drop = FALSE]
}
if (nzchar(args$virus_mapper)) {
  rvdb_files <- rvdb_files[rvdb_files$virus_mapper == args$virus_mapper, , drop = FALSE]
  scan_files <- scan_files[scan_files$virus_mapper == args$virus_mapper, , drop = FALSE]
}

if (!nrow(kraken_files)) {
  stop("No Kraken2 files matched the requested filters.")
}
if (!nrow(rvdb_files)) {
  stop("No RVDBv31 read_count files matched the requested filters.")
}
if (!nrow(scan_files)) {
  stop("No ViraQuant_scan ncbitaxon files matched the requested filters.")
}

kraken_rows <- bind_rows_fill(lapply(seq_len(nrow(kraken_files)), function(i) {
  parse_kraken2_report(kraken_files$path[[i]], kraken_files$sample_id[[i]], kraken_files$host_mapper[[i]])
}))
kraken_rows <- add_metadata(kraken_rows, metadata)
kraken_rows$taxid <- ifelse(kraken_rows$taxid == "", NA_character_, kraken_rows$taxid)
kraken_rows <- kraken_rows[!is.na(kraken_rows$taxid) & kraken_rows$taxid != "0", , drop = FALSE]

rvdb_parsed <- lapply(seq_len(nrow(rvdb_files)), function(i) {
  parse_rvdb_taxid_file(
    rvdb_files$path[[i]],
    rvdb_files$sample_id[[i]],
    rvdb_files$host_mapper[[i]],
    rvdb_files$virus_mapper[[i]],
    args$count_method
  )
})

rvdb_main <- bind_rows_fill(lapply(rvdb_parsed, `[[`, "main"))
rvdb_fallback <- bind_rows_fill(lapply(rvdb_parsed, `[[`, "fallback"))
rvdb_summary <- bind_rows_fill(lapply(rvdb_parsed, `[[`, "summary"))
rvdb_main <- add_metadata(rvdb_main, metadata)
rvdb_fallback <- add_metadata(rvdb_fallback, metadata)
rvdb_summary <- add_metadata(rvdb_summary, metadata)

scan_rows <- bind_rows_fill(lapply(seq_len(nrow(scan_files)), function(i) {
  parse_viraquant_scan_file(
    scan_files$path[[i]],
    scan_files$sample_id[[i]],
    scan_files$host_mapper[[i]],
    scan_files$virus_mapper[[i]]
  )
}))
scan_rows <- add_metadata(scan_rows, metadata)
scan_rows$taxid <- ifelse(scan_rows$taxid == "", NA_character_, scan_rows$taxid)
scan_rows <- scan_rows[!is.na(scan_rows$taxid) & scan_rows$taxid != "0", , drop = FALSE]

comparison_groups <- unique(rbind(
  rvdb_files[, c("sample_id", "host_mapper", "virus_mapper"), drop = FALSE],
  scan_files[, c("sample_id", "host_mapper", "virus_mapper"), drop = FALSE]
))

comparison_rows <- lapply(seq_len(nrow(comparison_groups)), function(i) {
  sample_id <- comparison_groups$sample_id[[i]]
  host_mapper <- comparison_groups$host_mapper[[i]]
  virus_mapper <- comparison_groups$virus_mapper[[i]]

  kraken_subset <- unique(kraken_rows[
    kraken_rows$sample_id == sample_id & kraken_rows$host_mapper == host_mapper,
    c("sample_id", "host_mapper", "taxid", "kraken_rank_code", "kraken_name", "kraken_percent", "kraken_clade_count", "kraken_direct_count"),
    drop = FALSE
  ])
  kraken_subset$virus_mapper <- virus_mapper
  kraken_subset <- kraken_subset[
    ,
    c("sample_id", "host_mapper", "virus_mapper", "taxid", "kraken_rank_code", "kraken_name", "kraken_percent", "kraken_clade_count", "kraken_direct_count"),
    drop = FALSE
  ]
  rvdb_subset <- unique(rvdb_main[
    rvdb_main$sample_id == sample_id &
      rvdb_main$host_mapper == host_mapper &
      rvdb_main$virus_mapper == virus_mapper,
    c("sample_id", "host_mapper", "virus_mapper", "taxid", "rvdb_count"),
    drop = FALSE
  ])
  scan_subset <- unique(scan_rows[
    scan_rows$sample_id == sample_id &
      scan_rows$host_mapper == host_mapper &
      scan_rows$virus_mapper == virus_mapper,
    c("sample_id", "host_mapper", "virus_mapper", "taxid", "scan_taxid_name", "scan_n_contigs", "scan_length",
      "scan_mapped_reads", "scan_mapped_per_bp", "scan_breadth_cov_gt0", "scan_mean_depth_all",
      "scan_frac_ge_1", "scan_mean_depth_ge_1", "scan_median_depth_ge_1", "scan_best_contig_id"),
    drop = FALSE
  ])

  joined <- Reduce(function(x, y) merge(x, y, by = c("sample_id", "host_mapper", "virus_mapper", "taxid"), all = TRUE, sort = FALSE), list(
    kraken_subset,
    rvdb_subset,
    scan_subset
  ))

  joined$virus_mapper[is.na(joined$virus_mapper)] <- virus_mapper
  joined$detected_by_kraken2 <- !is.na(joined$kraken_name)
  joined$detected_by_rvdb_read_count <- !is.na(joined$rvdb_count)
  joined$detected_by_viraquant_scan <- !is.na(joined$scan_mapped_reads)
  joined$detection_pattern <- vapply(
    seq_len(nrow(joined)),
    function(idx) build_detection_pattern(
      joined$detected_by_kraken2[[idx]],
      joined$detected_by_rvdb_read_count[[idx]],
      joined$detected_by_viraquant_scan[[idx]]
    ),
    character(1)
  )
  joined
})
comparison <- bind_rows_fill(comparison_rows)
comparison <- add_metadata(comparison, metadata)

summary_rows <- lapply(seq_len(nrow(comparison_groups)), function(i) {
  sample_id <- comparison_groups$sample_id[[i]]
  host_mapper <- comparison_groups$host_mapper[[i]]
  virus_mapper <- comparison_groups$virus_mapper[[i]]
  subset_cmp <- comparison[
    comparison$sample_id == sample_id &
      comparison$host_mapper == host_mapper &
      comparison$virus_mapper == virus_mapper,
    ,
    drop = FALSE
  ]
  subset_fallback <- rvdb_fallback[
    rvdb_fallback$sample_id == sample_id &
      rvdb_fallback$host_mapper == host_mapper &
      rvdb_fallback$virus_mapper == virus_mapper,
    ,
    drop = FALSE
  ]
  data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    virus_mapper = virus_mapper,
    kraken_taxids_total = sum(subset_cmp$detected_by_kraken2, na.rm = TRUE),
    rvdb_taxids_total = sum(subset_cmp$detected_by_rvdb_read_count, na.rm = TRUE),
    scan_taxids_total = sum(subset_cmp$detected_by_viraquant_scan, na.rm = TRUE),
    shared_all_three_taxids = sum(subset_cmp$detection_pattern == "kraken2;rvdb_read_count;viraquant_scan", na.rm = TRUE),
    kraken_rvdb_shared_taxids = sum(subset_cmp$detected_by_kraken2 & subset_cmp$detected_by_rvdb_read_count, na.rm = TRUE),
    kraken_scan_shared_taxids = sum(subset_cmp$detected_by_kraken2 & subset_cmp$detected_by_viraquant_scan, na.rm = TRUE),
    rvdb_scan_shared_taxids = sum(subset_cmp$detected_by_rvdb_read_count & subset_cmp$detected_by_viraquant_scan, na.rm = TRUE),
    kraken_only_taxids = sum(subset_cmp$detection_pattern == "kraken2", na.rm = TRUE),
    rvdb_only_taxids = sum(subset_cmp$detection_pattern == "rvdb_read_count", na.rm = TRUE),
    scan_only_taxids = sum(subset_cmp$detection_pattern == "viraquant_scan", na.rm = TRUE),
    kraken_rvdb_only_taxids = sum(subset_cmp$detection_pattern == "kraken2;rvdb_read_count", na.rm = TRUE),
    kraken_scan_only_taxids = sum(subset_cmp$detection_pattern == "kraken2;viraquant_scan", na.rm = TRUE),
    rvdb_scan_only_taxids = sum(subset_cmp$detection_pattern == "rvdb_read_count;viraquant_scan", na.rm = TRUE),
    shared_all_three_kraken_clade_count = sum(subset_cmp$kraken_clade_count[subset_cmp$detection_pattern == "kraken2;rvdb_read_count;viraquant_scan"], na.rm = TRUE),
    shared_all_three_rvdb_count = sum(subset_cmp$rvdb_count[subset_cmp$detection_pattern == "kraken2;rvdb_read_count;viraquant_scan"], na.rm = TRUE),
    shared_all_three_scan_mapped_reads = sum(subset_cmp$scan_mapped_reads[subset_cmp$detection_pattern == "kraken2;rvdb_read_count;viraquant_scan"], na.rm = TRUE),
    rvdb_total_count = sum(subset_cmp$rvdb_count, na.rm = TRUE),
    scan_total_mapped_reads = sum(subset_cmp$scan_mapped_reads, na.rm = TRUE),
    scan_total_n_contigs = sum(subset_cmp$scan_n_contigs, na.rm = TRUE),
    fallback_entries = nrow(subset_fallback),
    fallback_total_count = sum(subset_fallback$rvdb_count, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
comparison_summary <- add_metadata(bind_rows_fill(summary_rows), metadata)
comparison_summary$shared_all_three_rvdb_fraction <- with(
  comparison_summary,
  ifelse(rvdb_total_count > 0, shared_all_three_rvdb_count / rvdb_total_count, NA_real_)
)
comparison_summary$shared_all_three_scan_fraction <- with(
  comparison_summary,
  ifelse(scan_total_mapped_reads > 0, shared_all_three_scan_mapped_reads / scan_total_mapped_reads, NA_real_)
)

comparison_shared_all_three <- comparison[comparison$detection_pattern == "kraken2;rvdb_read_count;viraquant_scan", , drop = FALSE]
comparison_kraken_only <- comparison[comparison$detection_pattern == "kraken2", , drop = FALSE]
comparison_rvdb_only <- comparison[comparison$detection_pattern == "rvdb_read_count", , drop = FALSE]
comparison_scan_only <- comparison[comparison$detection_pattern == "viraquant_scan", , drop = FALSE]

write_tsv(kraken_rows, file.path(out_dirs$tables, "08_kraken2_taxid_rows.tsv"))
write_tsv(rvdb_main, file.path(out_dirs$tables, "08_rvdbv31_read_count_taxid_rows.tsv"))
write_tsv(rvdb_fallback, file.path(out_dirs$tables, "08_rvdbv31_fallback_rows.tsv"))
write_tsv(rvdb_summary, file.path(out_dirs$tables, "08_rvdbv31_summary_rows.tsv"))
write_tsv(scan_rows, file.path(out_dirs$tables, "08_viraquant_scan_ncbitaxon_rows.tsv"))
write_tsv(comparison, file.path(out_dirs$tables, "08_kraken2_rvdbv31_viraquant_scan_taxid_comparison.tsv"))
write_tsv(comparison_shared_all_three, file.path(out_dirs$tables, "08_kraken2_rvdbv31_viraquant_scan_taxid_shared_all_three.tsv"))
write_tsv(comparison_kraken_only, file.path(out_dirs$tables, "08_kraken2_rvdbv31_viraquant_scan_taxid_kraken_only.tsv"))
write_tsv(comparison_rvdb_only, file.path(out_dirs$tables, "08_kraken2_rvdbv31_viraquant_scan_taxid_rvdb_only.tsv"))
write_tsv(comparison_scan_only, file.path(out_dirs$tables, "08_kraken2_rvdbv31_viraquant_scan_taxid_scan_only.tsv"))
write_tsv(comparison_summary, file.path(out_dirs$tables, "08_kraken2_rvdbv31_viraquant_scan_taxid_summary.tsv"))

log_lines <- c(
  sprintf("Kraken2 files parsed: %d", nrow(kraken_files)),
  sprintf("RVDBv31 read_count files parsed: %d", nrow(rvdb_files)),
  sprintf("ViraQuant_scan ncbitaxon files parsed: %d", nrow(scan_files)),
  sprintf("Kraken2 taxid rows parsed: %d", nrow(kraken_rows)),
  sprintf("RVDBv31 read_count taxid rows parsed: %d", nrow(rvdb_main)),
  sprintf("ViraQuant_scan ncbitaxon rows parsed: %d", nrow(scan_rows)),
  sprintf("RVDBv31 fallback rows parsed: %d", nrow(rvdb_fallback)),
  sprintf("Comparison rows written: %d", nrow(comparison)),
  sprintf("Shared-all-three taxid rows: %d", nrow(comparison_shared_all_three)),
  sprintf("Count method used: %s", args$count_method)
)
write_log_lines(log_lines, file.path(out_dirs$logs, "08_compare_kraken2_rvdb_taxid.log"))
