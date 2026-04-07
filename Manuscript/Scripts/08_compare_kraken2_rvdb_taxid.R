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

kraken_files <- list_kraken2_files(args$results_dir)
rvdb_files <- list_rvdb_files(args$results_dir)

if (nzchar(args$sample_id)) {
  kraken_files <- kraken_files[kraken_files$sample_id == args$sample_id, , drop = FALSE]
  rvdb_files <- rvdb_files[rvdb_files$sample_id == args$sample_id, , drop = FALSE]
}
if (nzchar(args$host_mapper)) {
  kraken_files <- kraken_files[kraken_files$host_mapper == args$host_mapper, , drop = FALSE]
  rvdb_files <- rvdb_files[rvdb_files$host_mapper == args$host_mapper, , drop = FALSE]
}
if (nzchar(args$virus_mapper)) {
  rvdb_files <- rvdb_files[rvdb_files$virus_mapper == args$virus_mapper, , drop = FALSE]
}

if (!nrow(kraken_files)) {
  stop("No Kraken2 files matched the requested filters.")
}
if (!nrow(rvdb_files)) {
  stop("No RVDBv31 read_count files matched the requested filters.")
}

kraken_rows <- do.call(
  rbind,
  lapply(seq_len(nrow(kraken_files)), function(i) {
    parse_kraken2_report(kraken_files$path[[i]], kraken_files$sample_id[[i]], kraken_files$host_mapper[[i]])
  })
)
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

rvdb_main <- do.call(rbind, lapply(rvdb_parsed, `[[`, "main"))
rvdb_fallback <- do.call(rbind, lapply(rvdb_parsed, `[[`, "fallback"))
rvdb_summary <- do.call(rbind, lapply(rvdb_parsed, `[[`, "summary"))
rvdb_main <- add_metadata(rvdb_main, metadata)
rvdb_fallback <- add_metadata(rvdb_fallback, metadata)
rvdb_summary <- add_metadata(rvdb_summary, metadata)

comparison_rows <- lapply(seq_len(nrow(rvdb_files)), function(i) {
  sample_id <- rvdb_files$sample_id[[i]]
  host_mapper <- rvdb_files$host_mapper[[i]]
  virus_mapper <- rvdb_files$virus_mapper[[i]]

  kraken_subset <- kraken_rows[
    kraken_rows$sample_id == sample_id & kraken_rows$host_mapper == host_mapper,
    c("sample_id", "host_mapper", "taxid", "kraken_rank_code", "kraken_name", "kraken_percent", "kraken_clade_count", "kraken_direct_count"),
    drop = FALSE
  ]
  rvdb_subset <- rvdb_main[
    rvdb_main$sample_id == sample_id &
      rvdb_main$host_mapper == host_mapper &
      rvdb_main$virus_mapper == virus_mapper,
    c("sample_id", "host_mapper", "virus_mapper", "taxid", "rvdb_count"),
    drop = FALSE
  ]

  joined <- merge(
    kraken_subset,
    rvdb_subset,
    by = c("sample_id", "host_mapper", "taxid"),
    all = TRUE,
    sort = FALSE
  )
  joined$virus_mapper[is.na(joined$virus_mapper)] <- virus_mapper
  joined$match_status <- ifelse(
    !is.na(joined$kraken_name) & !is.na(joined$rvdb_count),
    "shared_taxid",
    ifelse(!is.na(joined$kraken_name), "kraken_only", "rvdb_only")
  )
  joined
})
comparison <- do.call(rbind, comparison_rows)
comparison <- add_metadata(comparison, metadata)

summary_rows <- lapply(seq_len(nrow(rvdb_files)), function(i) {
  sample_id <- rvdb_files$sample_id[[i]]
  host_mapper <- rvdb_files$host_mapper[[i]]
  virus_mapper <- rvdb_files$virus_mapper[[i]]
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
    kraken_taxids_total = sum(!is.na(subset_cmp$kraken_name)),
    rvdb_taxids_total = sum(!is.na(subset_cmp$rvdb_count)),
    shared_taxids = sum(subset_cmp$match_status == "shared_taxid"),
    kraken_only_taxids = sum(subset_cmp$match_status == "kraken_only"),
    rvdb_only_taxids = sum(subset_cmp$match_status == "rvdb_only"),
    shared_kraken_clade_count = sum(subset_cmp$kraken_clade_count[subset_cmp$match_status == "shared_taxid"], na.rm = TRUE),
    shared_rvdb_count = sum(subset_cmp$rvdb_count[subset_cmp$match_status == "shared_taxid"], na.rm = TRUE),
    rvdb_total_count = sum(subset_cmp$rvdb_count, na.rm = TRUE),
    fallback_entries = nrow(subset_fallback),
    fallback_total_count = sum(subset_fallback$rvdb_count, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
comparison_summary <- add_metadata(do.call(rbind, summary_rows), metadata)
comparison_summary$shared_rvdb_fraction <- with(
  comparison_summary,
  ifelse(rvdb_total_count > 0, shared_rvdb_count / rvdb_total_count, NA_real_)
)

comparison_shared <- comparison[comparison$match_status == "shared_taxid", , drop = FALSE]
comparison_kraken_only <- comparison[comparison$match_status == "kraken_only", , drop = FALSE]
comparison_rvdb_only <- comparison[comparison$match_status == "rvdb_only", , drop = FALSE]

write_tsv(kraken_rows, file.path(out_dirs$tables, "08_kraken2_taxid_rows.tsv"))
write_tsv(rvdb_main, file.path(out_dirs$tables, "08_rvdbv31_taxid_rows.tsv"))
write_tsv(rvdb_fallback, file.path(out_dirs$tables, "08_rvdbv31_fallback_rows.tsv"))
write_tsv(rvdb_summary, file.path(out_dirs$tables, "08_rvdbv31_summary_rows.tsv"))
write_tsv(comparison, file.path(out_dirs$tables, "08_kraken2_vs_rvdbv31_taxid_comparison.tsv"))
write_tsv(comparison_shared, file.path(out_dirs$tables, "08_kraken2_vs_rvdbv31_taxid_shared.tsv"))
write_tsv(comparison_kraken_only, file.path(out_dirs$tables, "08_kraken2_vs_rvdbv31_taxid_kraken_only.tsv"))
write_tsv(comparison_rvdb_only, file.path(out_dirs$tables, "08_kraken2_vs_rvdbv31_taxid_rvdb_only.tsv"))
write_tsv(comparison_summary, file.path(out_dirs$tables, "08_kraken2_vs_rvdbv31_taxid_summary.tsv"))

log_lines <- c(
  sprintf("Kraken2 files parsed: %d", nrow(kraken_files)),
  sprintf("RVDBv31 files parsed: %d", nrow(rvdb_files)),
  sprintf("Kraken2 taxid rows parsed: %d", nrow(kraken_rows)),
  sprintf("RVDBv31 taxid rows parsed: %d", nrow(rvdb_main)),
  sprintf("RVDBv31 fallback rows parsed: %d", nrow(rvdb_fallback)),
  sprintf("Comparison rows written: %d", nrow(comparison)),
  sprintf("Shared-taxid rows: %d", nrow(comparison_shared)),
  sprintf("Count method used: %s", args$count_method)
)
write_log_lines(log_lines, file.path(out_dirs$logs, "08_compare_kraken2_rvdb_taxid.log"))
