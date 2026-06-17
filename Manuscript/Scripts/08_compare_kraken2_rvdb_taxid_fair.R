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
  virus_mapper = "",
  kraken_background = "viral.inspect.with_rvdb.tsv.gz",
  background_ranks = "S,S1,S2",
  min_rpm = "1"
))

if (!(args$count_method %in% c("reads", "pairs"))) {
  stop("--count-method must be one of: reads, pairs")
}
if (!nzchar(args$kraken_background) || !file.exists(args$kraken_background)) {
  stop("--kraken-background must point to viral.inspect.with_rvdb.tsv.gz (or .tsv)")
}

normalize_portable_path <- function(path, mustWork = FALSE) {
  normalizePath(path, winslash = "/", mustWork = mustWork)
}

args$results_dir <- normalize_portable_path(args$results_dir, mustWork = FALSE)
args$sample_info <- normalize_portable_path(args$sample_info, mustWork = FALSE)
args$output_dir <- normalize_portable_path(args$output_dir, mustWork = FALSE)
args$kraken_background <- normalize_portable_path(args$kraken_background, mustWork = FALSE)


parse_rank_arg <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals <- vals[nzchar(vals)]
  unique(vals)
}
selected_ranks <- parse_rank_arg(args$background_ranks)
if (!length(selected_ranks)) {
  stop("--background-ranks must contain at least one rank code, e.g. S,S1,S2")
}

metadata <- read_sample_metadata(args$sample_info)
out_dirs <- prepare_output_dirs(args$output_dir)

# --- Library-size (RPM) normalization for the database-search methods ---
# To make the cross-method comparison fair, Kraken2 and RVDB read_count detections
# are normalized by library size (total reads after fastp filtering, which equals
# denom_total_reads) and required to reach a minimum reads-per-million (RPM)
# threshold, instead of counting any single mapped read as a detection.
# ViraQuant_scan keeps its breadth-aware coverage criterion (depth_at_50pct >= 1,
# applied at parse time) because depth and read counts are different quantities.
min_rpm <- suppressWarnings(as.numeric(args$min_rpm))
if (is.na(min_rpm)) min_rpm <- 1
norm_denoms <- build_normalization_denominators(args$results_dir)
sample_denom <- if (is.data.frame(norm_denoms) && nrow(norm_denoms)) {
  unique(norm_denoms[, c("sample_id", "denom_total_reads")])
} else {
  data.frame(sample_id = character(), denom_total_reads = numeric(), stringsAsFactors = FALSE)
}
lookup_denom_total <- function(sid) {
  if (!nrow(sample_denom)) return(NA_real_)
  sample_denom$denom_total_reads[match(sid, sample_denom$sample_id)]
}


read_tsv_maybe_gz <- function(path, ...) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

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
  scan_dir <- file.path(results_dir, "ViraQuant_scan", "ncbitaxon")
  if (!dir.exists(scan_dir)) {
    return(data.frame(path = character(), sample_id = character(), host_mapper = character(), virus_mapper = character(), stringsAsFactors = FALSE))
  }
  files <- list.files(scan_dir, pattern = "\\.RVDBv31\\.ViraQuant_scan\\.tsv\\.gz$", full.names = TRUE)
  if (!length(files)) {
    return(data.frame(path = character(), sample_id = character(), host_mapper = character(), virus_mapper = character(), stringsAsFactors = FALSE))
  }
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
  tab <- utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, fill = TRUE, stringsAsFactors = FALSE)
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
    taxid = trimws(tab$label[main_rows]),
    rvdb_count = tab$count_number[main_rows],
    stringsAsFactors = FALSE
  )
  main_df$rvdb_label <- main_df$taxid

  fallback_df <- data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    virus_mapper = virus_mapper,
    fallback_name = trimws(tab$label[fallback_rows]),
    rvdb_count = tab$count_number[fallback_rows],
    stringsAsFactors = FALSE
  )

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
  tab <- utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab)) {
    return(data.frame())
  }

  keep <- tab$level == "ncbitaxon"
  # ViraQuant_scan detection criterion: require breadth-aware support via
  # depth_at_50pct >= 1, i.e. at least half of the reference is covered at >= 1x.
  # The previous criterion (mean_depth_ge_1 >= 1) was effectively a non-filter,
  # because mean depth over covered positions is >= 1 whenever any position has
  # coverage, so it passed essentially every taxid with a single mapped read and
  # inflated the ViraQuant_scan counts. depth_at_50pct >= 1 removes low-breadth
  # spurious hits.
  if ("depth_at_50pct" %in% names(tab)) {
    tab$depth_at_50pct <- suppressWarnings(as.numeric(tab$depth_at_50pct))
    keep <- keep & !is.na(tab$depth_at_50pct) & tab$depth_at_50pct >= 1
  } else if ("mean_depth_ge_1" %in% names(tab)) {
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
  if (isTRUE(kraken_present)) pieces <- c(pieces, "kraken2")
  if (isTRUE(rvdb_present)) pieces <- c(pieces, "rvdb_read_count")
  if (isTRUE(scan_present)) pieces <- c(pieces, "viraquant_scan")
  if (!length(pieces)) return("none")
  paste(pieces, collapse = ";")
}

first_non_na <- function(x) {
  idx <- which(!is.na(x) & x != "")
  if (!length(idx)) return(NA)
  x[[idx[[1]]]]
}

safe_num_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_num_sum <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  sum(x)
}

collapse_unique_rows <- function(df, by_cols, agg_map) {
  if (!nrow(df)) return(df)
  key <- do.call(paste, c(df[by_cols], sep = "\r"))
  split_idx <- split(seq_len(nrow(df)), key)
  rows <- lapply(split_idx, function(idx) {
    sub <- df[idx, , drop = FALSE]
    out <- sub[1, by_cols, drop = FALSE]
    for (nm in names(agg_map)) {
      out[[nm]] <- agg_map[[nm]](sub[[nm]])
    }
    out
  })
  bind_rows_fill(rows)
}

parse_background_table <- function(path, selected_ranks) {
  bg <- read_tsv_maybe_gz(path)
  required_cols <- c("rank", "taxid", "in_rvdb")
  missing_cols <- setdiff(required_cols, names(bg))
  if (length(missing_cols)) {
    stop(sprintf("Background file missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  bg$rank <- trimws(as.character(bg$rank))
  bg$taxid <- trimws(as.character(bg$taxid))
  bg$in_rvdb <- trimws(as.character(bg$in_rvdb))
  bg$rvdb_accessions <- if ("rvdb_accessions" %in% names(bg)) as.character(bg$rvdb_accessions) else NA_character_
  bg <- bg[!is.na(bg$taxid) & nzchar(bg$taxid) & bg$taxid != "0", , drop = FALSE]
  bg <- bg[bg$in_rvdb %in% c("1", "TRUE", "True", "true"), , drop = FALSE]
  bg <- bg[bg$rank %in% selected_ranks, , drop = FALSE]
  bg <- unique(data.frame(
    taxid = bg$taxid,
    background_rank = bg$rank,
    background_name = if ("name" %in% names(bg)) as.character(bg$name) else NA_character_,
    background_rvdb_accessions = bg$rvdb_accessions,
    stringsAsFactors = FALSE
  ))
  bg[order(match(bg$background_rank, selected_ranks), bg$taxid), , drop = FALSE]
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

if (!nrow(kraken_files)) stop("No Kraken2 files matched the requested filters.")
if (!nrow(rvdb_files)) stop("No RVDBv31 read_count files matched the requested filters.")

background_table <- parse_background_table(args$kraken_background, selected_ranks)
if (!nrow(background_table)) {
  stop("No background taxids remained after filtering in_rvdb and selected ranks.")
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
rvdb_main$taxid <- ifelse(rvdb_main$taxid == "", NA_character_, rvdb_main$taxid)
rvdb_main <- rvdb_main[!is.na(rvdb_main$taxid) & rvdb_main$taxid != "0", , drop = FALSE]

scan_rows <- if (nrow(scan_files)) {
  bind_rows_fill(lapply(seq_len(nrow(scan_files)), function(i) {
    parse_viraquant_scan_file(
      scan_files$path[[i]],
      scan_files$sample_id[[i]],
      scan_files$host_mapper[[i]],
      scan_files$virus_mapper[[i]]
    )
  }))
} else {
  data.frame(
    sample_id = character(),
    host_mapper = character(),
    virus_mapper = character(),
    taxid = character(),
    scan_taxid_name = character(),
    scan_n_contigs = numeric(),
    scan_length = numeric(),
    scan_mapped_reads = numeric(),
    scan_mapped_per_bp = numeric(),
    scan_breadth_cov_gt0 = numeric(),
    scan_mean_depth_all = numeric(),
    scan_frac_ge_1 = numeric(),
    scan_mean_depth_ge_1 = numeric(),
    scan_median_depth_ge_1 = numeric(),
    scan_best_contig_id = character(),
    stringsAsFactors = FALSE
  )
}
scan_rows <- add_metadata(scan_rows, metadata)
if (nrow(scan_rows)) {
  scan_rows$taxid <- ifelse(scan_rows$taxid == "", NA_character_, scan_rows$taxid)
  scan_rows <- scan_rows[!is.na(scan_rows$taxid) & scan_rows$taxid != "0", , drop = FALSE]
}

kraken_rows_fair <- merge(kraken_rows, background_table, by = "taxid", all = FALSE, sort = FALSE)
rvdb_main_fair <- merge(rvdb_main, background_table, by = "taxid", all = FALSE, sort = FALSE)
scan_rows_fair <- if (nrow(scan_rows)) merge(scan_rows, background_table, by = "taxid", all = FALSE, sort = FALSE) else data.frame()

kraken_rows_fair <- collapse_unique_rows(
  kraken_rows_fair,
  by_cols = c("sample_id", "host_mapper", "taxid", "background_rank"),
  agg_map = list(
    kraken_percent = safe_num_max,
    kraken_clade_count = safe_num_max,
    kraken_direct_count = safe_num_max,
    kraken_rank_code = first_non_na,
    kraken_name = first_non_na,
    background_name = first_non_na,
    background_rvdb_accessions = first_non_na
  )
)
rvdb_main_fair <- collapse_unique_rows(
  rvdb_main_fair,
  by_cols = c("sample_id", "host_mapper", "virus_mapper", "taxid", "background_rank"),
  agg_map = list(
    rvdb_count = safe_num_sum,
    rvdb_label = first_non_na,
    background_name = first_non_na,
    background_rvdb_accessions = first_non_na
  )
)
scan_rows_fair <- if (nrow(scan_rows_fair)) {
  collapse_unique_rows(
    scan_rows_fair,
    by_cols = c("sample_id", "host_mapper", "virus_mapper", "taxid", "background_rank"),
    agg_map = list(
      scan_taxid_name = first_non_na,
      scan_n_contigs = safe_num_max,
      scan_length = safe_num_max,
      scan_mapped_reads = safe_num_max,
      scan_mapped_per_bp = safe_num_max,
      scan_breadth_cov_gt0 = safe_num_max,
      scan_mean_depth_all = safe_num_max,
      scan_frac_ge_1 = safe_num_max,
      scan_mean_depth_ge_1 = safe_num_max,
      scan_median_depth_ge_1 = safe_num_max,
      scan_best_contig_id = first_non_na,
      background_name = first_non_na,
      background_rvdb_accessions = first_non_na
    )
  )
} else {
  data.frame(
    sample_id = character(),
    host_mapper = character(),
    virus_mapper = character(),
    taxid = character(),
    background_rank = character(),
    scan_taxid_name = character(),
    scan_n_contigs = numeric(),
    scan_length = numeric(),
    scan_mapped_reads = numeric(),
    scan_mapped_per_bp = numeric(),
    scan_breadth_cov_gt0 = numeric(),
    scan_mean_depth_all = numeric(),
    scan_frac_ge_1 = numeric(),
    scan_mean_depth_ge_1 = numeric(),
    scan_median_depth_ge_1 = numeric(),
    scan_best_contig_id = character(),
    background_name = character(),
    background_rvdb_accessions = character(),
    stringsAsFactors = FALSE
  )
}

comparison_groups <- unique(rbind(
  rvdb_files[, c("sample_id", "host_mapper", "virus_mapper"), drop = FALSE],
  if (nrow(scan_files)) scan_files[, c("sample_id", "host_mapper", "virus_mapper"), drop = FALSE] else data.frame(sample_id = character(), host_mapper = character(), virus_mapper = character(), stringsAsFactors = FALSE)
))

comparison_rows <- lapply(seq_len(nrow(comparison_groups)), function(i) {
  sample_id <- comparison_groups$sample_id[[i]]
  host_mapper <- comparison_groups$host_mapper[[i]]
  virus_mapper <- comparison_groups$virus_mapper[[i]]

  base_bg <- background_table
  base_bg$sample_id <- sample_id
  base_bg$host_mapper <- host_mapper
  base_bg$virus_mapper <- virus_mapper
  base_bg <- base_bg[, c("sample_id", "host_mapper", "virus_mapper", "taxid", "background_rank", "background_name", "background_rvdb_accessions"), drop = FALSE]

  kraken_subset <- kraken_rows_fair[
    kraken_rows_fair$sample_id == sample_id & kraken_rows_fair$host_mapper == host_mapper,
    c("sample_id", "host_mapper", "taxid", "background_rank", "kraken_rank_code", "kraken_name", "kraken_percent", "kraken_clade_count", "kraken_direct_count"),
    drop = FALSE
  ]
  if (nrow(kraken_subset)) {
    kraken_subset$virus_mapper <- virus_mapper
    kraken_subset <- kraken_subset[, c("sample_id", "host_mapper", "virus_mapper", "taxid", "background_rank", "kraken_rank_code", "kraken_name", "kraken_percent", "kraken_clade_count", "kraken_direct_count"), drop = FALSE]
  } else {
    kraken_subset <- data.frame(sample_id = character(), host_mapper = character(), virus_mapper = character(), taxid = character(), background_rank = character(), kraken_rank_code = character(), kraken_name = character(), kraken_percent = numeric(), kraken_clade_count = numeric(), kraken_direct_count = numeric(), stringsAsFactors = FALSE)
  }

  rvdb_subset <- rvdb_main_fair[
    rvdb_main_fair$sample_id == sample_id & rvdb_main_fair$host_mapper == host_mapper & rvdb_main_fair$virus_mapper == virus_mapper,
    c("sample_id", "host_mapper", "virus_mapper", "taxid", "background_rank", "rvdb_count", "rvdb_label"),
    drop = FALSE
  ]

  scan_subset <- if (nrow(scan_rows_fair)) {
    scan_rows_fair[
      scan_rows_fair$sample_id == sample_id & scan_rows_fair$host_mapper == host_mapper & scan_rows_fair$virus_mapper == virus_mapper,
      c("sample_id", "host_mapper", "virus_mapper", "taxid", "background_rank", "scan_taxid_name", "scan_n_contigs", "scan_length",
        "scan_mapped_reads", "scan_mapped_per_bp", "scan_breadth_cov_gt0", "scan_mean_depth_all",
        "scan_frac_ge_1", "scan_mean_depth_ge_1", "scan_median_depth_ge_1", "scan_best_contig_id"),
      drop = FALSE
    ]
  } else {
    data.frame(sample_id = character(), host_mapper = character(), virus_mapper = character(), taxid = character(), background_rank = character(), stringsAsFactors = FALSE)
  }

  joined <- Reduce(
    function(x, y) merge(x, y, by = c("sample_id", "host_mapper", "virus_mapper", "taxid", "background_rank"), all.x = TRUE, all.y = FALSE, sort = FALSE),
    list(base_bg, kraken_subset, rvdb_subset, scan_subset)
  )

  denom_total <- lookup_denom_total(sample_id)
  joined$denom_total_reads <- denom_total
  joined$kraken_rpm <- ifelse(
    !is.na(joined$kraken_clade_count) & !is.na(denom_total) & denom_total > 0,
    joined$kraken_clade_count * 1e6 / denom_total, NA_real_
  )
  joined$rvdb_rpm <- ifelse(
    !is.na(joined$rvdb_count) & !is.na(denom_total) & denom_total > 0,
    joined$rvdb_count * 1e6 / denom_total, NA_real_
  )
  # Detection: Kraken2 and RVDB require RPM >= min_rpm (library-size normalized);
  # ViraQuant_scan requires its breadth filter (already applied at parse time).
  joined$detected_by_kraken2 <- !is.na(joined$kraken_rpm) & joined$kraken_rpm >= min_rpm
  joined$detected_by_rvdb_read_count <- !is.na(joined$rvdb_rpm) & joined$rvdb_rpm >= min_rpm
  joined$detected_by_viraquant_scan <- if ("scan_mapped_reads" %in% names(joined)) !is.na(joined$scan_mapped_reads) else FALSE
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
if (!nrow(comparison) && !length(names(comparison))) {
  comparison <- data.frame(
    sample_id = character(),
    host_mapper = character(),
    virus_mapper = character(),
    taxid = character(),
    background_rank = character(),
    background_name = character(),
    background_rvdb_accessions = character(),
    kraken_rank_code = character(),
    kraken_name = character(),
    kraken_percent = numeric(),
    kraken_clade_count = numeric(),
    kraken_direct_count = numeric(),
    rvdb_count = numeric(),
    rvdb_label = character(),
    scan_taxid_name = character(),
    scan_n_contigs = numeric(),
    scan_length = numeric(),
    scan_mapped_reads = numeric(),
    scan_mapped_per_bp = numeric(),
    scan_breadth_cov_gt0 = numeric(),
    scan_mean_depth_all = numeric(),
    scan_frac_ge_1 = numeric(),
    scan_mean_depth_ge_1 = numeric(),
    scan_median_depth_ge_1 = numeric(),
    scan_best_contig_id = character(),
    detected_by_kraken2 = logical(),
    detected_by_rvdb_read_count = logical(),
    detected_by_viraquant_scan = logical(),
    detection_pattern = character(),
    stringsAsFactors = FALSE
  )
}
comparison <- add_metadata(comparison, metadata)

summarise_group_rank <- function(subset_cmp) {
  data.frame(
    background_taxids_total = nrow(subset_cmp),
    kraken_taxids_total = sum(subset_cmp$detected_by_kraken2, na.rm = TRUE),
    rvdb_taxids_total = sum(subset_cmp$detected_by_rvdb_read_count, na.rm = TRUE),
    scan_taxids_total = sum(subset_cmp$detected_by_viraquant_scan, na.rm = TRUE),
    kraken_rvdb_shared_taxids = sum(subset_cmp$detected_by_kraken2 & subset_cmp$detected_by_rvdb_read_count, na.rm = TRUE),
    kraken_only_taxids = sum(subset_cmp$detected_by_kraken2 & !subset_cmp$detected_by_rvdb_read_count, na.rm = TRUE),
    rvdb_only_taxids = sum(!subset_cmp$detected_by_kraken2 & subset_cmp$detected_by_rvdb_read_count, na.rm = TRUE),
    neither_taxids = sum(!subset_cmp$detected_by_kraken2 & !subset_cmp$detected_by_rvdb_read_count, na.rm = TRUE),
    kraken_total_clade_count = sum(subset_cmp$kraken_clade_count, na.rm = TRUE),
    rvdb_total_count = sum(subset_cmp$rvdb_count, na.rm = TRUE),
    scan_total_mapped_reads = sum(subset_cmp$scan_mapped_reads, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
idx <- 1L
for (i in seq_len(nrow(comparison_groups))) {
  sample_id <- comparison_groups$sample_id[[i]]
  host_mapper <- comparison_groups$host_mapper[[i]]
  virus_mapper <- comparison_groups$virus_mapper[[i]]
  group_cmp <- comparison[
    comparison$sample_id == sample_id & comparison$host_mapper == host_mapper & comparison$virus_mapper == virus_mapper,
    , drop = FALSE
  ]
  if (!nrow(group_cmp)) next

  for (rank_val in c(selected_ranks, "ALL_SELECTED")) {
    subset_cmp <- if (identical(rank_val, "ALL_SELECTED")) group_cmp else group_cmp[group_cmp$background_rank == rank_val, , drop = FALSE]
    if (!nrow(subset_cmp) && !identical(rank_val, "ALL_SELECTED")) next
    one <- summarise_group_rank(subset_cmp)
    one$sample_id <- sample_id
    one$host_mapper <- host_mapper
    one$virus_mapper <- virus_mapper
    one$background_rank <- rank_val
    summary_rows[[idx]] <- one
    idx <- idx + 1L
  }
}
comparison_summary <- bind_rows_fill(summary_rows)
if (!nrow(comparison_summary) && !length(names(comparison_summary))) {
  comparison_summary <- data.frame(
    sample_id = character(),
    host_mapper = character(),
    virus_mapper = character(),
    background_rank = character(),
    background_taxids_total = numeric(),
    kraken_taxids_total = numeric(),
    rvdb_taxids_total = numeric(),
    scan_taxids_total = numeric(),
    kraken_rvdb_shared_taxids = numeric(),
    kraken_only_taxids = numeric(),
    rvdb_only_taxids = numeric(),
    neither_taxids = numeric(),
    kraken_total_clade_count = numeric(),
    rvdb_total_count = numeric(),
    scan_total_mapped_reads = numeric(),
    stringsAsFactors = FALSE
  )
}
comparison_summary <- comparison_summary[, c("sample_id", "host_mapper", "virus_mapper", "background_rank", setdiff(names(comparison_summary), c("sample_id", "host_mapper", "virus_mapper", "background_rank"))), drop = FALSE]
comparison_summary <- add_metadata(comparison_summary, metadata)
comparison_summary$kraken_detection_fraction <- with(comparison_summary, ifelse(background_taxids_total > 0, kraken_taxids_total / background_taxids_total, NA_real_))
comparison_summary$rvdb_detection_fraction <- with(comparison_summary, ifelse(background_taxids_total > 0, rvdb_taxids_total / background_taxids_total, NA_real_))
comparison_summary$shared_detection_fraction <- with(comparison_summary, ifelse(background_taxids_total > 0, kraken_rvdb_shared_taxids / background_taxids_total, NA_real_))
comparison_summary$kraken_precision_vs_rvdb <- with(comparison_summary, ifelse(kraken_taxids_total > 0, kraken_rvdb_shared_taxids / kraken_taxids_total, NA_real_))
comparison_summary$rvdb_precision_vs_kraken <- with(comparison_summary, ifelse(rvdb_taxids_total > 0, kraken_rvdb_shared_taxids / rvdb_taxids_total, NA_real_))

comparison_kraken_only <- comparison[comparison$detected_by_kraken2 & !comparison$detected_by_rvdb_read_count, , drop = FALSE]
comparison_rvdb_only <- comparison[!comparison$detected_by_kraken2 & comparison$detected_by_rvdb_read_count, , drop = FALSE]
comparison_shared_kraken_rvdb <- comparison[comparison$detected_by_kraken2 & comparison$detected_by_rvdb_read_count, , drop = FALSE]
comparison_neither <- comparison[!comparison$detected_by_kraken2 & !comparison$detected_by_rvdb_read_count, , drop = FALSE]

background_summary <- data.frame(
  background_rank = c(selected_ranks, "ALL_SELECTED"),
  background_taxids_total = c(vapply(selected_ranks, function(r) sum(background_table$background_rank == r), numeric(1)), nrow(background_table)),
  stringsAsFactors = FALSE
)

write_tsv(background_table, file.path(out_dirs$tables, "08_fair_background_taxids.tsv"))
write_tsv(background_summary, file.path(out_dirs$tables, "08_fair_background_summary.tsv"))
write_tsv(kraken_rows_fair, file.path(out_dirs$tables, "08_kraken2_taxid_rows_fair_background.tsv"))
write_tsv(rvdb_main_fair, file.path(out_dirs$tables, "08_rvdbv31_read_count_taxid_rows_fair_background.tsv"))
if (nrow(scan_rows_fair)) write_tsv(scan_rows_fair, file.path(out_dirs$tables, "08_viraquant_scan_ncbitaxon_rows_fair_background.tsv"))
write_tsv(comparison, file.path(out_dirs$tables, "08_kraken2_rvdbv31_taxid_comparison_fair_background.tsv"))
write_tsv(comparison_shared_kraken_rvdb, file.path(out_dirs$tables, "08_kraken2_rvdbv31_taxid_shared_fair_background.tsv"))
write_tsv(comparison_kraken_only, file.path(out_dirs$tables, "08_kraken2_taxid_only_fair_background.tsv"))
write_tsv(comparison_rvdb_only, file.path(out_dirs$tables, "08_rvdbv31_taxid_only_fair_background.tsv"))
write_tsv(comparison_neither, file.path(out_dirs$tables, "08_neither_detected_taxid_fair_background.tsv"))
write_tsv(comparison_summary, file.path(out_dirs$tables, "08_kraken2_rvdbv31_taxid_summary_fair_background.tsv"))
write_tsv(rvdb_fallback, file.path(out_dirs$tables, "08_rvdbv31_fallback_rows.tsv"))
write_tsv(rvdb_summary, file.path(out_dirs$tables, "08_rvdbv31_summary_rows.tsv"))

log_lines <- c(
  sprintf("Kraken2 files parsed: %d", nrow(kraken_files)),
  sprintf("RVDBv31 read_count files parsed: %d", nrow(rvdb_files)),
  sprintf("ViraQuant scan files parsed: %d", nrow(scan_files)),
  sprintf("Background ranks selected: %s", paste(selected_ranks, collapse = ",")),
  sprintf("Background taxids total: %d", nrow(background_table)),
  sprintf("Kraken fair-background rows: %d", nrow(kraken_rows_fair)),
  sprintf("RVDB fair-background rows: %d", nrow(rvdb_main_fair)),
  sprintf("Comparison rows: %d", nrow(comparison)),
  sprintf("Kraken2/RVDB detection threshold: RPM >= %s of denom_total_reads (after_filtering_total_reads)", format(min_rpm)),
  "ViraQuant_scan detection: depth_at_50pct >= 1 (breadth-aware)"
)
write_log_lines(log_lines, file.path(out_dirs$logs, "08_compare_kraken2_rvdb_taxid_fair.log"))
