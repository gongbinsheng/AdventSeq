options(stringsAsFactors = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required but not installed.")
}

HOST_MAPPERS <- c("bwa_mem", "HISAT2", "STAR")
VIRUS_MAPPERS <- c("bwa_mem", "bowtie2")
ALL_VIRUS_MAPPERS <- c(VIRUS_MAPPERS, "minimap2")
VIRUS_REFERENCES <- c("7viruses", "RVDBv31")
VIRAQUANT_FIELDS <- c(
  "mean_depth_ge_1",
  "mapped_reads",
  "mapped_per_bp",
  "mean_depth_all",
  "median_depth_ge_1",
  "frac_ge_1"
)
COVERAGE_ABUNDANCE_FIELDS <- c(
  "mean_depth_ge_1",
  "mapped_per_bp",
  "mean_depth_all",
  "median_depth_ge_1",
  "frac_ge_1"
)

TARGETED_VIRUS_NAME_MAP <- c(
  "Feline leukemia virus Kawakami-Theilen" = "FeLV",
  "Hepatitis E virus clone p6" = "HEV",
  "Human respiratory syncytial virus" = "RSV",
  "Mammalian orthoreovirus 1 Lang" = "REO1",
  "Mammalian orthoreovirus 3 Dearing" = "REO3",
  "Porcine circovirus 1" = "PCV1",
  "Zika virus" = "Zika",
  "human gammaherpesvirus 4 (Epstein-Barr virus)" = "EBV"
)
TARGETED_VIRUS_LEVELS <- unname(TARGETED_VIRUS_NAME_MAP)

current_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(file_arg)) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
}

parse_cli_args <- function(defaults = list()) {
  args <- defaults
  tokens <- commandArgs(trailingOnly = TRUE)
  idx <- 1
  while (idx <= length(tokens)) {
    key <- tokens[[idx]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key_name <- gsub("-", "_", sub("^--", "", key))
    if (idx == length(tokens) || startsWith(tokens[[idx + 1]], "--")) {
      args[[key_name]] <- TRUE
      idx <- idx + 1
    } else {
      args[[key_name]] <- tokens[[idx + 1]]
      idx <- idx + 2
    }
  }
  args
}

resolve_common_args <- function(defaults = list()) {
  args <- parse_cli_args(defaults)
  args$results_dir <- normalizePath(args$results_dir, mustWork = FALSE)
  if (is.null(args$sample_info) || !nzchar(args$sample_info)) {
    args$sample_info <- file.path(args$results_dir, "sample_info.tsv")
  }
  args$sample_info <- normalizePath(args$sample_info, mustWork = FALSE)
  args$output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
  args
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

prepare_output_dirs <- function(output_dir) {
  out <- list(
    root = output_dir,
    tables = file.path(output_dir, "tables"),
    plots = file.path(output_dir, "plots"),
    logs = file.path(output_dir, "logs")
  )
  for (path in out) {
    ensure_dir(path)
  }
  out
}

write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

write_gz_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

write_log_lines <- function(lines, path) {
  ensure_dir(dirname(path))
  writeLines(lines, con = path, useBytes = TRUE)
}

bind_rows_fill <- function(dfs) {
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) >= 0, dfs)
  if (!length(dfs)) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(df) {
    missing_cols <- setdiff(cols, names(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    df[, cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}

with_pdf <- function(path, expr, width = 10, height = 7) {
  ensure_dir(dirname(path))
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub("%", "", x, fixed = TRUE)))
}

sanitize_metric_name <- function(x) {
  clean <- tolower(trimws(x))
  clean <- gsub("[^a-z0-9]+", "_", clean)
  clean <- gsub("^_+|_+$", "", clean)
  clean
}

combine_flags <- function(...) {
  flags <- unlist(list(...), use.names = FALSE)
  flags <- flags[!is.na(flags) & nzchar(flags) & flags != "ok"]
  if (!length(flags)) {
    return("ok")
  }
  paste(unique(flags), collapse = ";")
}

compose_normalization_status <- function(denom_total_reads, denom_non_host_reads) {
  total_flag <- ifelse(
    is.na(denom_total_reads),
    "missing_total_reads",
    ifelse(denom_total_reads <= 0, "non_positive_total_reads", "ok")
  )
  non_host_flag <- ifelse(
    is.na(denom_non_host_reads),
    "missing_non_host_denominator",
    ifelse(denom_non_host_reads <= 0, "non_positive_non_host_denominator", "ok")
  )
  mapply(combine_flags, total_flag, non_host_flag, USE.NAMES = FALSE)
}

normalize_abundance <- function(abundance_raw, denominator, scale_factor = 1e6) {
  out <- rep(NA_real_, length(abundance_raw))
  valid <- !is.na(abundance_raw) & !is.na(denominator) & denominator > 0
  out[valid] <- abundance_raw[valid] / denominator[valid] * scale_factor
  out
}

standardize_targeted_virus_name <- function(x) {
  mapped <- unname(TARGETED_VIRUS_NAME_MAP[x])
  mapped[is.na(mapped)] <- x[is.na(mapped)]
  mapped
}

read_sample_metadata <- function(sample_info_path) {
  meta <- utils::read.delim(sample_info_path, sep = "\t", check.names = FALSE)
  required_cols <- c("library_ID", "title")
  missing_cols <- setdiff(required_cols, names(meta))
  if (length(missing_cols)) {
    stop("Missing required columns in sample info: ", paste(missing_cols, collapse = ", "))
  }

  parsed <- lapply(meta$title, function(title) {
    if (startsWith(title, "U937")) {
      return(list(
        library_type = NA_character_,
        technical_replicate = NA_integer_,
        experiment_protocol = NA_character_,
        truth_group = "excluded_u937"
      ))
    }
    parts <- strsplit(title, "-", fixed = TRUE)[[1]]
    list(
      library_type = parts[[1]],
      technical_replicate = as.integer(parts[[2]]),
      experiment_protocol = paste(parts[-c(1, 2)], collapse = "-"),
      truth_group = if (parts[[1]] == "Control") "negative_control" else "positive_titration"
    )
  })

  meta$sample_id <- meta$library_ID
  meta$library_type <- vapply(parsed, `[[`, character(1), "library_type")
  meta$technical_replicate <- vapply(parsed, `[[`, integer(1), "technical_replicate")
  meta$experiment_protocol <- vapply(parsed, `[[`, character(1), "experiment_protocol")
  meta$truth_group <- vapply(parsed, `[[`, character(1), "truth_group")
  meta$include_analysis <- !startsWith(meta$title, "U937")
  meta$sample_group <- ifelse(
    meta$include_analysis,
    paste(meta$library_type, meta$experiment_protocol, sep = "__"),
    NA_character_
  )
  meta
}

list_fastp_files <- function(results_dir) {
  folder <- file.path(results_dir, "fastp")
  files <- list.files(folder, pattern = "\\.fastp\\.json$", full.names = TRUE)
  if (!length(files)) {
    return(data.frame())
  }
  data.frame(
    source = "fastp",
    sample_id = sub("\\.fastp\\.json$", "", basename(files)),
    path = files,
    stringsAsFactors = FALSE
  )
}

parse_named_file_info <- function(paths, pattern, fields) {
  if (!length(paths)) {
    return(data.frame())
  }
  matched <- regexec(pattern, basename(paths), perl = TRUE)
  parts <- regmatches(basename(paths), matched)
  keep <- lengths(parts) > 0
  if (!any(keep)) {
    return(data.frame())
  }
  extracted <- do.call(rbind, lapply(parts[keep], function(x) x[-1]))
  out <- as.data.frame(extracted, stringsAsFactors = FALSE)
  names(out) <- fields
  out$path <- paths[keep]
  out
}

list_stats_files <- function(results_dir) {
  files <- list.files(file.path(results_dir, "stats"), pattern = "\\.stats$", full.names = TRUE)
  out <- parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.unsorted\\.bam\\.stats$",
    c("sample_id", "host_mapper")
  )
  if (!nrow(out)) {
    return(out)
  }
  out$source <- "stats"
  out
}

list_featurecounts_matrices <- function(results_dir) {
  files <- list.files(file.path(results_dir, "featureCounts"), pattern = "\\.tsv\\.gz$", full.names = TRUE)
  out <- parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.featureCounts\\.(hg38_ncbiRefSeq_exon)\\.(count_reads|count_pairs)\\.tsv\\.gz$",
    c("sample_id", "host_mapper", "model", "count_kind")
  )
  if (!nrow(out)) {
    return(out)
  }
  out$source <- "featureCounts_matrix"
  out
}

list_featurecounts_summaries <- function(results_dir) {
  files <- list.files(file.path(results_dir, "featureCounts", "summary"), pattern = "\\.summary$", full.names = TRUE)
  out <- parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.featureCounts\\.(hg38_exon|hg38_ncbiRefSeq_exon)\\.(count_reads|count_pairs)\\.tsv\\.summary$",
    c("sample_id", "host_mapper", "model", "count_kind")
  )
  if (!nrow(out)) {
    return(out)
  }
  out$source <- "featureCounts_summary"
  out$in_scope <- out$model == "hg38_ncbiRefSeq_exon"
  out
}

list_read_count_files <- function(results_dir, in_scope_only = FALSE) {
  files <- list.files(file.path(results_dir, "read_count"), pattern = "\\.read_count\\.txt$", full.names = TRUE)
  out <- parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.unmapped2host\\.(bwa_mem|bowtie2|minimap2)\\.(7viruses|RVDBv31)\\.unsorted\\.read_count\\.txt$",
    c("sample_id", "host_mapper", "virus_mapper", "virus_reference")
  )
  if (!nrow(out)) {
    return(out)
  }
  out$source <- "read_count"
  out$in_scope <- out$virus_mapper %in% VIRUS_MAPPERS & out$virus_reference %in% VIRUS_REFERENCES
  if (in_scope_only) {
    out <- out[out$in_scope, , drop = FALSE]
  }
  out
}

list_viraquant_files <- function(results_dir, in_scope_only = FALSE) {
  files <- list.files(file.path(results_dir, "ViraQuant"), pattern = "\\.ViraQuant\\.tsv$", full.names = TRUE)
  out <- parse_named_file_info(
    files,
    "^(.*)\\.(bwa_mem|HISAT2|STAR)\\.hg38\\.unmapped2host\\.(bwa_mem|bowtie2|minimap2)\\.7viruses\\.ViraQuant\\.tsv$",
    c("sample_id", "host_mapper", "virus_mapper")
  )
  if (!nrow(out)) {
    return(out)
  }
  out$source <- "ViraQuant"
  out$virus_reference <- "7viruses"
  out$in_scope <- out$virus_mapper %in% VIRUS_MAPPERS
  if (in_scope_only) {
    out <- out[out$in_scope, , drop = FALSE]
  }
  out
}

parse_fastp_json <- function(path) {
  obj <- jsonlite::fromJSON(path)
  sample_id <- sub("\\.fastp\\.json$", "", basename(path))
  data.frame(
    sample_id = sample_id,
    path = path,
    before_filtering_total_reads = obj$summary$before_filtering$total_reads,
    after_filtering_total_reads = obj$summary$after_filtering$total_reads,
    before_filtering_total_bases = obj$summary$before_filtering$total_bases,
    after_filtering_total_bases = obj$summary$after_filtering$total_bases,
    after_filtering_q20_rate = obj$summary$after_filtering$q20_rate,
    after_filtering_q30_rate = obj$summary$after_filtering$q30_rate,
    after_filtering_gc_content = obj$summary$after_filtering$gc_content,
    duplication_rate = if (!is.null(obj$duplication$rate)) obj$duplication$rate else NA_real_,
    stringsAsFactors = FALSE
  )
}

read_fastp_qc <- function(results_dir) {
  files <- list_fastp_files(results_dir)
  if (!nrow(files)) {
    return(data.frame())
  }
  rows <- lapply(files$path, parse_fastp_json)
  do.call(rbind, rows)
}

parse_stats_file <- function(path) {
  info <- list_stats_files(dirname(dirname(path)))
  matched <- info[info$path == path, , drop = FALSE]
  sample_id <- if (nrow(matched)) matched$sample_id[[1]] else sub("\\..*$", "", basename(path))
  host_mapper <- if (nrow(matched)) matched$host_mapper[[1]] else NA_character_

  lines <- readLines(path, warn = FALSE)
  current_section <- NULL
  view_metrics <- list()
  flag_metrics <- list()

  for (line in lines) {
    if (!nzchar(line) || line == "====================") {
      next
    }
    if (line == "# samtools view") {
      current_section <- "view"
      next
    }
    if (line == "# samtools flagstat") {
      current_section <- "flagstat"
      next
    }
    if (startsWith(line, "# samtools stats")) {
      current_section <- "ignore"
      next
    }
    if (startsWith(line, "#")) {
      next
    }

    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (current_section == "view" && length(fields) >= 2) {
      view_metrics[[sanitize_metric_name(fields[[1]])]] <- suppressWarnings(as.numeric(fields[[2]]))
    }
    if (current_section == "flagstat" && length(fields) >= 3) {
      flag_metrics[[sanitize_metric_name(fields[[3]])]] <- safe_numeric(fields[[1]])
    }
  }

  get_flag <- function(name) {
    value <- flag_metrics[[name]]
    if (is.null(value)) NA_real_ else value
  }
  get_view <- function(name) {
    value <- view_metrics[[name]]
    if (is.null(value)) NA_real_ else value
  }

  data.frame(
    sample_id = sample_id,
    host_mapper = host_mapper,
    path = path,
    total_reads = get_flag("total_qc_passed_reads_qc_failed_reads"),
    mapped_reads = get_flag("mapped"),
    primary_mapped_reads = get_flag("primary_mapped"),
    paired_in_sequencing = get_flag("paired_in_sequencing"),
    properly_paired_reads = get_flag("properly_paired"),
    with_itself_and_mate_mapped = get_flag("with_itself_and_mate_mapped"),
    singletons = get_flag("singletons"),
    mate_mapped_different_chr = get_flag("with_mate_mapped_to_a_different_chr"),
    mate_mapped_different_chr_mapq5 = get_flag("with_mate_mapped_to_a_different_chr_map_q_5"),
    view_primary_mapped_pairs = get_view("primary_mapped_pairs"),
    view_primary_mapped_properly_paired = get_view("primary_mapped_properly_paired"),
    stringsAsFactors = FALSE
  )
}

read_stats_qc <- function(results_dir) {
  files <- list_stats_files(results_dir)
  if (!nrow(files)) {
    return(data.frame())
  }
  rows <- lapply(files$path, parse_stats_file)
  out <- do.call(rbind, rows)
  out$mapped_rate <- out$mapped_reads / out$total_reads
  out$primary_mapped_rate <- out$primary_mapped_reads / out$total_reads
  out$properly_paired_rate <- out$properly_paired_reads / out$paired_in_sequencing
  out$singleton_rate <- out$singletons / out$paired_in_sequencing
  out
}

build_normalization_denominators <- function(results_dir) {
  fastp <- read_fastp_qc(results_dir)
  stats <- read_stats_qc(results_dir)
  if (!nrow(fastp)) {
    return(data.frame())
  }

  base_grid <- expand.grid(
    sample_id = unique(fastp$sample_id),
    host_mapper = HOST_MAPPERS,
    stringsAsFactors = FALSE
  )
  norm_df <- merge(
    base_grid,
    fastp[, c("sample_id", "after_filtering_total_reads")],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  norm_df <- merge(
    norm_df,
    stats[, c("sample_id", "host_mapper", "with_itself_and_mate_mapped")],
    by = c("sample_id", "host_mapper"),
    all.x = TRUE,
    sort = FALSE
  )
  norm_df$denom_total_reads <- norm_df$after_filtering_total_reads
  norm_df$denom_non_host_reads <- norm_df$after_filtering_total_reads - norm_df$with_itself_and_mate_mapped
  norm_df$normalization_status <- compose_normalization_status(
    norm_df$denom_total_reads,
    norm_df$denom_non_host_reads
  )
  norm_df
}

parse_featurecounts_summary_file <- function(path) {
  info <- list_featurecounts_summaries(dirname(dirname(dirname(path))))
  matched <- info[info$path == path, , drop = FALSE]
  tab <- utils::read.delim(path, sep = "\t", check.names = FALSE)
  status_values <- setNames(tab[[2]], sanitize_metric_name(tab$Status))
  get_status <- function(name) {
    value <- status_values[[name]]
    if (is.null(value)) NA_real_ else suppressWarnings(as.numeric(value))
  }

  data.frame(
    sample_id = matched$sample_id[[1]],
    host_mapper = matched$host_mapper[[1]],
    model = matched$model[[1]],
    count_kind = matched$count_kind[[1]],
    path = path,
    assigned = get_status("assigned"),
    unassigned_unmapped = get_status("unassigned_unmapped"),
    unassigned_singleton = get_status("unassigned_singleton"),
    unassigned_chimera = get_status("unassigned_chimera"),
    unassigned_multimapping = get_status("unassigned_multimapping"),
    unassigned_nofeatures = get_status("unassigned_nofeatures"),
    unassigned_ambiguity = get_status("unassigned_ambiguity"),
    stringsAsFactors = FALSE
  )
}

read_featurecounts_summaries <- function(results_dir, in_scope_only = TRUE) {
  files <- list_featurecounts_summaries(results_dir)
  if (!nrow(files)) {
    return(data.frame())
  }
  if (in_scope_only) {
    files <- files[files$in_scope, , drop = FALSE]
  }
  rows <- lapply(files$path, parse_featurecounts_summary_file)
  out <- do.call(rbind, rows)
  out$assigned_fraction <- with(out, assigned / rowSums(out[, c(
    "assigned",
    "unassigned_unmapped",
    "unassigned_singleton",
    "unassigned_chimera",
    "unassigned_multimapping",
    "unassigned_nofeatures",
    "unassigned_ambiguity"
  )], na.rm = TRUE))
  out
}

read_featurecounts_matrix_file <- function(path) {
  info <- list_featurecounts_matrices(dirname(dirname(path)))
  matched <- info[info$path == path, , drop = FALSE]
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, sep = "\t", check.names = FALSE, comment.char = "#")
}

extract_featurecounts_counts <- function(path) {
  info <- list_featurecounts_matrices(dirname(dirname(path)))
  matched <- info[info$path == path, , drop = FALSE]
  tab <- read_featurecounts_matrix_file(path)
  out <- data.frame(
    Geneid = tab$Geneid,
    Length = tab$Length,
    count = tab[[ncol(tab)]],
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "count"] <- matched$sample_id[[1]]
  out
}

merge_featurecounts_group <- function(file_info) {
  count_tables <- lapply(file_info$path, extract_featurecounts_counts)
  merged <- Reduce(function(x, y) merge(x, y, by = c("Geneid", "Length"), all = TRUE, sort = FALSE), count_tables)
  count_cols <- setdiff(names(merged), c("Geneid", "Length"))
  merged[count_cols][is.na(merged[count_cols])] <- 0
  merged
}

parse_read_count_file <- function(path, count_method = "reads") {
  info <- list_read_count_files(dirname(dirname(path)), in_scope_only = FALSE)
  matched <- info[info$path == path, , drop = FALSE]
  tab <- utils::read.delim(path, sep = "\t", check.names = FALSE)
  names(tab) <- c("label", "count_method", "count_number")
  tab$count_number <- suppressWarnings(as.numeric(tab$count_number))

  detail <- tab[!startsWith(tab$label, "#") & tab$count_method == count_method, , drop = FALSE]
  summary_rows <- tab[startsWith(tab$label, "#"), , drop = FALSE]
  if (matched$virus_reference[[1]] == "7viruses") {
    detail$canonical_virus <- standardize_targeted_virus_name(detail$label)
  } else {
    detail$canonical_virus <- detail$label
  }

  detail$sample_id <- rep(matched$sample_id[[1]], nrow(detail))
  detail$host_mapper <- rep(matched$host_mapper[[1]], nrow(detail))
  detail$virus_mapper <- rep(matched$virus_mapper[[1]], nrow(detail))
  detail$virus_reference <- rep(matched$virus_reference[[1]], nrow(detail))
  detail$path <- rep(path, nrow(detail))

  total_virus <- summary_rows[
    summary_rows$label == "#Virus" & summary_rows$count_method == count_method,
    "count_number",
    drop = TRUE
  ]
  total_mapped <- summary_rows[
    summary_rows$label == "#Mapped" & summary_rows$count_method == count_method,
    "count_number",
    drop = TRUE
  ]
  total_unmapped <- summary_rows[
    summary_rows$label == "#Unmapped" & summary_rows$count_method == count_method,
    "count_number",
    drop = TRUE
  ]
  total_records <- summary_rows[
    summary_rows$label == "#Total" & summary_rows$count_method == count_method,
    "count_number",
    drop = TRUE
  ]

  summary_df <- data.frame(
    sample_id = matched$sample_id[[1]],
    host_mapper = matched$host_mapper[[1]],
    virus_mapper = matched$virus_mapper[[1]],
    virus_reference = matched$virus_reference[[1]],
    count_method = count_method,
    total_records = if (length(total_records)) total_records[[1]] else NA_real_,
    total_mapped = if (length(total_mapped)) total_mapped[[1]] else NA_real_,
    total_unmapped = if (length(total_unmapped)) total_unmapped[[1]] else NA_real_,
    total_virus = if (length(total_virus)) total_virus[[1]] else NA_real_,
    path = path,
    stringsAsFactors = FALSE
  )

  list(detail = detail, summary = summary_df)
}

read_read_count_data <- function(results_dir, count_method = "reads") {
  files <- list_read_count_files(results_dir, in_scope_only = TRUE)
  if (!nrow(files)) {
    return(list(detail = data.frame(), summary = data.frame()))
  }
  parsed <- lapply(files$path, parse_read_count_file, count_method = count_method)
  detail <- do.call(rbind, lapply(parsed, `[[`, "detail"))
  summary_df <- do.call(rbind, lapply(parsed, `[[`, "summary"))
  detail$count_method <- count_method
  list(detail = detail, summary = summary_df)
}

parse_viraquant_file <- function(path) {
  info <- list_viraquant_files(dirname(dirname(path)), in_scope_only = FALSE)
  matched <- info[info$path == path, , drop = FALSE]
  tab <- utils::read.delim(path, sep = "\t", check.names = FALSE)
  virus_rows <- tab[tab$level == "virus", , drop = FALSE]
  virus_rows$sample_id <- matched$sample_id[[1]]
  virus_rows$host_mapper <- matched$host_mapper[[1]]
  virus_rows$virus_mapper <- matched$virus_mapper[[1]]
  virus_rows$virus_reference <- matched$virus_reference[[1]]
  virus_rows$path <- path
  virus_rows
}

read_viraquant_data <- function(results_dir) {
  files <- list_viraquant_files(results_dir, in_scope_only = TRUE)
  if (!nrow(files)) {
    return(data.frame())
  }
  rows <- lapply(files$path, parse_viraquant_file)
  do.call(rbind, rows)
}

add_metadata <- function(df, metadata) {
  merge(df, metadata, by = "sample_id", all.x = TRUE, sort = FALSE)
}

add_normalization_columns <- function(df, normalization_df, scale_factor = 1e6, abundance_col = "abundance_raw") {
  out <- merge(
    df,
    normalization_df[, c(
      "sample_id",
      "host_mapper",
      "denom_total_reads",
      "denom_non_host_reads",
      "normalization_status"
    )],
    by = c("sample_id", "host_mapper"),
    all.x = TRUE,
    sort = FALSE
  )
  out$norm_total_reads <- normalize_abundance(out[[abundance_col]], out$denom_total_reads, scale_factor)
  out$norm_non_host_reads <- normalize_abundance(out[[abundance_col]], out$denom_non_host_reads, scale_factor)
  out
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

top_n_by_group <- function(df, group_cols, value_col, n = 5) {
  split_df <- split(df, interaction(df[group_cols], drop = TRUE, lex.order = TRUE))
  out <- lapply(split_df, function(part) {
    part <- part[order(part[[value_col]], decreasing = TRUE), , drop = FALSE]
    head(part, n)
  })
  do.call(rbind, out)
}

summarize_group_metrics <- function(df, group_cols, value_cols) {
  if (!nrow(df)) {
    return(data.frame())
  }
  parts <- split(df, interaction(df[group_cols], drop = TRUE, lex.order = TRUE))
  rows <- lapply(parts, function(part) {
    row <- as.list(part[1, group_cols, drop = FALSE])
    for (value_col in value_cols) {
      values <- suppressWarnings(as.numeric(part[[value_col]]))
      row[[paste0(value_col, "_median")]] <- if (all(is.na(values))) NA_real_ else stats::median(values, na.rm = TRUE)
      row[[paste0(value_col, "_mean")]] <- if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
    }
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  bind_rows_fill(rows)
}

plot_grouped_boxplot <- function(values, groups, path, main, ylab) {
  with_pdf(path, {
    graphics::boxplot(values ~ groups, las = 2, main = main, ylab = ylab, col = "grey85")
    graphics::grid()
  })
}

plot_simple_bar <- function(height, names_vec, path, main, ylab, col = "steelblue") {
  with_pdf(path, {
    graphics::barplot(height, names.arg = names_vec, las = 2, col = col, main = main, ylab = ylab)
    graphics::grid()
  })
}
