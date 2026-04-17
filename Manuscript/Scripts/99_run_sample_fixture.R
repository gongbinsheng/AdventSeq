#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[[1]]), mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

args <- parse_cli_args(list(
  source_results_dir = "Results",
  sample_results_dir = file.path("sample_data", "Results"),
  output_dir = file.path("Analyses", "sample_test"),
  sample_info = "",
  count_method = "reads",
  abundance_field = "mean_depth_ge_1",
  normalization_mode = "both"
))
args$source_results_dir <- normalizePath(args$source_results_dir, mustWork = FALSE)
args$sample_results_dir <- normalizePath(args$sample_results_dir, mustWork = FALSE)
args$output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
if (is.null(args$sample_info) || !nzchar(args$sample_info)) {
  args$sample_info <- file.path(args$sample_results_dir, "sample_info.tsv")
}
args$sample_info <- normalizePath(args$sample_info, mustWork = FALSE)

ensure_dir(args$output_dir)
ensure_dir(file.path(args$output_dir, "logs"))

run_script <- function(script_name, extra_args) {
  cmd <- c(file.path(script_dir, script_name), extra_args)
  output <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (status != 0L) {
    stop(
      sprintf("Script failed: %s\n%s", script_name, paste(output, collapse = "\n")),
      call. = FALSE
    )
  }
  output
}

driver_log <- character(0)

driver_log <- c(driver_log, "Running 90_make_sample_data.R")
driver_log <- c(
  driver_log,
  run_script("90_make_sample_data.R", c(
    "--results-dir", args$source_results_dir,
    "--sample-info", file.path(args$source_results_dir, "sample_info.tsv"),
    "--output-results-dir", args$sample_results_dir
  ))
)

stage_args <- c(
  "--results-dir", args$sample_results_dir,
  "--sample-info", args$sample_info,
  "--output-dir", args$output_dir
)

driver_log <- c(driver_log, "Running 01_metadata_inventory.R")
driver_log <- c(driver_log, run_script("01_metadata_inventory.R", stage_args))
driver_log <- c(driver_log, "Running 02_fastp_qc.R")
driver_log <- c(driver_log, run_script("02_fastp_qc.R", stage_args))
driver_log <- c(driver_log, "Running 03_host_mapping_qc.R")
driver_log <- c(driver_log, run_script("03_host_mapping_qc.R", stage_args))
driver_log <- c(driver_log, "Running 04_host_gene_quant.R")
driver_log <- c(driver_log, run_script("04_host_gene_quant.R", stage_args))
driver_log <- c(driver_log, "Running 05_virus_read_count.R")
driver_log <- c(driver_log, run_script("05_virus_read_count.R", c(
  stage_args,
  "--count-method", args$count_method
)))
driver_log <- c(driver_log, "Running 06_viraquant.R")
driver_log <- c(driver_log, run_script("06_viraquant.R", c(
  stage_args,
  "--abundance-field", args$abundance_field,
  "--normalization-mode", args$normalization_mode
)))
driver_log <- c(driver_log, "Running 07_integrate_compare.R")
driver_log <- c(driver_log, run_script("07_integrate_compare.R", stage_args))

write_log_lines(driver_log, file.path(args$output_dir, "logs", "99_run_sample_fixture_driver.log"))
