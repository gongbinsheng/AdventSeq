# Scripts Overview

This directory contains the main R analysis pipeline, helper utilities, fixture-generation utilities, and one standalone Kraken2/RVDBv31 taxid comparison script.

## Core Helpers

- `helpers.R`
  - Shared utilities for CLI parsing, metadata parsing, file discovery, table writing, plotting, normalization helpers, and parsers used by the main pipeline.

## Main Analysis Pipeline

- `01_metadata_inventory.R`
  - Parses `sample_info.tsv`, excludes `U937*`, validates the expected sample design, inventories prompt-scope files, and writes completeness tables.

- `02_fastp_qc.R`
  - Extracts `fastp` QC metrics and writes normalization denominators based on `after_filtering.total_reads` and host-mapper-specific non-host estimates.

- `03_host_mapping_qc.R`
  - Parses the first two sections of `samtools` stats/flagstat outputs and summarizes host-mapping metrics by mapper and protocol.

- `04_host_gene_quant.R`
  - Builds host gene-count matrices from `featureCounts` `hg38_ncbiRefSeq_exon` outputs and reports detection, replicate concordance, mapper concordance, and PCA summaries.

- `05_virus_read_count.R`
  - Parses prompt-scope `read_count` files, writes raw and normalized abundance tables, and summarizes targeted and untargeted virus signals.
  - Key args:
    - `--count-method reads|pairs`
    - `--scale-factor 1000000`

- `06_viraquant.R`
  - Parses `ViraQuant` virus-level rows, applies the same normalization denominators as `05`, and supports user-selectable abundance fields.
  - Key args:
    - `--abundance-field mean_depth_ge_1|mapped_reads|mapped_per_bp|mean_depth_all|median_depth_ge_1|frac_ge_1`
    - `--normalization-mode both|total_reads|non_host_reads`

- `07_integrate_compare.R`
  - Integrates host, `read_count`, and `ViraQuant` outputs into matched comparison tables, detection summaries, and host-versus-virus tradeoff summaries.

## Kraken2 / RVDBv31 Taxid Comparison

- `08_compare_kraken2_rvdb_taxid.R`
  - Standalone comparison script for Kraken2 viral reports versus RVDBv31 `read_count` results using taxonomy IDs.
  - Kraken2 uses taxid from report column 5.
  - RVDBv31 main taxid rows are compared directly by `sample_id + host_mapper + virus_mapper + taxid`.
  - RVDBv31 rows below `#Fallback to organism` are parsed separately and reported as fallback rows because they do not have taxids for direct joining.
  - Outputs:
    - parsed Kraken2 taxid rows
    - parsed RVDBv31 taxid rows
    - parsed RVDBv31 fallback rows
    - shared / Kraken-only / RVDB-only comparison tables
    - per-combination summary table
  - Key args:
    - `--results-dir`
    - `--sample-info`
    - `--output-dir`
    - `--count-method reads|pairs`
    - optional filters: `--sample-id`, `--host-mapper`, `--virus-mapper`

## Fixture / Smoke-Test Utilities

- `90_make_sample_data.R`
  - Builds the representative subset in `sample_data/Results` and truncates large `featureCounts` matrices for fast smoke testing.

- `99_run_sample_fixture.R`
  - Refreshes `sample_data/Results` and runs stages `01` through `07` only on the sample fixture.

## Notes

- The main analysis pipeline currently runs through `07`.
- `08_compare_kraken2_rvdb_taxid.R` is standalone and does not change or extend the existing `01`-`07` workflow.
- Most pipeline outputs are written under `Analyses/<run_name>/{tables,plots,logs}`.

## Usage

Run commands from the manuscript root directory:

```bash
cd /Users/bgong/workspace/AdVentSeq/Manuscript
```

### 1. Run the sample fixture smoke test

This rebuilds `sample_data/Results` and runs stages `01` to `07` on the representative subset.

```bash
Rscript Scripts/99_run_sample_fixture.R
```

### 2. Run the main pipeline stage by stage on full results

```bash
Rscript Scripts/01_metadata_inventory.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run
Rscript Scripts/02_fastp_qc.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run
Rscript Scripts/03_host_mapping_qc.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run
Rscript Scripts/04_host_gene_quant.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run
Rscript Scripts/05_virus_read_count.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run --count-method reads
Rscript Scripts/06_viraquant.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run --abundance-field mean_depth_ge_1 --normalization-mode both
Rscript Scripts/07_integrate_compare.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run
Rscript Scripts/08_compare_kraken2_rvdb_taxid.R --results-dir Results --sample-info Results/sample_info.tsv --output-dir Analyses/full_run --count-method reads
```

### 3. Build the representative sample-data subset only

```bash
Rscript Scripts/90_make_sample_data.R
```

### 4. Run the Kraken2 vs RVDBv31 taxid comparison

All available samples and mapper combinations:

```bash
Rscript Scripts/08_compare_kraken2_rvdb_taxid.R \
  --results-dir Results \
  --sample-info Results/sample_info.tsv \
  --output-dir Analyses/kraken2_rvdb_taxid_comparison \
  --count-method reads
```

One sample and one host mapper only:

```bash
Rscript Scripts/08_compare_kraken2_rvdb_taxid.R \
  --results-dir Results \
  --sample-info Results/sample_info.tsv \
  --output-dir Analyses/kraken2_rvdb_taxid_comparison_1S1 \
  --count-method reads \
  --sample-id 1_S1 \
  --host-mapper bwa_mem
```

One sample, one host mapper, and one virus mapper only:

```bash
Rscript Scripts/08_compare_kraken2_rvdb_taxid.R \
  --results-dir Results \
  --sample-info Results/sample_info.tsv \
  --output-dir Analyses/kraken2_rvdb_taxid_comparison_1S1_bowtie2 \
  --count-method reads \
  --sample-id 1_S1 \
  --host-mapper bwa_mem \
  --virus-mapper bowtie2
```

### 5. Common output locations

- Main pipeline: `Analyses/<run_name>/tables`, `Analyses/<run_name>/plots`, `Analyses/<run_name>/logs`
- Kraken2 taxid comparison: `Analyses/<run_name>/tables`, `Analyses/<run_name>/logs`
- Sample fixture: `sample_data/Results`
