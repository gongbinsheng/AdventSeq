# collect_results.py

## Overview

Collects the per-sample outputs of an AdVentSeq pipeline run into consolidated
tables. The pipeline writes each sample's files into a per-sample working directory
(`Pipeline.WD/$SID`); this tool walks that output tree, matches the files belonging
to every sample listed in the SRA metadata sheet, and assembles them into tables.

Installed as the console command `collect-results` (see the main
[README.md](../README.md) for setup). Run it with the `AdVentSeq` environment active.

```bash
conda activate AdVentSeq
collect-results --help
```

## Usage

```bash
collect-results \
  --sra SRA_metadata_example.xlsx \
  --results-dir ./Results \
  --out-dir ./Collected_Tables \
  --scan-ref RVDBv29 \
  --targeted-refs 7viruses
```

## Arguments

- `-s, --sra` — SRA metadata workbook (`.xlsx`). Required. Uses `library_ID`
  (sample ID) and `title` (sample title).
- `--sheet` — worksheet name (default `SRA_data`).
- `-r, --results-dir` — pipeline output directory (the `-o` directory of the run).
  Required.
- `-o, --out-dir` — output directory for the tables (default `./collected_tables`).
- `--scan-ref` — reference name of the full-BAM scan / direct-count step, matched
  as a substring of the file names (default `RVDBv29`). The RVDB release used at run
  time determines this (e.g. `RVDBv31`).
- `--targeted-refs` — comma-separated list of targeted virus-list reference names,
  e.g. `7viruses,9viruses` (default `7viruses`). One set of targeted tables and one
  selectable reference in the interactive HTML is produced per name.
- `--host-mapper` / `--virus-mapper` — aligner names. Auto-detected when exactly one
  is present in the run; required only when the output mixes multiple mappers.
- `--taxonomy-map` — optional `ncbi_taxonomy_map.json[.gz]`; adds a readable `name`
  column to the scan-ref taxonomy tables (which are keyed by NCBI taxon id).

## File discovery

Discovery is robust to both output layouts:

- per-sample folders — `<results_dir>/<library_ID>/<library_ID>.*`
- flattened per-type — `<results_dir>/<type>/<library_ID>.*`

Files are matched by basename prefix `<library_ID>.` (the trailing dot prevents
`1_S1` from matching `10_S10`). The host/virus aligner and the viral reference are
parsed from the embedded `ref_genome` tokens in the file names.

## Outputs

Each table is written as a `.tsv` and bundled as a sheet in `AdVentSeq_results.xlsx`:

| Table | Source | Shape |
|---|---|---|
| `sample_info` | `*.fastp.json` | sample ID, title, total reads before/after filtering, passed-filter reads, after/before fraction |
| `host_featureCounts_by_reads` | `*.featureCounts.*.count_reads.tsv.gz` | gene × sample counts |
| `host_featureCounts_by_pairs` | `*.featureCounts.*.count_pairs.tsv.gz` | gene × sample counts |
| `kraken2` | `*.Kraken2_viral.report` | taxon × sample (clade-assigned reads) |
| `scan_taxonomy_reads` / `scan_taxonomy_pairs` | `*<scan-ref>*.read_count.txt` | virus × sample direct counts |
| `scan_summary` | same | per-sample mapped / primary / discordant / virus totals |
| `viraquant_scan_coverage` | `*<scan-ref>*.ViraQuant_scan.tsv` | long/tidy coverage (variable contigs) |
| `targeted_taxonomy_reads_<ref>` / `_pairs_<ref>` | `*<ref>*.read_count.txt` | virus × sample direct counts |
| `targeted_summary_<ref>` | same | per-sample totals |
| `viraquant_targeted_coverage_<ref>` | `*<ref>*.ViraQuant.tsv` | long/tidy coverage |

### Interactive table

`viraquant_targeted_interactive.html` is a self-contained page (data embedded, no
external dependencies). It shows a **virus (rows) × sample (columns)** table whose
cells update when a measurement is chosen from the drop-down; a second drop-down
selects the targeted reference when more than one is given. Virus-level ViraQuant
rows are shown (falling back to contig-level rows for any sample without virus rows).

## Related files

- `AdVentSeq_Pipeline.py` — builds the per-sample pipeline whose outputs are collected.
- `ViraQuant.py`, `TaxonomyClassifier.py` — produce the coverage and read-count inputs.
