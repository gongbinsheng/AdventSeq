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
  as a substring of the file names (required), e.g. `RVDBv29` (the RVDB release used
  at run time determines this, e.g. `RVDBv31`).
- `--targeted-refs` — comma-separated list of targeted virus-list reference names
  (required), e.g. `7viruses,9viruses`. One set of targeted tables and one selectable
  reference in the interactive HTML is produced per name.
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
external dependencies). It shows a **sample (rows) × virus (columns)** table — each row
is a sample, with its `library_ID` and `title` in two sticky left columns — whose cells
update when a measurement is chosen from the drop-down; a second drop-down selects the
targeted reference when more than one is given. Virus-level ViraQuant values are shown
(falling back to contig-level for any sample without virus rows). `length` is a contig
length and is not offered as a measurement.

Controls:

- **Normalized** toggle (default on) — divides each value by the sample's filtered
  (fastp after-filtering) total reads per million (CPM-style). `pass_*` booleans are
  shown as green/empty blocks and are never normalized.
- **Threshold + Apply** — a per-measurement value; any cell whose **raw** value is
  below it is shown empty. Defaults are `1` for integer measurements and `0.1` for
  floats, seeded from `viraquant_targeted_thresholds.yml` written beside the HTML. The
  threshold always compares the raw value (the toggle only changes what number is
  displayed). Click **Apply** (or press Enter in the box) to update the table with the
  new threshold and save all thresholds back to `viraquant_targeted_thresholds.yml`
  (a browser download). Edits also persist in the browser via `localStorage`; on open
  the loading order is `localStorage` > fetched sidecar YAML > embedded defaults (the
  sidecar `fetch` works when the page is served over http, e.g. `python -m http.server`).
- **Download table (.tsv)** — saves two TSVs of raw values for the current measurement
  and reference: `viraquant_<ref>_<measurement>_raw.tsv` (all cells) and
  `…_filtered.tsv` (cells below the threshold left empty, i.e. blank/NA in Excel).

## Related files

- `AdVentSeq_Pipeline.py` — builds the per-sample pipeline whose outputs are collected.
- `ViraQuant.py`, `TaxonomyClassifier.py` — produce the coverage and read-count inputs.
