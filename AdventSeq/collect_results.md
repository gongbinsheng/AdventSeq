# collect_results.py

## Overview

Collects the per-sample outputs of an AdventSeq pipeline run into consolidated
tables. The pipeline writes each sample's files into a per-sample working directory
(`Pipeline.WD/$SID`); this tool walks that output tree, matches the files belonging
to every sample listed in the SRA metadata sheet, and assembles them into tables.

Installed as the console command `collect-results` (see the main
[README.md](../README.md) for setup). Run it with the `AdventSeq` environment active.

```bash
conda activate AdventSeq
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
- `--top-n` — comma-separated top-N value(s) for the scan top-N interactive table,
  e.g. `10,20,50` (default `10`). One tab is produced per value.

## File discovery

Discovery is robust to both output layouts:

- per-sample folders — `<results_dir>/<library_ID>/<library_ID>.*`
- flattened per-type — `<results_dir>/<type>/<library_ID>.*`

Files are matched by basename prefix `<library_ID>.` (the trailing dot prevents
`1_S1` from matching `10_S10`). The host/virus aligner and the viral reference are
parsed from the embedded `ref_genome` tokens in the file names.

## Outputs

Each table is written as a `.tsv` and bundled as a sheet in `AdventSeq_results.xlsx`:

| Table | Source | Shape |
|---|---|---|
| `sample_info` | `*.fastp.json` | sample ID, title, total reads before/after filtering, passed-filter reads, after/before fraction |
| `host_featureCounts_by_reads` | `*.featureCounts.*.count_reads.tsv.gz` | gene × sample counts |
| `host_featureCounts_by_pairs` | `*.featureCounts.*.count_pairs.tsv.gz` | gene × sample counts |
| `kraken2` | `*.Kraken2_viral.report` | taxon × sample (clade-assigned reads) |
| `scan_taxonomy_reads` / `scan_taxonomy_pairs` | `*<scan-ref>*.read_count.txt` | virus × sample direct counts |
| `scan_summary` | same | per-sample mapped / primary / discordant / virus totals |
| `viraquant_scan_coverage` | `*<scan-ref>*.ViraQuant_scan_ncbitaxon.tsv` (preferred), else `*<scan-ref>*.ViraQuant_scan.tsv` | long/tidy coverage (variable contigs / ncbitaxon groups); prefers the ncbitaxon-grouped scan output and falls back to the raw scan TSV per sample; written **gzip-compressed** (`.tsv.gz`) and **not** included in the Excel workbook (it is large) |
| `targeted_taxonomy_reads_<ref>` / `_pairs_<ref>` | `*<ref>*.read_count.txt` | virus × sample direct counts |
| `targeted_summary_<ref>` | same | per-sample totals |
| `viraquant_targeted_coverage_<ref>` | `*<ref>*.ViraQuant.tsv` | long/tidy coverage |

### CPM-normalized workbook

A companion `AdventSeq_results_CPM.xlsx` is written alongside `AdventSeq_results.xlsx`
with the **same sheets** normalized to **counts per million filtered reads** (CPM):
each value is divided by its sample's fastp after-filtering total reads, per million
(`value / (filtered_reads / 1e6)`). A sample whose filtered-read count is missing or
zero is left blank in the normalized sheet.

- **Count matrices** (`host_featureCounts_*`, `kraken2`, `scan_taxonomy_*`,
  `targeted_taxonomy_*`) are normalized **column-wise** — each per-sample column is
  divided by that sample's reads-per-million.
- **`viraquant_targeted_coverage_<ref>`** is a long table (one row per sample), so it
  is normalized **row-wise**, and **only on its measurement columns**. The id/info
  columns — `level`, `virus`, `ncbitaxonname`, `seq_id`, `n_contigs`,
  `best_contig_id`, `length` — and the boolean `pass_*` flags are **left unchanged**.

The `sample_info`, `scan_summary`, and `targeted_summary_<ref>` sheets have no
per-sample count columns and are **copied through unchanged**.
`viraquant_scan_coverage` is excluded from both workbooks (it is large; written as
`.tsv.gz`).

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
- **Threshold + Apply** — each measurement has **two** thresholds: a **raw** one and a
  **normalized** one. The threshold applies to the value as displayed, and the input
  switches with the Normalized toggle. Any cell whose displayed value is below the
  active threshold is shown empty. Both the raw and normalized thresholds default to
  `0` for every measurement (no filtering — tune them). Defaults come from the script's
  `SUGGESTED_THRESHOLDS` / `SUGGESTED_THRESHOLDS_NORM` dicts (add entries and regenerate
  to pre-set non-zero built-in defaults). Stepping the value
  with the up/down spinner arrows (or the mouse wheel / Up–Down keys) updates the table
  **immediately**; a number typed by hand still applies on **Apply** (or Enter). The
  spinner step adapts to the measurement type and view mode: raw uses `1` (integer) /
  `0.1` (float), normalized uses the finer `0.1` (integer) / `0.05` (float) since
  per-million values are much smaller. Edits persist per-browser via `localStorage`
  (works even when the page is opened directly as a `file://` document).
- **Download table (.tsv)** — downloads two TSVs for the current measurement and
  reference (to your browser's downloads folder): `viraquant_<ref>_<measurement>_raw.tsv`
  (raw values, all cells) and `…_filtered.tsv` (the displayed values — raw or normalized
  — with the active threshold applied; cells below it left empty, i.e. blank/NA in Excel).
- **Highlight a column** — click a virus/entity column header to highlight that whole
  column; click it again to clear. Highlights follow the entity across measurement and
  top-N tab changes.

### Scan top-N interactive table

`viraquant_scan_topn_interactive.html` is the same self-contained viewer applied to
the scan coverage (raw or ncbitaxon-grouped). Because a scan covers far more entities
than the targeted lists, it shows only the **top N** per measurement: for the selected
measurement each sample's top-N entities are pooled into a **union** (so a tab may show
more than N columns), and the columns are ordered by how many samples' top-N include
each entity (most-shared first). The **`--top-n` tabs** (e.g. Top 10 / Top 20 / Top 50)
switch N; the measurement drop-down, Normalized toggle, thresholds, and Download work as
above. Because normalization is a per-sample constant scaling, the top-N membership is
the same in raw and normalized views.

Column headers are rendered **vertically** (labels are long) and show the
`ncbitaxonname` for grouped data or the `seq_id` (contig id) for raw scan data; the
NCBI taxon id (grouped) is available as a tooltip. `pass_*` booleans, `n_contigs`, and
`length` are excluded from the measurement list (they cannot be ranked).

## Related files

- `AdventSeq_Pipeline.py` — builds the per-sample pipeline whose outputs are collected.
- `ViraQuant.py`, `TaxonomyClassifier.py` — produce the coverage and read-count inputs.
