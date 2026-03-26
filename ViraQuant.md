# ViraQuant.py

Quantify viruses from Illumina paired-end read mappings (BWA → sorted/indexed BAM) using **per-base depth** summaries.

The script reports both:

1) **Contig-level** depth summaries (each reference sequence / RVDB entry / FASTA record).
2) **Virus-level** summaries when your `virus_list.txt` provides a virus name that groups multiple contigs (e.g., strains, incomplete genomes, segments).

Output is a **single TSV** with a `level` column (`contig` or `virus`).

---

## Requirements

- Python 3.12+
- `pysam`
- Input BAM must be:
  - coordinate-sorted
  - indexed (`.bai` present)

Install `pysam` (example):

```bash
pip install pysam
```

## Inputs
### 1) Sorted, indexed BAM

Example:
- sample.sorted.bam
- sample.sorted.bam.bai

### 2) Optional virus list (--viruses)

You can provide either:
- A path to a text file (virus_list.txt), 1 or 2 columns (whitespace-delimited):

    1 column
    ```
    contig_id
    ```
    
    2 columns
    ```
    contig_id   virus_name
    ```

> Rules:
> - Column 1 must be unique (contig IDs / reference names in BAM header).
> - Column 2 may repeat (multiple contigs can belong to the same virus name).

Or:

- A comma-separated list of contig IDs:
```bash
--viruses contigA,contigB,contigC
```

If `--viruses` is not provided, the script automatically reports the top N contigs ranked by:
`mapped_reads / contig_length` (default N=10; adjust with --top-n).

## Output (TSV)
One row per contig, plus optional aggregated virus rows.
Key output columns:

### Identification
- level: contig or virus
- virus:
  - contig rows: label (virus_name if 2-col list, otherwise contig_id)
  - virus rows: virus_name (group)
- seq_id:
  - contig rows: contig ID
  - virus rows: comma-separated list of contig IDs in the group
- n_contigs:
  - contig rows: 1
  - virus rows: number of contigs in the group

###Basic mapping / depth
- `length`: contig length (virus rows: sum of lengths across contigs)
- `mapped_reads`: from samtools idxstats/pysam.idxstats (virus rows: sum across contigs)
- `mapped_per_bp`: mapped_reads / length
- `breadth_cov_gt0`: fraction of positions with depth > 0
- `mean_depth_all`: mean depth across all positions, including zeros

### Method A: Depth quantiles (include zeros)
- For each P in --percents (default 50..100 step 10):
  - `depth_at_80pct` = depth d such that 80% of positions have depth ≥ d
> This matches statements like:
> - “5 counts for 80% of genome positions”
> - “3 counts for 90% of genome positions”

### Method B: Threshold/breadth + mean/median among covered bases
For chosen --n-min (default 1):
- `frac_ge_1`: fraction of positions with depth ≥ 1
- `mean_depth_ge_1`: mean depth among positions with depth ≥ 1
- `median_depth_ge_1`: median depth among positions with depth ≥ 1
- `pass_k_80pct_ge_1`: 1 if frac_ge_1 >= 0.80, else 0

### Virus-level aggregation details (when 2-col virus list is used)
Virus rows are computed by pooling per-base depth histograms across member contigs (position-weighted), i.e. equivalent to concatenating contigs end-to-end.

This is recommended because it:
- naturally handles different contig lengths
- avoids “averaging quantiles” (which is statistically awkward)
- yields interpretable statements like:
  - “Across all contigs labeled X, 90% of positions have depth ≥ d”

### Best-contig fields (included in TSV)
To support detection vs quantification:
- `best_contig_id`: contig with highest mapped_per_bp within the virus group
- `best_contig_mapped_per_bp`
- `best_contig_depth_at_90pct`
- `best_contig_frac_ge_1` (or chosen n)

- For contig rows, these fields are self-referential (best contig = the contig itself) to keep the TSV rectangular.

## Usage Examples

### 1) Analyze a virus list (with virus-level grouping)
```bash
python ViraQuant.py \
  --bam sample.sorted.bam \
  --viruses virus_list.txt \
  --out virus_metrics.tsv
```

### 2) Analyze a single contig
```bash
python ViraQuant.py \
  --bam sample.sorted.bam \
  --viruses NC_045512.2 \
  --out sars2_metrics.tsv
```

### 3) No virus list: report top 10 contigs by mapped_reads/length
```bash
python ViraQuant.py \
  --bam sample.sorted.bam \
  --viruses virus_list.txt \
  --percents 80,90,95,100 \
  --kpercents 70,80,90 \
  --n-min 3 \
  --min-mapq 20 \
  --out virus_metrics.tsv
```

### 4) Customize quantiles and breadth thresholds; require depth ≥ 3
```bash
python ViraQuant.py \
  --bam sample.sorted.bam \
  --viruses virus_list.txt \
  --percents 80,90,95,100 \
  --kpercents 70,80,90 \
  --n-min 3 \
  --min-mapq 20 \
  --out virus_metrics.tsv
```

## Notes / Caveats

- `mapped_reads` comes from idxstats and can be influenced by multi-mapping across highly similar references (common in RVDB).
  - Depth-based metrics are typically more robust than raw mapped read totals.
- Virus-level pooling assumes contigs sharing the same virus name are reasonably comparable (segments, strains, partials).
  - If contigs represent many different strains, pooled virus-level metrics describe “support for the virus name” rather than a single strain.

## Recommended primary metrics

For quantification-like interpretation:
- Virus rows:
  - `depth_at_80pct` / `depth_at_90pct` (pooled)
  - `frac_ge_n` + `mean_depth_ge_n`
- For detection:
  - `best_contig_depth_at_90pct`
  - `best_contig_frac_ge_n`
---
If you want the virus-level “best contig” to be chosen by something *other than* `mapped_per_bp` (e.g., highest `depth_at_90pct` or highest `frac_ge_n`), tell me your preference and I’ll swap the ranking key.







