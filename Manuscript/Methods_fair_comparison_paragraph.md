# Fair cross-method comparison of untargeted viral screening (Methods text)

To compare Kraken2, RVDB read-count screening, and ViraQuant scan mode on equal terms, all
detections were mapped to species-level NCBI Taxonomy identifiers (ranks S, S1, S2) and
restricted to a shared background of taxa represented in both the Kraken2 viral database and
RVDB, and agreement between the two database-search methods was summarized as taxids detected
by both, by Kraken2 only, and by RVDB only. To make thresholds comparable across libraries,
read counts were normalized to reads per million using a single per-sample factor
(1,000,000 / total reads passing fastp filtering) applied to every taxon for both Kraken2 and
RVDB, and a taxon was scored as detected at RPM >= 1. ViraQuant scan mode was held to a
stricter, coverage-based criterion, requiring at least half of the reference covered at >= 1x
(depth_at_50pct >= 1); it therefore reports fewer but breadth-confirmed detections, whereas the
read-count methods retain many taxa supported only by low-breadth, locally piled-up reads.

---

*Notes (not for the manuscript): per-sample factor `1e6 / after_filtering_total_reads`; knobs in `Scripts/08_compare_kraken2_rvdb_taxid_fair.R` are `--min-rpm` (default 1), `--background-ranks` (S,S1,S2), and the `depth_at_50pct >= 1` breadth cutoff. The low ViraQuant_scan counts in Figure 2B are expected: across 1_S1, only 1-3 of several hundred-to-thousand scan taxa per sample reach 50% breadth, while ~40% of scan taxa are present in the background (so the collapse is driven by the breadth criterion, not the background restriction).*
