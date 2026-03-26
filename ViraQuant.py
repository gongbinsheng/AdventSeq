#!/usr/bin/env python3
"""
ViraQuant.py

Quantify viruses from a sorted+indexed BAM by per-base depth summaries.

Outputs a single TSV containing:
  - contig-level rows (level=contig)
  - virus-level rows (level=virus) when a 2-column virus list is provided (col2 groups contigs)

Virus list file formats:
  1 column:  contig_id
  2 columns: contig_id <whitespace> virus_name   (contig_id must be unique; virus_name may repeat)

If --viruses is omitted:
  report top N contigs by mapped_reads/length (default N=10).

Two metric families:
A) Depth quantiles (include zeros):
   depth_at_P: depth d such that at least P% of positions have depth >= d.

B) Threshold/breadth method:
   For chosen n (default 1):
     frac_ge_n = fraction of positions with depth >= n
     mean_depth_ge_n, median_depth_ge_n computed over positions with depth >= n only.
   pass_k_* flags indicate whether frac_ge_n >= k% for each k in --kpercents.

Virus-level aggregation (when 2-col virus list is provided):
  - pooled per-base depth histogram across member contigs (position-weighted)
  - pooled/weighted depth quantiles, breadth, frac_ge_n, mean/median over depth>=n positions
  - plus "best contig" detection-friendly fields in TSV
"""

import sys
import pysam
import argparse

from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_PCTS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def parse_percents(s: str):
    parts = [p.strip() for p in s.split(",") if p.strip()]
    vals = [float(p) for p in parts]
    # If user supplies 50,60,... treat as percent
    if any(v > 1.0 for v in vals):
        vals = [v / 100.0 for v in vals]
    for v in vals:
        if v < 0 or v > 1:
            raise ValueError(f"Percent/fraction out of range [0,1]: {v}")
    return sorted(vals)


def load_virus_map_and_grouping(virus_arg: str):
    """
    Returns:
      ids: list of contig IDs (order preserved)
      id_to_name: contigID -> output label (virus_name if 2-col, else contigID)
      has_two_cols: bool
      name_to_ids: virus_name -> list of contig IDs (only if has_two_cols)
    """
    if virus_arg is None:
        return None, None, False, None

    # Try file first
    try:
        ids = []
        id_to_name = {}
        name_to_ids = defaultdict(list)
        seen = set()
        has_two_cols = False

        with open(virus_arg, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if not parts:
                    continue
                vid = parts[0]
                if vid in seen:
                    raise ValueError(
                        f"Duplicate contig ID in virus list (col1 must be unique): {vid}"
                    )
                seen.add(vid)

                if len(parts) >= 2:
                    has_two_cols = True
                    vname = parts[1]
                else:
                    vname = vid

                ids.append(vid)
                id_to_name[vid] = vname
                if has_two_cols:
                    name_to_ids[vname].append(vid)

        if ids:
            if not has_two_cols:
                name_to_ids = None
            return ids, id_to_name, has_two_cols, name_to_ids
    except FileNotFoundError:
        pass

    # Otherwise treat as comma-separated list
    ids = [x.strip() for x in virus_arg.split(",") if x.strip()]
    id_to_name = {vid: vid for vid in ids}
    return ids, id_to_name, False, None


def depth_at_fraction_positions(depth_counts: Counter, genome_len: int, frac: float) -> int:
    """
    Returns d such that at least frac of positions have depth >= d.
    Equivalent to (1-frac) quantile INCLUDING zeros.
    """
    if genome_len <= 0:
        return 0

    below_target = int((1.0 - frac) * genome_len)

    depths_sorted = sorted(depth_counts.items(), key=lambda kv: kv[0])
    cum_below = 0
    for depth, npos in depths_sorted:
        cum_below += npos
        if cum_below > below_target:
            return depth

    return depths_sorted[-1][0] if depths_sorted else 0


def get_idxstats_table(bam_path: str):
    """
    Returns dict contig -> (length, mapped_reads, unmapped_reads)
    (idxstats also returns '*' for unmapped; caller should ignore it)
    """
    stats_txt = pysam.idxstats(bam_path)
    table = {}
    for line in stats_txt.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        contig = parts[0]
        length = int(parts[1])
        mapped = int(parts[2])
        unmapped = int(parts[3])
        table[contig] = (length, mapped, unmapped)
    return table


def choose_top_contigs_by_norm_mapped(lengths, mapped_reads, top_n: int):
    scored = []
    for c, L in lengths.items():
        if c == "*" or L <= 0:
            continue
        m = mapped_reads.get(c, 0)
        score = m / L
        scored.append((score, m, L, c))
    scored.sort(reverse=True, key=lambda x: (x[0], x[1], x[2]))
    return [c for _, _, _, c in scored[: max(0, top_n)]]


def compute_depth_hist_and_stats(
    bam: pysam.AlignmentFile,
    contig: str,
    length: int,
    min_mapq: int,
    flag_filter: int,
    max_depth: int,
):
    """
    Build depth histogram for all positions (including zeros) plus mean depth and breadth>0.
    """
    depth_counts = Counter()
    covered_positions = 0

    for col in bam.pileup(
        contig=contig,
        start=0,
        end=length,
        truncate=True,
        stepper="all",
        min_mapping_quality=min_mapq,
        flag_filter=flag_filter,
        max_depth=max_depth,
    ):
        d = col.nsegments
        depth_counts[d] += 1
        covered_positions += 1

    zeros = max(0, length - covered_positions)
    if zeros:
        depth_counts[0] += zeros

    total_depth = sum(depth * npos for depth, npos in depth_counts.items())
    mean_depth = (total_depth / length) if length > 0 else 0.0
    breadth = 1.0 - (depth_counts.get(0, 0) / length if length > 0 else 1.0)

    return depth_counts, mean_depth, breadth


def mean_median_depth_ge_n(depth_counts: Counter, n_min: int):
    """
    Compute frac_ge_n, mean_ge_n, median_ge_n over positions with depth >= n_min.

    Median computed exactly from histogram (no per-base array).
    Returns:
      frac_ge, mean_ge (float or None), median_ge (float or None)
    """
    total_positions = sum(depth_counts.values())
    if total_positions == 0:
        return 0.0, None, None

    ge_counts = [(d, c) for d, c in depth_counts.items() if d >= n_min]
    ge_npos = sum(c for _, c in ge_counts)
    frac_ge = ge_npos / total_positions if total_positions else 0.0
    if ge_npos == 0:
        return frac_ge, None, None

    total_depth_ge = sum(d * c for d, c in ge_counts)
    mean_ge = total_depth_ge / ge_npos

    # median over ge positions
    mid1 = (ge_npos - 1) // 2
    mid2 = ge_npos // 2
    ge_counts_sorted = sorted(ge_counts, key=lambda x: x[0])

    cum = 0
    med1 = med2 = None
    for d, c in ge_counts_sorted:
        cum += c
        if med1 is None and mid1 < cum:
            med1 = d
        if med2 is None and mid2 < cum:
            med2 = d
        if med1 is not None and med2 is not None:
            break

    median_ge = (med1 + med2) / 2.0 if med1 is not None else None
    return frac_ge, mean_ge, median_ge


def fmt(x):
    return "NA" if x is None else x


def main():
    ap = argparse.ArgumentParser(
        description="Quantify virus contigs (and virus groups) from BAM using depth quantiles and breadth/threshold metrics."
    )
    ap.add_argument("--bam", required=True, help="Sorted BAM (must be indexed with .bai).")
    ap.add_argument(
        "--viruses",
        default=None,
        help=(
            "Either: (a) virus_list.txt with 1-2 columns: contig_id [virus_name], "
            "or (b) comma-separated contig IDs. If omitted, use top N contigs by mapped_reads/length."
        ),
    )
    ap.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="When --viruses is omitted, report top N contigs by mapped_reads/length. Default 10.",
    )
    ap.add_argument(
        "--percents",
        default=None,
        help="Quantile percents list for depth_at_P, e.g. '50,60,70,80,90,100' or '0.5,0.6,...'.",
    )
    ap.add_argument(
        "--kpercents",
        default=None,
        help="Breadth thresholds k% for pass_k_* flags. Default: same as --percents.",
    )
    ap.add_argument("--n-min", type=int, default=1, help="n in 'depth >= n reads'. Default 1.")
    ap.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ for pileup. Default 0.")
    ap.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include PCR/optical duplicates. Default filters duplicates.",
    )
    ap.add_argument("--max-depth", type=int, default=1000000, help="Pileup max depth. Default 1,000,000.")
    ap.add_argument("--out", default="-", help="Output TSV path (default stdout).")

    args = ap.parse_args()

    if args.n_min < 1:
        print("ERROR: --n-min must be >= 1", file=sys.stderr)
        sys.exit(2)

    percents = parse_percents(args.percents) if args.percents else DEFAULT_PCTS
    kpercents = parse_percents(args.kpercents) if args.kpercents else list(percents)

    bam = pysam.AlignmentFile(args.bam, "rb")
    lengths = dict(zip(bam.references, bam.lengths))

    idx_table = get_idxstats_table(args.bam)
    mapped_reads = {c: v[1] for c, v in idx_table.items() if c != "*"}

    # Filter flags: unmapped, secondary, QCfail, supplementary; and duplicates unless included
    flag_filter = 0x4 | 0x100 | 0x200 | 0x800
    if not args.include_duplicates:
        flag_filter |= 0x400

    virus_ids, id_to_name, has_two_cols, name_to_ids = load_virus_map_and_grouping(args.viruses)

    if virus_ids is None:
        if args.top_n <= 0:
            print("ERROR: --top-n must be >= 1 when --viruses is omitted.", file=sys.stderr)
            sys.exit(2)
        virus_ids = choose_top_contigs_by_norm_mapped(lengths, mapped_reads, args.top_n)
        id_to_name = {vid: vid for vid in virus_ids}
        has_two_cols = False
        name_to_ids = None

    # Validate contigs exist
    missing = [v for v in virus_ids if v not in lengths]
    if missing:
        print("WARNING: contigs not in BAM header will be skipped:", file=sys.stderr)
        for v in missing:
            print(f"  {v}", file=sys.stderr)
        virus_ids = [v for v in virus_ids if v in lengths]

    if args.out == "-":
        out_fh = sys.stdout
    else:
        out_path = Path(args.out).expanduser().resolve()
        print(f"Writing output to\n  {out_path}", file=sys.stderr)
        out_fh = open(out_path, "w", encoding="utf-8")

    # TSV header
    header = (
        [
            "level",  # contig or virus
            "virus",  # contig row: label (virus_name or contig_id); virus row: virus_name
            "seq_id",  # contig row: contig_id; virus row: comma-separated contig_ids
            "n_contigs",  # contig row: 1; virus row: number of contigs grouped
            "length",
            "mapped_reads",
            "mapped_per_bp",
            "breadth_cov_gt0",
            "mean_depth_all",
            f"frac_ge_{args.n_min}",
            f"mean_depth_ge_{args.n_min}",
            f"median_depth_ge_{args.n_min}",
        ]
        + [f"pass_k_{int(k*100)}pct_ge_{args.n_min}" for k in kpercents]
        + [f"depth_at_{int(p*100)}pct" for p in percents]
        + [
            # Added to TSV per your request
            "best_contig_id",
            "best_contig_mapped_per_bp",
            "best_contig_depth_at_90pct",
            f"best_contig_frac_ge_{args.n_min}",
        ]
    )
    print("\t".join(header), file=out_fh)

    # Collect per-contig depth histograms and stats to support virus-level pooling
    contig_depth_counts = {}  # vid -> Counter(depth -> npos)
    contig_stats = {}  # vid -> dict

    for vid in virus_ids:
        L = lengths[vid]
        m = mapped_reads.get(vid, 0)
        mapped_per_bp = (m / L) if L > 0 else 0.0

        dcounts, mean_all, breadth_gt0 = compute_depth_hist_and_stats(
            bam=bam,
            contig=vid,
            length=L,
            min_mapq=args.min_mapq,
            flag_filter=flag_filter,
            max_depth=args.max_depth,
        )

        frac_ge, mean_ge, median_ge = mean_median_depth_ge_n(dcounts, args.n_min)

        contig_depth_counts[vid] = dcounts
        contig_stats[vid] = {
            "L": L,
            "m": m,
            "mapped_per_bp": mapped_per_bp,
            "breadth_gt0": breadth_gt0,
            "mean_all": mean_all,
            "frac_ge": frac_ge,
            "mean_ge": mean_ge,
            "median_ge": median_ge,
        }

        label = id_to_name.get(vid, vid)

        # contig row: best_contig_* = self (makes TSV rectangular)
        best_depth90 = depth_at_fraction_positions(dcounts, L, 0.9)

        row = [
            "contig",
            label,
            vid,
            "1",
            str(L),
            str(m),
            f"{mapped_per_bp:.6g}",
            f"{breadth_gt0:.6f}",
            f"{mean_all:.6f}",
            f"{frac_ge:.6f}",
            (f"{mean_ge:.6f}" if mean_ge is not None else "NA"),
            (f"{median_ge:.6f}" if median_ge is not None else "NA"),
        ]
        for k in kpercents:
            row.append("1" if frac_ge >= k else "0")
        for p in percents:
            row.append(str(depth_at_fraction_positions(dcounts, L, p)))

        row += [
            vid,
            f"{mapped_per_bp:.6g}",
            str(best_depth90),
            f"{frac_ge:.6f}",
        ]

        print("\t".join(row), file=out_fh)

    # Virus-level rows (only if 2-column list provided)
    if has_two_cols and name_to_ids:
        for vname, vids in name_to_ids.items():
            vids = [v for v in vids if v in contig_stats]
            if not vids:
                continue

            virus_len = sum(contig_stats[v]["L"] for v in vids)
            virus_mapped = sum(contig_stats[v]["m"] for v in vids)
            virus_mapped_per_bp = (virus_mapped / virus_len) if virus_len > 0 else 0.0

            # Pool histograms across contigs (position-weighted)
            pooled = Counter()
            for v in vids:
                pooled.update(contig_depth_counts[v])

            total_depth_all = sum(d * c for d, c in pooled.items())
            virus_mean_all = (total_depth_all / virus_len) if virus_len > 0 else 0.0
            virus_breadth_gt0 = 1.0 - (pooled.get(0, 0) / virus_len if virus_len > 0 else 1.0)

            virus_frac_ge, virus_mean_ge, virus_median_ge = mean_median_depth_ge_n(pooled, args.n_min)

            # Best contig within this virus group (detection-friendly)
            best_vid = max(vids, key=lambda x: contig_stats[x]["mapped_per_bp"])
            best_mapped_per_bp = contig_stats[best_vid]["mapped_per_bp"]
            best_depth90 = depth_at_fraction_positions(
                contig_depth_counts[best_vid], contig_stats[best_vid]["L"], 0.9
            )
            best_frac_ge = contig_stats[best_vid]["frac_ge"]

            row = [
                "virus",
                vname,
                ",".join(vids),
                str(len(vids)),
                str(virus_len),
                str(virus_mapped),
                f"{virus_mapped_per_bp:.6g}",
                f"{virus_breadth_gt0:.6f}",
                f"{virus_mean_all:.6f}",
                f"{virus_frac_ge:.6f}",
                (f"{virus_mean_ge:.6f}" if virus_mean_ge is not None else "NA"),
                (f"{virus_median_ge:.6f}" if virus_median_ge is not None else "NA"),
            ]
            for k in kpercents:
                row.append("1" if virus_frac_ge >= k else "0")
            for p in percents:
                row.append(str(depth_at_fraction_positions(pooled, virus_len, p)))

            row += [
                best_vid,
                f"{best_mapped_per_bp:.6g}",
                str(best_depth90),
                f"{best_frac_ge:.6f}",
            ]

            print("\t".join(row), file=out_fh)

    if out_fh is not sys.stdout:
        out_fh.close()
    bam.close()


if __name__ == "__main__":
    main()
