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
import json
import operator
import re

from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_PCTS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
SCAN_EXPR_RE = re.compile(
    r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(>=|<=|==|>|<)\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s*$"
)
SCAN_N_MIN_METRIC_RE = re.compile(r"^(frac_ge|mean_depth_ge|median_depth_ge)_(\d+)$")
SCAN_DEPTH_AT_METRIC_RE = re.compile(r"^depth_at_(\d+)pct$")
SCAN_PASS_K_METRIC_RE = re.compile(r"^pass_k_(\d+)pct_ge_(\d+)$")
SCAN_OPERATORS = {
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
}
STATIC_FILTERABLE_METRICS = {
    "length",
    "mapped_reads",
    "mapped_per_bp",
    "breadth_cov_gt0",
    "mean_depth_all",
}


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


def unique_sorted(values):
    return sorted(set(values))


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


def option_was_provided(argv, option: str):
    return any(arg == option or arg.startswith(f"{option}=") for arg in argv)


def build_flag_filter(include_duplicates: bool, include_secondary: bool = False):
    # Always exclude unmapped, QC-fail, and supplementary alignments.
    flag_filter = 0x4 | 0x200 | 0x800
    if not include_secondary:
        flag_filter |= 0x100
    if not include_duplicates:
        flag_filter |= 0x400
    return flag_filter


def get_contig_metric_field_names(n_min: int, percents, kpercents):
    frac_field = f"frac_ge_{n_min}"
    mean_field = f"mean_depth_ge_{n_min}"
    median_field = f"median_depth_ge_{n_min}"
    pass_fields = [f"pass_k_{int(k*100)}pct_ge_{n_min}" for k in kpercents]
    depth_fields = [f"depth_at_{int(p*100)}pct" for p in percents]
    return frac_field, mean_field, median_field, pass_fields, depth_fields


def get_filterable_metric_names(n_min: int, percents, kpercents):
    frac_field, mean_field, median_field, pass_fields, depth_fields = get_contig_metric_field_names(
        n_min, percents, kpercents
    )
    return {
        *STATIC_FILTERABLE_METRICS,
        frac_field,
        mean_field,
        median_field,
        *pass_fields,
        *depth_fields,
    }


def parse_scan_by_expression(expr: str):
    match = SCAN_EXPR_RE.match(expr)
    if not match:
        raise ValueError(
            "Invalid --scan-by expression. Expected exactly one comparison like "
            "'mean_depth_ge_1>=1' or 'depth_at_90pct > 5'."
        )
    metric_name, op_symbol, raw_value = match.groups()
    return metric_name, op_symbol, float(raw_value)


def validate_scan_metric_name(metric_name: str, valid_metric_names):
    if metric_name not in valid_metric_names:
        valid = ", ".join(sorted(valid_metric_names))
        raise ValueError(
            f"Unknown --scan-by metric '{metric_name}'. Valid metrics for this run: {valid}"
        )


def infer_scan_metric_requirements(metric_name: str):
    match = SCAN_N_MIN_METRIC_RE.match(metric_name)
    if match:
        return {"n_min": int(match.group(2)), "percents": [], "kpercents": []}

    match = SCAN_DEPTH_AT_METRIC_RE.match(metric_name)
    if match:
        return {"n_min": None, "percents": [int(match.group(1)) / 100.0], "kpercents": []}

    match = SCAN_PASS_K_METRIC_RE.match(metric_name)
    if match:
        return {
            "n_min": int(match.group(2)),
            "percents": [],
            "kpercents": [int(match.group(1)) / 100.0],
        }

    return {"n_min": None, "percents": [], "kpercents": []}


def format_percent_token(value: float):
    pct = value * 100.0
    if abs(pct - round(pct)) < 1e-9:
        return str(int(round(pct)))
    return f"{pct:g}"


def format_percent_list(values):
    return ",".join(format_percent_token(v) for v in values)


def normalize_percent_list(values):
    return unique_sorted(values)


def make_adjustment_entry(original, effective, reason: str):
    return {
        "original": original,
        "effective": effective,
        "reason": reason,
    }


def print_adjustment_warning(option_name: str, original, effective, reason: str):
    print(
        f"WARNING: {reason}; overriding {option_name} from {original!r} to {effective!r}.",
        file=sys.stderr,
    )


def resolve_effective_runtime_args(args, argv):
    explicit_n_min = option_was_provided(argv, "--n-min")
    explicit_percents = option_was_provided(argv, "--percents")
    explicit_kpercents = option_was_provided(argv, "--kpercents")

    effective_n_min = args.n_min
    effective_percents = normalize_percent_list(parse_percents(args.percents) if args.percents else DEFAULT_PCTS)
    auto_adjustments = {}
    scan_filter = None
    implied = {"n_min": None, "percents": [], "kpercents": []}

    if args.scan_by is not None:
        metric_name, op_symbol, threshold = parse_scan_by_expression(args.scan_by)
        implied = infer_scan_metric_requirements(metric_name)
        scan_filter = (metric_name, op_symbol, threshold)

        if implied["n_min"] is not None and effective_n_min != implied["n_min"]:
            reason = f"--scan-by metric '{metric_name}' implies --n-min {implied['n_min']}"
            if explicit_n_min:
                print_adjustment_warning("--n-min", args.n_min, implied["n_min"], reason)
            auto_adjustments["n_min"] = make_adjustment_entry(args.n_min, implied["n_min"], reason)
            effective_n_min = implied["n_min"]

        merged_percents = normalize_percent_list([*effective_percents, *implied["percents"]])
        if merged_percents != effective_percents:
            reason = (
                f"--scan-by metric '{metric_name}' requires --percents to include "
                f"{format_percent_list(implied['percents'])}"
            )
            if explicit_percents:
                print_adjustment_warning(
                    "--percents",
                    format_percent_list(effective_percents),
                    format_percent_list(merged_percents),
                    reason,
                )
            auto_adjustments["percents"] = make_adjustment_entry(
                list(effective_percents),
                list(merged_percents),
                reason,
            )
            effective_percents = merged_percents

    if args.kpercents:
        effective_kpercents = normalize_percent_list(parse_percents(args.kpercents))
    else:
        effective_kpercents = list(effective_percents)

    if args.scan_by is not None:
        metric_name = scan_filter[0]
        merged_kpercents = normalize_percent_list([*effective_kpercents, *implied["kpercents"]])
        if merged_kpercents != effective_kpercents:
            reason = (
                f"--scan-by metric '{metric_name}' requires --kpercents to include "
                f"{format_percent_list(implied['kpercents'])}"
            )
            if explicit_kpercents:
                print_adjustment_warning(
                    "--kpercents",
                    format_percent_list(effective_kpercents),
                    format_percent_list(merged_kpercents),
                    reason,
                )
            auto_adjustments["kpercents"] = make_adjustment_entry(
                list(effective_kpercents),
                list(merged_kpercents),
                reason,
            )
            effective_kpercents = merged_kpercents

        if not explicit_kpercents:
            baseline_kpercents = normalize_percent_list(parse_percents(args.percents) if args.percents else DEFAULT_PCTS)
            if effective_kpercents != baseline_kpercents and "kpercents" not in auto_adjustments:
                auto_adjustments["kpercents"] = make_adjustment_entry(
                    list(baseline_kpercents),
                    list(effective_kpercents),
                    "--kpercents default was recomputed after --scan-by adjustments",
                )

        validate_scan_metric_name(
            scan_filter[0],
            get_filterable_metric_names(effective_n_min, effective_percents, effective_kpercents),
        )

    return effective_n_min, effective_percents, effective_kpercents, scan_filter, auto_adjustments


def yaml_scalar(value):
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return f"{value:g}" if isinstance(value, float) else str(value)
    return json.dumps(str(value))


def yaml_lines(value, indent: int = 0):
    prefix = " " * indent
    if isinstance(value, dict):
        if not value:
            return [f"{prefix}{{}}"]
        lines = []
        for key, item in value.items():
            if isinstance(item, dict) and not item:
                lines.append(f"{prefix}{key}: {{}}")
            elif isinstance(item, list) and not item:
                lines.append(f"{prefix}{key}: []")
            elif isinstance(item, (dict, list)):
                lines.append(f"{prefix}{key}:")
                lines.extend(yaml_lines(item, indent + 2))
            else:
                lines.append(f"{prefix}{key}: {yaml_scalar(item)}")
        return lines

    if isinstance(value, list):
        if not value:
            return [f"{prefix}[]"]
        lines = []
        for item in value:
            if isinstance(item, (dict, list)):
                lines.append(f"{prefix}-")
                lines.extend(yaml_lines(item, indent + 2))
            else:
                lines.append(f"{prefix}- {yaml_scalar(item)}")
        return lines

    return [f"{prefix}{yaml_scalar(value)}"]


def write_yaml_file(path: Path, payload: dict):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(yaml_lines(payload)))
        fh.write("\n")


def build_run_metadata(args, out_path: Path, n_min: int, percents, kpercents, scan_mode: bool):
    payload = {
        "bam": str(Path(args.bam).expanduser().resolve()),
        "percents": list(percents),
        "kpercents": list(kpercents),
        "n_min": n_min,
        "min_mapq": args.min_mapq,
        "include_duplicates": args.include_duplicates,
        "max_depth": args.max_depth,
        "out": str(out_path),
    }
    if scan_mode:
        payload["scan_by"] = args.scan_by
    elif args.viruses is not None:
        payload["viruses"] = args.viruses
    else:
        payload["top_n"] = args.top_n
    return payload


def passes_scan_filter(metrics, scan_filter):
    if scan_filter is None:
        return True

    metric_name, op_symbol, threshold = scan_filter
    value = metrics.get(metric_name)
    if value is None:
        return False
    return SCAN_OPERATORS[op_symbol](value, threshold)


def compute_contig_result(
    bam: pysam.AlignmentFile,
    contig: str,
    length: int,
    mapped_reads: dict,
    min_mapq: int,
    flag_filter: int,
    max_depth: int,
    n_min: int,
    percents,
    kpercents,
):
    mapped = mapped_reads.get(contig, 0)
    mapped_per_bp = (mapped / length) if length > 0 else 0.0

    depth_counts, mean_all, breadth_gt0 = compute_depth_hist_and_stats(
        bam=bam,
        contig=contig,
        length=length,
        min_mapq=min_mapq,
        flag_filter=flag_filter,
        max_depth=max_depth,
    )

    frac_ge, mean_ge, median_ge = mean_median_depth_ge_n(depth_counts, n_min)
    best_depth90 = depth_at_fraction_positions(depth_counts, length, 0.9)

    frac_field, mean_field, median_field, pass_fields, depth_fields = get_contig_metric_field_names(
        n_min, percents, kpercents
    )
    metrics = {
        "length": length,
        "mapped_reads": mapped,
        "mapped_per_bp": mapped_per_bp,
        "breadth_cov_gt0": breadth_gt0,
        "mean_depth_all": mean_all,
        frac_field: frac_ge,
        mean_field: mean_ge,
        median_field: median_ge,
    }
    for k, field in zip(kpercents, pass_fields):
        metrics[field] = int(frac_ge >= k)
    for p, field in zip(percents, depth_fields):
        metrics[field] = depth_at_fraction_positions(depth_counts, length, p)

    return {
        "depth_counts": depth_counts,
        "length": length,
        "mapped_reads": mapped,
        "mapped_per_bp": mapped_per_bp,
        "breadth_gt0": breadth_gt0,
        "mean_all": mean_all,
        "frac_ge": frac_ge,
        "mean_ge": mean_ge,
        "median_ge": median_ge,
        "best_depth90": best_depth90,
        "metrics": metrics,
    }


def build_header(n_min: int, percents, kpercents):
    frac_field, _, _, pass_fields, depth_fields = get_contig_metric_field_names(n_min, percents, kpercents)
    return (
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
            frac_field,
            f"mean_depth_ge_{n_min}",
            f"median_depth_ge_{n_min}",
        ]
        + pass_fields
        + depth_fields
        + [
            "best_contig_id",
            "best_contig_mapped_per_bp",
            "best_contig_depth_at_90pct",
            f"best_contig_frac_ge_{n_min}",
        ]
    )


def build_contig_row(contig: str, label: str, result, n_min: int, percents, kpercents):
    frac_field, mean_field, median_field, pass_fields, depth_fields = get_contig_metric_field_names(
        n_min, percents, kpercents
    )
    metrics = result["metrics"]
    row = [
        "contig",
        label,
        contig,
        "1",
        str(result["length"]),
        str(result["mapped_reads"]),
        f"{result['mapped_per_bp']:.6g}",
        f"{result['breadth_gt0']:.6f}",
        f"{result['mean_all']:.6f}",
        f"{metrics[frac_field]:.6f}",
        (f"{metrics[mean_field]:.6f}" if metrics[mean_field] is not None else "NA"),
        (f"{metrics[median_field]:.6f}" if metrics[median_field] is not None else "NA"),
    ]
    for field in pass_fields:
        row.append(str(metrics[field]))
    for field in depth_fields:
        row.append(str(metrics[field]))

    row += [
        contig,
        f"{result['mapped_per_bp']:.6g}",
        str(result["best_depth90"]),
        f"{metrics[frac_field]:.6f}",
    ]
    return row


def main():
    argv = sys.argv[1:]
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
        "--scan-by",
        default=None,
        help=(
            "Full-BAM contig scan mode. Compute depth metrics for every contig, ignore --viruses and --top-n, "
            "and report only contigs passing a single numeric comparison such as 'mean_depth_ge_2>=2'. "
            "Dynamic metric names auto-adjust --n-min, --percents, and --kpercents as needed."
        ),
    )
    ap.add_argument(
        "--percents",
        default=None,
        help="Quantile percents list for depth_at_P, e.g. '50,60,70,80,90,100' or '0.5,0.6,...'.",
    )
    ap.add_argument(
        "--kpercents",
        default=None,
        help="Breadth thresholds k%% for pass_k_* flags. Default: same as --percents.",
    )
    ap.add_argument("--n-min", type=int, default=1, help="n in 'depth >= n reads'. Default 1.")
    ap.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ for pileup. Default 0.")
    ap.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include PCR/optical duplicates. Default filters duplicates.",
    )
    ap.add_argument("--max-depth", type=int, default=1000000, help="Pileup max depth. Default 1,000,000.")
    ap.add_argument("--out", required=True, help="Output TSV path. A companion .yml file is always written.")

    args = ap.parse_args(argv)

    if args.n_min < 1:
        print("ERROR: --n-min must be >= 1", file=sys.stderr)
        sys.exit(2)
    if args.out == "-":
        print("ERROR: --out must be a real file path because a companion .yml file is always written.", file=sys.stderr)
        sys.exit(2)

    scan_mode = args.scan_by is not None
    try:
        n_min, percents, kpercents, scan_filter, _auto_adjustments = resolve_effective_runtime_args(
            args, argv
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    bam = pysam.AlignmentFile(args.bam, "rb")
    lengths = dict(zip(bam.references, bam.lengths))

    idx_table = get_idxstats_table(args.bam)
    mapped_reads = {c: v[1] for c, v in idx_table.items() if c != "*"}

    flag_filter = build_flag_filter(
        include_duplicates=args.include_duplicates,
        include_secondary=scan_mode,
    )

    if scan_mode:
        virus_ids = [c for c, length in lengths.items() if c != "*" and length > 0]
        id_to_name = {vid: vid for vid in virus_ids}
        has_two_cols = False
        name_to_ids = None
    else:
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

    out_path = Path(args.out).expanduser().resolve()
    yml_path = out_path.with_suffix(".yml")
    print(f"Writing output to\n  {out_path}", file=sys.stderr)
    print(f"Writing run arguments to\n  {yml_path}", file=sys.stderr)
    write_yaml_file(
        yml_path,
        build_run_metadata(args, out_path, n_min, percents, kpercents, scan_mode),
    )
    out_fh = open(out_path, "w", encoding="utf-8")

    header = build_header(n_min, percents, kpercents)
    print("\t".join(header), file=out_fh)

    # Collect per-contig depth histograms and stats to support virus-level pooling
    contig_depth_counts = {}  # vid -> Counter(depth -> npos)
    contig_stats = {}  # vid -> dict

    for vid in virus_ids:
        result = compute_contig_result(
            bam=bam,
            contig=vid,
            length=lengths[vid],
            mapped_reads=mapped_reads,
            min_mapq=args.min_mapq,
            flag_filter=flag_filter,
            max_depth=args.max_depth,
            n_min=n_min,
            percents=percents,
            kpercents=kpercents,
        )

        if scan_mode:
            if passes_scan_filter(result["metrics"], scan_filter):
                row = build_contig_row(vid, vid, result, n_min, percents, kpercents)
                print("\t".join(row), file=out_fh)
            continue

        contig_depth_counts[vid] = result["depth_counts"]
        contig_stats[vid] = {
            "L": result["length"],
            "m": result["mapped_reads"],
            "mapped_per_bp": result["mapped_per_bp"],
            "breadth_gt0": result["breadth_gt0"],
            "mean_all": result["mean_all"],
            "frac_ge": result["frac_ge"],
            "mean_ge": result["mean_ge"],
            "median_ge": result["median_ge"],
        }

        row = build_contig_row(vid, id_to_name.get(vid, vid), result, n_min, percents, kpercents)
        print("\t".join(row), file=out_fh)

    # Virus-level rows (only if 2-column list provided)
    if not scan_mode and has_two_cols and name_to_ids:
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

            virus_frac_ge, virus_mean_ge, virus_median_ge = mean_median_depth_ge_n(pooled, n_min)

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

    out_fh.close()
    bam.close()


if __name__ == "__main__":
    main()
