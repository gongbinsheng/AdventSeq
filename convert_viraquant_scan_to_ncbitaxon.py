#!/usr/bin/env python3
"""
Convert ViraQuant scan summaries from contig level to ncbitaxon level.

This script aggregates contig rows from a ViraQuant_scan TSV/TSV.GZ file using an
accession-keyed NCBI taxonomy mapping. Exact pooled metrics are recomputed where
possible from the existing summary columns, while quantile- and median-like
metrics use documented TSV-only approximations.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path
import re

from contig_info_utils import load_json_mapping


NA_TOKEN = "NA"
REQUIRED_BASE_COLUMNS = {
    "level",
    "virus",
    "seq_id",
    "n_contigs",
    "length",
    "mapped_reads",
    "mapped_per_bp",
    "breadth_cov_gt0",
    "mean_depth_all",
    "best_contig_id",
    "best_contig_mapped_per_bp",
}

FRAC_RE = re.compile(r"^frac_ge_(\d+)$")
MEAN_GE_RE = re.compile(r"^mean_depth_ge_(\d+)$")
MEDIAN_GE_RE = re.compile(r"^median_depth_ge_(\d+)$")
PASS_RE = re.compile(r"^pass_k_(\d+)pct_ge_(\d+)$")
DEPTH_AT_RE = re.compile(r"^depth_at_(\d+)pct$")
BEST_FRAC_RE = re.compile(r"^best_contig_frac_ge_(\d+)$")


def get_entrez_id(seq_id: str) -> str:
    parts = seq_id.strip().split("|")
    if len(parts) > 2:
        return parts[2]
    return parts[0]


def open_text(path: str | Path, mode: str):
    path = Path(path)
    if "b" in mode:
        raise ValueError("Binary mode is not supported")
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def parse_int(value: str | None) -> int | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.upper() == NA_TOKEN:
        return None
    return int(value)


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.upper() == NA_TOKEN:
        return None
    return float(value)


def format_int(value: int | None) -> str:
    return NA_TOKEN if value is None else str(value)


def format_fixed(value: float | None) -> str:
    return NA_TOKEN if value is None else f"{value:.6f}"


def format_compact(value: float | None) -> str:
    return NA_TOKEN if value is None else f"{value:.6g}"


def format_pass(value: int | None) -> str:
    return NA_TOKEN if value is None else str(value)


def weighted_mean(pairs: list[tuple[float, float]]) -> float | None:
    usable = [(value, weight) for value, weight in pairs if value is not None and weight > 0]
    if not usable:
        return None
    total_weight = sum(weight for _, weight in usable)
    if total_weight <= 0:
        return None
    return sum(value * weight for value, weight in usable) / total_weight


def weighted_median(pairs: list[tuple[float, float]]) -> float | None:
    usable = sorted((value, weight) for value, weight in pairs if value is not None and weight > 0)
    if not usable:
        return None
    total_weight = sum(weight for _, weight in usable)
    if total_weight <= 0:
        return None
    threshold = total_weight / 2.0
    cumulative = 0.0
    for value, weight in usable:
        cumulative += weight
        if cumulative >= threshold:
            return value
    return usable[-1][0]


def build_output_header(fieldnames: list[str]) -> list[str]:
    if "ncbitaxonname" in fieldnames:
        return list(fieldnames)
    insert_after = fieldnames.index("virus") + 1
    return fieldnames[:insert_after] + ["ncbitaxonname"] + fieldnames[insert_after:]


def validate_header(fieldnames: list[str]) -> None:
    missing = sorted(REQUIRED_BASE_COLUMNS - set(fieldnames))
    if missing:
        raise ValueError(
            "Input is not a ViraQuant_scan-style table; missing required columns: "
            + ", ".join(missing)
        )

    frac_ns = {int(match.group(1)) for name in fieldnames if (match := FRAC_RE.match(name))}
    for name in fieldnames:
        if match := PASS_RE.match(name):
            n_value = int(match.group(2))
            if n_value not in frac_ns:
                raise ValueError(f"Column '{name}' requires matching 'frac_ge_{n_value}'")
        if match := MEAN_GE_RE.match(name):
            n_value = int(match.group(1))
            if n_value not in frac_ns:
                raise ValueError(f"Column '{name}' requires matching 'frac_ge_{n_value}'")
        if match := MEDIAN_GE_RE.match(name):
            n_value = int(match.group(1))
            if n_value not in frac_ns:
                raise ValueError(f"Column '{name}' requires matching 'frac_ge_{n_value}'")


def discover_metric_columns(fieldnames: list[str]) -> dict[str, dict | list]:
    return {
        "frac": {int(match.group(1)): name for name in fieldnames if (match := FRAC_RE.match(name))},
        "mean_ge": {int(match.group(1)): name for name in fieldnames if (match := MEAN_GE_RE.match(name))},
        "median_ge": {int(match.group(1)): name for name in fieldnames if (match := MEDIAN_GE_RE.match(name))},
        "pass": [
            (name, int(match.group(1)), int(match.group(2)))
            for name in fieldnames
            if (match := PASS_RE.match(name))
        ],
        "depth_at": [
            (name, int(match.group(1)))
            for name in fieldnames
            if (match := DEPTH_AT_RE.match(name))
        ],
        "best_frac": {
            int(match.group(1)): name for name in fieldnames if (match := BEST_FRAC_RE.match(name))
        },
    }


def choose_best_contig(rows: list[dict]) -> dict:
    return max(
        rows,
        key=lambda item: (
            item["mapped_per_bp_num"],
            item["mapped_reads_num"],
            item["length_num"],
            -item["input_order"],
        ),
    )


def resolve_ncbitaxon(seq_id: str, taxonomy_map: dict) -> tuple[str | None, str | None]:
    accession = get_entrez_id(seq_id)
    record = taxonomy_map.get(accession, {})
    ncbitaxon = record.get("ncbitaxon")
    if ncbitaxon is not None:
        ncbitaxon = str(ncbitaxon).strip()
    if not ncbitaxon:
        return None, None

    tax_name = record.get("ncbitaxonname")
    if tax_name is not None:
        tax_name = str(tax_name).strip()
    return ncbitaxon, (tax_name or None)


def load_scan_rows(input_path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(input_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input scan file is empty or missing a header")
        fieldnames = list(reader.fieldnames)
        validate_header(fieldnames)
        rows = [row for row in reader if row.get("level") == "contig"]
    return fieldnames, rows


def enrich_row(row: dict[str, str], input_order: int) -> dict:
    enriched = dict(row)
    enriched["input_order"] = input_order
    enriched["length_num"] = parse_int(row.get("length")) or 0
    enriched["mapped_reads_num"] = parse_int(row.get("mapped_reads")) or 0
    enriched["mapped_per_bp_num"] = parse_float(row.get("mapped_per_bp")) or 0.0
    return enriched


def aggregate_group(
    ncbitaxon: str,
    ncbitaxonname: str | None,
    rows: list[dict],
    output_header: list[str],
    metrics: dict[str, dict | list],
) -> dict[str, str]:
    best_row = choose_best_contig(rows)
    total_length = sum(row["length_num"] for row in rows)
    total_mapped_reads = sum(row["mapped_reads_num"] for row in rows)
    group_frac: dict[int, float | None] = {}

    result = {column: NA_TOKEN for column in output_header}
    result["level"] = "ncbitaxon"
    result["virus"] = ncbitaxon
    result["ncbitaxonname"] = ncbitaxonname or NA_TOKEN
    result["seq_id"] = ",".join(row["seq_id"] for row in rows)
    result["n_contigs"] = str(len(rows))
    result["length"] = str(total_length)
    result["mapped_reads"] = str(total_mapped_reads)
    result["mapped_per_bp"] = format_compact(
        (total_mapped_reads / total_length) if total_length > 0 else 0.0
    )
    result["breadth_cov_gt0"] = format_fixed(
        weighted_mean(
            [(parse_float(row.get("breadth_cov_gt0")), row["length_num"]) for row in rows]
        )
    )
    result["mean_depth_all"] = format_fixed(
        weighted_mean(
            [(parse_float(row.get("mean_depth_all")), row["length_num"]) for row in rows]
        )
    )

    for n_value, column in metrics["frac"].items():
        frac_value = weighted_mean(
            [(parse_float(row.get(column)), row["length_num"]) for row in rows]
        )
        group_frac[n_value] = frac_value
        result[column] = format_fixed(frac_value)

    for n_value, column in metrics["mean_ge"].items():
        result[column] = format_fixed(
            weighted_mean(
                [
                    (
                        parse_float(row.get(column)),
                        row["length_num"] * (parse_float(row.get(metrics["frac"][n_value])) or 0.0),
                    )
                    for row in rows
                ]
            )
        )

    for n_value, column in metrics["median_ge"].items():
        result[column] = format_fixed(
            weighted_median(
                [
                    (
                        parse_float(row.get(column)),
                        row["length_num"] * (parse_float(row.get(metrics["frac"][n_value])) or 0.0),
                    )
                    for row in rows
                ]
            )
        )

    for column, k_percent, n_value in metrics["pass"]:
        frac_value = group_frac.get(n_value)
        if frac_value is None:
            result[column] = NA_TOKEN
        else:
            result[column] = format_pass(1 if frac_value >= (k_percent / 100.0) else 0)

    for column, _percent in metrics["depth_at"]:
        result[column] = format_int(
            weighted_median(
                [(parse_int(row.get(column)), row["length_num"]) for row in rows]
            )
        )

    result["best_contig_id"] = best_row["seq_id"]
    result["best_contig_mapped_per_bp"] = format_compact(best_row["mapped_per_bp_num"])

    if "best_contig_depth_at_90pct" in result:
        best_depth_90 = parse_int(best_row.get("depth_at_90pct"))
        if best_depth_90 is None:
            best_depth_90 = parse_int(best_row.get("best_contig_depth_at_90pct"))
        result["best_contig_depth_at_90pct"] = format_int(best_depth_90)

    for n_value, column in metrics["best_frac"].items():
        best_frac = parse_float(best_row.get(metrics["frac"].get(n_value, "")))
        if best_frac is None:
            best_frac = parse_float(best_row.get(column))
        result[column] = format_fixed(best_frac)

    return result


def convert_scan(input_path: str | Path, taxonomy_map_path: str | Path, output_path: str | Path) -> None:
    taxonomy_map = load_json_mapping(taxonomy_map_path)
    fieldnames, contig_rows = load_scan_rows(input_path)
    output_header = build_output_header(fieldnames)
    metrics = discover_metric_columns(fieldnames)

    grouped_rows: dict[str, list[dict]] = defaultdict(list)
    group_tax_names: dict[str, str | None] = {}
    unresolved_rows: list[dict[str, str]] = []
    warned_names: set[str] = set()

    for input_order, row in enumerate(contig_rows):
        enriched = enrich_row(row, input_order)
        ncbitaxon, ncbitaxonname = resolve_ncbitaxon(enriched["seq_id"], taxonomy_map)
        if ncbitaxon is None:
            passthrough = dict(row)
            passthrough["level"] = "contig"
            passthrough["ncbitaxonname"] = NA_TOKEN
            unresolved_rows.append(passthrough)
            continue

        grouped_rows[ncbitaxon].append(enriched)
        current_name = group_tax_names.get(ncbitaxon)
        if current_name is None and ncbitaxonname:
            group_tax_names[ncbitaxon] = ncbitaxonname
        elif current_name and ncbitaxonname and current_name != ncbitaxonname and ncbitaxon not in warned_names:
            print(
                f"WARNING: conflicting ncbitaxonname values for ncbitaxon {ncbitaxon}: "
                f"keeping '{current_name}', ignoring '{ncbitaxonname}'",
                file=sys.stderr,
            )
            warned_names.add(ncbitaxon)
        elif ncbitaxon not in group_tax_names:
            group_tax_names[ncbitaxon] = None

    output_rows = [
        aggregate_group(ncbitaxon, group_tax_names.get(ncbitaxon), rows, output_header, metrics)
        for ncbitaxon, rows in grouped_rows.items()
    ]
    output_rows.extend(unresolved_rows)

    with open_text(output_path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a ViraQuant_scan TSV from contig summaries to ncbitaxon-level grouped summaries."
    )
    parser.add_argument("--input", required=True, help="Input ViraQuant_scan TSV or TSV.GZ file")
    parser.add_argument(
        "--taxonomy-map",
        required=True,
        dest="taxonomy_map",
        help="Accession-keyed taxonomy map in JSON or JSON.GZ format",
    )
    parser.add_argument("--output", required=True, help="Output TSV or TSV.GZ path")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        convert_scan(args.input, args.taxonomy_map, args.output)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
