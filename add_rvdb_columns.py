#!/usr/bin/env python3
import json
import gzip
import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Append RVDB membership and accession-list columns to a Kraken2 "
            "inspect report while preserving the original order and content."
        )
    )
    parser.add_argument(
        "--kraken-inspect",
        required=True,
        help="Path to kraken2-inspect output file, e.g. viral.inspect.report",
    )
    parser.add_argument(
        "--rvdb-json",
        required=True,
        help="Path to RVDB JSON or JSON.gz file",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output TSV file",
    )
    return parser.parse_args()


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def load_rvdb_taxid_to_accessions(rvdb_json_path):
    taxid_to_accessions = defaultdict(list)

    with open_maybe_gzip(rvdb_json_path) as f:
        rvdb_data = json.load(f)

    for accession, meta in rvdb_data.items():
        if not isinstance(meta, dict):
            continue
        taxid = meta.get("ncbitaxon")
        if taxid is None:
            continue
        taxid = str(taxid).strip()
        if taxid and taxid.isdigit():
            taxid_to_accessions[taxid].append(str(accession).strip())

    for taxid in taxid_to_accessions:
        taxid_to_accessions[taxid] = sorted(set(taxid_to_accessions[taxid]))

    return taxid_to_accessions


def main():
    args = parse_args()
    taxid_to_accessions = load_rvdb_taxid_to_accessions(args.rvdb_json)

    with open(args.kraken_inspect, "r", encoding="utf-8") as fin, \
         open(args.output, "w", encoding="utf-8") as fout:

        # Write header
        fout.write(
            "percent\tclade_count\tdirect_count\trank\ttaxid\tname\tin_rvdb\trvdb_accessions\n"
        )

        for line in fin:
            line = line.rstrip("\n")

            if not line.strip():
                continue

            parts = line.split("\t")

            # Kraken2 inspect/report format:
            # 1 percent
            # 2 clade count
            # 3 direct count
            # 4 rank
            # 5 taxid
            # 6 name
            if len(parts) >= 6:
                taxid = parts[4].strip()
                if taxid.isdigit():
                    accessions = taxid_to_accessions.get(taxid, [])
                    in_rvdb = "1" if accessions else "0"
                    accession_list = ",".join(accessions)
                    fout.write(f"{line}\t{in_rvdb}\t{accession_list}\n")
                    continue

            # Fallback for unexpected lines
            fout.write(f"{line}\tNA\tNA\n")


if __name__ == "__main__":
    main()

