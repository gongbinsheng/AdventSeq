#!/usr/bin/env python3

import argparse
import gzip
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

from contig_info_utils import load_contig_info


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "AdVentSeq"


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a gzipped JSON mapping from accession to NCBI taxonomy ID, "
            "NCBI taxonomy name, and organism."
        )
    )
    parser.add_argument(
        "--contig_info",
        required=True,
        type=Path,
        help="Contig metadata input in JSON or SQLite format.",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email address required by NCBI E-utilities.",
    )
    parser.add_argument(
        "--out",
        default=Path("cbi_taxonomy_map.json.gz"),
        type=Path,
        help="Output gzipped JSON file. Defaults to cbi_taxonomy_map.json.gz.",
    )
    return parser.parse_args(argv)


def chunked(values, size):
    for index in range(0, len(values), size):
        yield values[index:index + size]


def build_eutils_url(endpoint: str, params: dict[str, str]) -> str:
    query = {
        "tool": TOOL_NAME,
        **params,
    }
    return f"{EUTILS_BASE}/{endpoint}?{urlencode(query)}"


def open_xml(url: str, opener=urlopen):
    with opener(url) as response:
        payload = response.read()
    return ET.fromstring(payload)


def _docsum_item_text(docsum, item_name):
    for item in docsum.findall("./Item"):
        if item.get("Name") == item_name:
            return (item.text or "").strip() or None
    return None


def _match_requested_accession(candidate, requested_batch):
    if not candidate:
        return None

    if candidate in requested_batch:
        return candidate

    candidate_stem = candidate.split(".", 1)[0]
    for accession in requested_batch:
        if accession.split(".", 1)[0] == candidate_stem:
            return accession

    return None


def fetch_nuccore_taxids(accessions, email, batch_size=200, opener=urlopen):
    accession_to_taxid = {}

    for batch in chunked(accessions, batch_size):
        url = build_eutils_url(
            "esummary.fcgi",
            {
                "db": "nuccore",
                "id": ",".join(batch),
                "retmode": "xml",
                "email": email,
            },
        )
        root = open_xml(url, opener=opener)

        for docsum in root.findall("./DocSum"):
            accession = _match_requested_accession(
                _docsum_item_text(docsum, "AccessionVersion") or _docsum_item_text(docsum, "Caption"),
                batch,
            )
            if not accession:
                continue

            taxid = _docsum_item_text(docsum, "TaxId")
            accession_to_taxid[accession] = taxid

    return accession_to_taxid


def fetch_taxonomy_names(taxids, email, batch_size=200, opener=urlopen):
    taxid_to_name = {}
    normalized_taxids = [str(taxid) for taxid in taxids if taxid]

    for batch in chunked(normalized_taxids, batch_size):
        url = build_eutils_url(
            "efetch.fcgi",
            {
                "db": "taxonomy",
                "id": ",".join(batch),
                "retmode": "xml",
                "email": email,
            },
        )
        root = open_xml(url, opener=opener)

        for taxon in root.findall(".//Taxon"):
            taxid = taxon.findtext("TaxId")
            scientific_name = taxon.findtext("ScientificName")
            if taxid and scientific_name:
                taxid_to_name[str(taxid)] = scientific_name.strip()

    return taxid_to_name


def build_taxonomy_map(contig_info, email, batch_size=200, opener=urlopen):
    accessions = sorted(contig_info.keys())
    accession_to_taxid = fetch_nuccore_taxids(
        accessions,
        email=email,
        batch_size=batch_size,
        opener=opener,
    )
    taxid_to_name = fetch_taxonomy_names(
        sorted({taxid for taxid in accession_to_taxid.values() if taxid}),
        email=email,
        batch_size=batch_size,
        opener=opener,
    )

    taxonomy_map = {}
    for accession in accessions:
        organism = contig_info[accession].get("organism")
        taxid = accession_to_taxid.get(accession)
        taxonomy_map[accession] = {
            "ncbitaxon": taxid if taxid else None,
            "ncbitaxonname": taxid_to_name.get(taxid) if taxid else None,
            "organism": organism,
        }

    return taxonomy_map


def write_taxonomy_map(taxonomy_map, output_path: Path) -> None:
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        json.dump(taxonomy_map, handle, ensure_ascii=False, sort_keys=True)


def main(argv=None) -> None:
    args = parse_args(argv)
    contig_info = load_contig_info(args.contig_info)
    taxonomy_map = build_taxonomy_map(contig_info, email=args.email)
    write_taxonomy_map(taxonomy_map, args.out)
    print(args.out)


if __name__ == "__main__":
    main()
