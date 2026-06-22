#!/usr/bin/env python3

import argparse
import gzip
import hashlib
import json
import shutil
import time
from http.client import HTTPException, IncompleteRead
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

from AdVentSeq.contig_info_utils import load_contig_info, load_json_mapping

try:
    from tqdm import tqdm
except ModuleNotFoundError:  # pragma: no cover - fallback for minimal environments
    def tqdm(iterable=None, **kwargs):
        return iterable


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
        default=Path("ncbi_taxonomy_map.json.gz"),
        type=Path,
        help="Output gzipped JSON file. Defaults to ncbi_taxonomy_map.json.gz.",
    )
    parser.add_argument(
        "--retry",
        default=3,
        type=int,
        help="Number of retries for each failed NCBI request. Defaults to 3.",
    )
    parser.add_argument(
        "--prev_map",
        default=None,
        type=Path,
        help=(
            "Optional previous taxonomy map (.json or .json.gz). Accessions already "
            "resolved there are reused, so only new/unresolved accessions are queried "
            "against NCBI."
        ),
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


def open_xml_with_retry(url: str, retry: int, opener=urlopen, sleep_base_seconds: float = 1.0):
    attempts = retry + 1
    last_error = None

    for attempt in range(1, attempts + 1):
        try:
            return open_xml(url, opener=opener)
        except (IncompleteRead, HTTPException, URLError, OSError, ET.ParseError, ValueError) as exc:
            last_error = exc
            if attempt == attempts:
                break
            time.sleep(sleep_base_seconds * attempt)

    raise RuntimeError(f"Failed to fetch NCBI response after {attempts} attempts: {url}") from last_error


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


def get_cache_dir(output_path: Path) -> Path:
    name = output_path.name
    if name.endswith(".json.gz"):
        base_name = name[:-8]
    else:
        base_name = output_path.stem
    return output_path.parent / f"{base_name}.cache"


def get_batch_cache_path(cache_dir: Path, prefix: str, batch) -> Path:
    batch_key = hashlib.sha1("\n".join(str(item) for item in batch).encode("utf-8")).hexdigest()
    return cache_dir / prefix / f"{batch_key}.json"


def load_cached_batch(cache_path: Path):
    if not cache_path.exists():
        return None

    with cache_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_cached_batch(cache_path: Path, data) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, sort_keys=True)


def fetch_nuccore_batch(batch, email, retry, opener):
    url = build_eutils_url(
        "esummary.fcgi",
        {
            "db": "nuccore",
            "id": ",".join(batch),
            "retmode": "xml",
            "email": email,
        },
    )
    root = open_xml_with_retry(url, retry=retry, opener=opener)
    batch_result = {accession: None for accession in batch}

    for docsum in root.findall("./DocSum"):
        accession = _match_requested_accession(
            _docsum_item_text(docsum, "AccessionVersion") or _docsum_item_text(docsum, "Caption"),
            batch,
        )
        if not accession:
            continue

        batch_result[accession] = _docsum_item_text(docsum, "TaxId")

    return batch_result


def fetch_taxonomy_batch(batch, email, retry, opener):
    url = build_eutils_url(
        "efetch.fcgi",
        {
            "db": "taxonomy",
            "id": ",".join(batch),
            "retmode": "xml",
            "email": email,
        },
    )
    root = open_xml_with_retry(url, retry=retry, opener=opener)
    batch_result = {}

    for taxon in root.findall(".//Taxon"):
        taxid = taxon.findtext("TaxId")
        scientific_name = taxon.findtext("ScientificName")
        if taxid and scientific_name:
            batch_result[str(taxid)] = scientific_name.strip()

    return batch_result


def fetch_nuccore_taxids(accessions, email, batch_size=200, retry=3, cache_dir=None, opener=urlopen):
    accession_to_taxid = {}
    batches = list(chunked(accessions, batch_size))

    for batch in tqdm(batches, desc="NCBI nuccore batches", unit="batch"):
        cache_path = get_batch_cache_path(cache_dir, "nuccore", batch) if cache_dir else None
        batch_result = load_cached_batch(cache_path) if cache_path else None
        if batch_result is None:
            batch_result = fetch_nuccore_batch(batch, email=email, retry=retry, opener=opener)
            if cache_path:
                write_cached_batch(cache_path, batch_result)
        accession_to_taxid.update(batch_result)

    return accession_to_taxid


def fetch_taxonomy_names(taxids, email, batch_size=200, retry=3, cache_dir=None, opener=urlopen):
    taxid_to_name = {}
    normalized_taxids = [str(taxid) for taxid in taxids if taxid]
    batches = list(chunked(normalized_taxids, batch_size))

    for batch in tqdm(batches, desc="NCBI taxonomy batches", unit="batch"):
        cache_path = get_batch_cache_path(cache_dir, "taxonomy", batch) if cache_dir else None
        batch_result = load_cached_batch(cache_path) if cache_path else None
        if batch_result is None:
            batch_result = fetch_taxonomy_batch(batch, email=email, retry=retry, opener=opener)
            if cache_path:
                write_cached_batch(cache_path, batch_result)
        taxid_to_name.update(batch_result)

    return taxid_to_name


def build_taxonomy_map(contig_info, email, batch_size=200, retry=3, cache_dir=None, opener=urlopen, prev_map=None):
    prev_map = prev_map or {}
    accessions = sorted(contig_info.keys())

    # Reuse non-null entries from the previous map; re-query everything else
    # (new accessions and previously unresolved ones).
    reused = {}
    to_lookup = []
    for accession in accessions:
        prev_entry = prev_map.get(accession)
        if prev_entry and prev_entry.get("ncbitaxon"):
            reused[accession] = (prev_entry["ncbitaxon"], prev_entry.get("ncbitaxonname"))
        else:
            to_lookup.append(accession)

    accession_to_taxid = fetch_nuccore_taxids(
        to_lookup,
        email=email,
        batch_size=batch_size,
        retry=retry,
        cache_dir=cache_dir,
        opener=opener,
    )

    # Names already known from the previous map; only query taxonomy for taxids
    # we have not seen a name for before.
    known_names = {
        str(entry["ncbitaxon"]): entry["ncbitaxonname"]
        for entry in prev_map.values()
        if entry.get("ncbitaxon") and entry.get("ncbitaxonname")
    }
    taxids_needing_names = sorted(
        {str(taxid) for taxid in accession_to_taxid.values() if taxid and str(taxid) not in known_names}
    )
    name_lookup = {
        **known_names,
        **fetch_taxonomy_names(
            taxids_needing_names,
            email=email,
            batch_size=batch_size,
            retry=retry,
            cache_dir=cache_dir,
            opener=opener,
        ),
    }

    taxonomy_map = {}
    for accession in accessions:
        organism = contig_info[accession].get("organism")
        if accession in reused:
            taxid, name = reused[accession]
            taxonomy_map[accession] = {
                "ncbitaxon": taxid,
                "ncbitaxonname": name,
                "organism": organism,
            }
            continue

        taxid = accession_to_taxid.get(accession)
        taxonomy_map[accession] = {
            "ncbitaxon": taxid if taxid else None,
            "ncbitaxonname": name_lookup.get(str(taxid)) if taxid else None,
            "organism": organism,
        }

    return taxonomy_map


def write_taxonomy_map(taxonomy_map, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        json.dump(taxonomy_map, handle, ensure_ascii=False, sort_keys=True)


def cleanup_cache_dir(cache_dir: Path) -> None:
    shutil.rmtree(cache_dir, ignore_errors=True)


def main(argv=None) -> None:
    args = parse_args(argv)
    if args.retry < 0:
        raise SystemExit("--retry must be 0 or greater")

    contig_info = load_contig_info(args.contig_info)
    prev_map = load_json_mapping(args.prev_map) if args.prev_map else None
    cache_dir = get_cache_dir(args.out)
    taxonomy_map = build_taxonomy_map(
        contig_info,
        email=args.email,
        retry=args.retry,
        cache_dir=cache_dir,
        prev_map=prev_map,
    )
    write_taxonomy_map(taxonomy_map, args.out)
    cleanup_cache_dir(cache_dir)
    print(args.out)


if __name__ == "__main__":
    main()
