import gzip
import json
import shutil
import sqlite3
import tempfile
from pathlib import Path


def normalize_record(columns, row):
    record = dict(zip(columns, row, strict=True))
    accession = record.pop("accs")

    seqlen = record.get("seqlen")
    if isinstance(seqlen, str):
        try:
            record["seqlen"] = int(seqlen)
        except ValueError:
            pass

    return accession, record


def load_contig_info_from_json(contig_info_path: Path):
    suffixes = contig_info_path.suffixes
    if suffixes[-2:] == [".json", ".gz"]:
        with gzip.open(contig_info_path, "rt", encoding="utf-8") as handle:
            return json.load(handle)

    with contig_info_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_contig_info_from_sqlite_db(db_path: Path):
    contig_info = {}
    with sqlite3.connect(db_path) as connection:
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM rvdb")
        columns = [description[0] for description in cursor.description]

        for row in cursor:
            accession, record = normalize_record(columns, row)
            contig_info[accession] = record

    return contig_info


def load_contig_info_from_sqlite(contig_info_path: Path):
    suffixes = contig_info_path.suffixes
    if suffixes[-1:] == [".gz"]:
        with tempfile.NamedTemporaryFile(suffix=".sqlite.db", delete=False) as tmp_handle:
            tmp_path = Path(tmp_handle.name)

        try:
            with gzip.open(contig_info_path, "rb") as source, tmp_path.open("wb") as destination:
                shutil.copyfileobj(source, destination)
            return load_contig_info_from_sqlite_db(tmp_path)
        finally:
            tmp_path.unlink(missing_ok=True)

    return load_contig_info_from_sqlite_db(contig_info_path)


def load_contig_info(contig_info_filepath):
    contig_info_path = Path(contig_info_filepath)
    suffixes = contig_info_path.suffixes

    if suffixes[-2:] == [".json", ".gz"] or suffixes[-1:] == [".json"]:
        return load_contig_info_from_json(contig_info_path)

    if suffixes[-3:] == [".sqlite", ".db", ".gz"] or suffixes[-2:] == [".db", ".gz"]:
        return load_contig_info_from_sqlite(contig_info_path)

    if suffixes[-2:] == [".sqlite", ".db"] or suffixes[-1:] == [".db"]:
        return load_contig_info_from_sqlite(contig_info_path)

    raise ValueError(
        f"Unsupported file type for {contig_info_filepath}. "
        "Supported formats: .json, .json.gz, .db, .db.gz, .sqlite.db, .sqlite.db.gz"
    )


def load_json_mapping(mapping_filepath):
    mapping_path = Path(mapping_filepath)
    suffixes = mapping_path.suffixes

    if suffixes[-2:] == [".json", ".gz"]:
        with gzip.open(mapping_path, "rt", encoding="utf-8") as handle:
            return json.load(handle)

    if suffixes[-1:] == [".json"]:
        with mapping_path.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    raise ValueError(
        f"Unsupported mapping file type for {mapping_filepath}. "
        "Supported formats: .json, .json.gz"
    )
