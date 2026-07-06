#!/usr/bin/env python3

import argparse
import shutil
import gzip
import json
import sqlite3
from pathlib import Path

from tqdm import tqdm
from AdventSeq.contig_info_utils import normalize_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert an RVDB SQLite database into a gzipped JSON mapping of "
            "accession to contig metadata. "
            "RVDB: https://rvdb.dbi.udel.edu"
        )
    )
    parser.add_argument(
        "--RVDB_ROOT",
        dest="rvdb_root",
        type=Path,
        default=Path.cwd(),
        help="RVDB root directory that contains the release folders. Defaults to the current directory.",
    )
    parser.add_argument(
        "release",
        help='RVDB release directory to process, for example "v31.0".',
    )
    return parser.parse_args()

def ensure_db_path(db_path: Path) -> None:
    if db_path.exists():
        return

    gz_path = Path(f"{db_path}.gz")
    if not gz_path.exists():
        raise FileNotFoundError(f"SQLite database not found: {db_path}")

    # Materialize the SQLite file only when needed so we can work from either
    # the plain database or the packaged .gz archive without deleting the input.
    with gzip.open(gz_path, "rb") as source, db_path.open("wb") as destination:
        shutil.copyfileobj(source, destination)


def convert_release(rvdb_root: Path, release: str) -> Path:
    release_dir = rvdb_root / release
    db_path = release_dir / f"U-RVDB{release}.sqlite.db"
    output_path = release_dir / "contig_info.json.gz"

    release_dir.mkdir(parents=True, exist_ok=True)
    ensure_db_path(db_path)

    with sqlite3.connect(db_path) as connection:
        connection.row_factory = None
        count_cursor = connection.cursor()
        count_cursor.execute("SELECT COUNT(*) FROM rvdb")
        total_rows = count_cursor.fetchone()[0]

        cursor = connection.cursor()
        cursor.execute("SELECT * FROM rvdb")
        columns = [description[0] for description in cursor.description]

        with gzip.open(output_path, "wt", encoding="utf-8") as handle:
            handle.write("{\n")

            first = True
            with tqdm(
                total=total_rows,
                unit="row",
                unit_scale=True,
                desc=f"Converting {release}",
            ) as progress:
                # Stream rows directly from SQLite into the gzip output so the
                # conversion does not need to hold the full dataset in memory.
                for row in cursor:
                    accession, record = normalize_record(columns, row)
                    if not first:
                        handle.write(",\n")
                    first = False

                    # Emit one top-level JSON entry per accession, matching the
                    # contig-info layout expected by downstream tooling.
                    key_json = json.dumps(accession, ensure_ascii=False)
                    record_json = json.dumps(record, ensure_ascii=False, sort_keys=True)
                    handle.write(f"  {key_json}: {record_json}")
                    progress.update()

            handle.write("\n}\n")

    return output_path


def main() -> None:
    args = parse_args()
    output_path = convert_release(args.rvdb_root, args.release)
    print(output_path)


if __name__ == "__main__":
    main()
