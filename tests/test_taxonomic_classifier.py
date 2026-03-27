import gzip
import json
import sqlite3
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

from contig_info_utils import load_contig_info
from taxonomic_classifier import main, resolve_virus_group, write_read_count_report


class ContigInfoLoaderTests(unittest.TestCase):
    def test_load_contig_info_from_json_gz(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            contig_path = Path(tmpdir) / "contig_info.json.gz"
            expected = {
                "HM118273.1": {
                    "source": "GENBANK",
                    "description": "Example",
                    "seqlen": 309,
                    "organism": "Human immunodeficiency virus 1",
                }
            }

            with gzip.open(contig_path, "wt", encoding="utf-8") as handle:
                json.dump(expected, handle)

            self.assertEqual(load_contig_info(contig_path), expected)

    def test_load_contig_info_from_sqlite_db(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "contig_info.sqlite.db"
            with sqlite3.connect(db_path) as connection:
                connection.execute(
                    "CREATE TABLE rvdb (accs TEXT, source TEXT, description TEXT, seqlen TEXT, organism TEXT)"
                )
                connection.execute(
                    "INSERT INTO rvdb VALUES (?, ?, ?, ?, ?)",
                    ("HM118273.1", "GENBANK", "Example", "309", "Human immunodeficiency virus 1"),
                )
                connection.commit()

            loaded = load_contig_info(db_path)
            self.assertEqual(loaded["HM118273.1"]["seqlen"], 309)
            self.assertEqual(loaded["HM118273.1"]["organism"], "Human immunodeficiency virus 1")


class VirusGroupingTests(unittest.TestCase):
    def setUp(self):
        self.contig_info = {
            "ACC1.1": {"organism": "Alpha virus"},
            "ACC2.1": {"organism": "Beta virus"},
            "ACC3.1": {"organism": "Shared Name"},
            "ACC4.1": {"organism": "Shared Name"},
        }
        self.taxonomy_map = {
            "ACC1.1": {
                "ncbitaxon": "111",
                "ncbitaxonname": "Taxon Alpha",
                "organism": "Alpha virus",
            },
            "ACC2.1": {
                "ncbitaxon": "111",
                "ncbitaxonname": "Taxon Alpha",
                "organism": "Beta virus",
            },
            "ACC3.1": {
                "ncbitaxon": "222",
                "ncbitaxonname": "Shared Name",
                "organism": "Shared Name",
            },
            "ACC4.1": {
                "ncbitaxon": None,
                "ncbitaxonname": None,
                "organism": "Shared Name",
            },
        }

    def test_resolve_organism_mode(self):
        group = resolve_virus_group("ACC1.1", self.contig_info, self.taxonomy_map, "organism")
        self.assertEqual(group.label, "Alpha virus")
        self.assertFalse(group.is_fallback)
        self.assertEqual(group.comparison_key, ("organism", "Alpha virus"))

    def test_same_taxid_collapses_even_with_different_organisms(self):
        first = resolve_virus_group("ACC1.1", self.contig_info, self.taxonomy_map, "ncbitaxon")
        second = resolve_virus_group("ACC2.1", self.contig_info, self.taxonomy_map, "ncbitaxon")

        self.assertEqual(first.label, "111")
        self.assertEqual(first.comparison_key, second.comparison_key)
        self.assertFalse(first.is_fallback)
        self.assertFalse(second.is_fallback)

    def test_missing_taxonomy_falls_back_to_organism(self):
        group = resolve_virus_group("ACC4.1", self.contig_info, self.taxonomy_map, "ncbitaxonname")
        self.assertEqual(group.label, "Shared Name")
        self.assertTrue(group.is_fallback)
        self.assertEqual(group.comparison_key, ("fallback_organism", "Shared Name"))

    def test_fallback_and_true_taxonomy_name_stay_separate_in_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_id = str(Path(tmpdir) / "sample")
            read_count = defaultdict(int, {
                "total": 10,
                "mapped": 6,
                "unmapped": 4,
                "secondary": 0,
                "supplementary": 0,
                "qcfail": 0,
                "human": 0,
                "viruses": 6,
            })
            read_pair_count = defaultdict(int, {
                "total": 5,
                "mapped": 3,
                "unmapped": 2,
                "primary": 2,
                "discordant": 1,
                "human": 0,
                "viruses": 2,
            })
            virus_count = defaultdict(lambda: defaultdict(int))
            fallback_virus_count = defaultdict(lambda: defaultdict(int))

            virus_count["Shared Name"]["read"] = 4
            virus_count["Shared Name"]["pair"] = 2
            fallback_virus_count["Shared Name"]["read"] = 1
            fallback_virus_count["Shared Name"]["pair"] = 1

            write_read_count_report(
                sample_id,
                read_count,
                read_pair_count,
                virus_count,
                fallback_virus_count,
            )

            report = Path(f"{sample_id}.read_count.txt").read_text(encoding="utf-8")
            self.assertIn("Shared Name\treads\t4", report)
            self.assertIn("#Fallback to organism\tCount Method\tCount Number", report)

            main_section, fallback_section = report.split("#Fallback to organism\tCount Method\tCount Number\n")
            self.assertIn("Shared Name\tpairs\t2", main_section)
            self.assertNotIn("Shared Name\treads\t1", main_section)
            self.assertIn("Shared Name\treads\t1", fallback_section)
            self.assertIn("Shared Name\tpairs\t1", fallback_section)

    def test_taxonomy_grouping_requires_taxonomy_map_argument(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            contig_path = Path(tmpdir) / "contig_info.json"
            contig_path.write_text(
                json.dumps({"ACC1.1": {"organism": "Alpha virus"}}),
                encoding="utf-8",
            )
            argv = [
                "taxonomic_classifier.py",
                "--bam",
                "sample.bam",
                "--combined_genome",
                "combined",
                "--contig_info",
                str(contig_path),
                "--virus_group_by",
                "ncbitaxon",
            ]

            with patch("sys.argv", argv):
                with self.assertRaises(SystemExit) as exc:
                    main()

            self.assertEqual(
                str(exc.exception),
                "--taxonomy_map is required when --virus_group_by is ncbitaxon or ncbitaxonname",
            )


if __name__ == "__main__":
    unittest.main()
