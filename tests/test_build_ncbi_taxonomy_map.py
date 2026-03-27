import contextlib
import io
import unittest
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from build_ncbi_taxonomy_map import build_taxonomy_map, parse_args


class FakeResponse:
    def __init__(self, payload: str):
        self.payload = payload.encode("utf-8")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read(self):
        return self.payload


class BuildNcbiTaxonomyMapTests(unittest.TestCase):
    def test_parse_args_defaults_out_path(self):
        args = parse_args([
            "--contig_info",
            "contig_info.json.gz",
            "--email",
            "Binsheng.Gong@fda.hhs.gov",
        ])
        self.assertEqual(args.out, Path("cbi_taxonomy_map.json.gz"))

    def test_email_argument_is_required(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parse_args(["--contig_info", "contig_info.json.gz"])

    def test_build_taxonomy_map_uses_email_and_preserves_unresolved_accessions(self):
        requested_urls = []
        contig_info = {
            "ACC1.1": {"organism": "Alpha virus"},
            "ACC2.1": {"organism": "Unknown virus"},
        }

        def fake_opener(url):
            requested_urls.append(url)
            parsed = urlparse(url)
            query = parse_qs(parsed.query)

            if parsed.path.endswith("/esummary.fcgi"):
                self.assertEqual(query["email"], ["Binsheng.Gong@fda.hhs.gov"])
                return FakeResponse(
                    """
                    <eSummaryResult>
                      <DocSum>
                        <Id>1</Id>
                        <Item Name="AccessionVersion" Type="String">ACC1.1</Item>
                        <Item Name="TaxId" Type="Integer">111</Item>
                      </DocSum>
                      <DocSum>
                        <Id>2</Id>
                        <Item Name="AccessionVersion" Type="String">ACC2.1</Item>
                      </DocSum>
                    </eSummaryResult>
                    """
                )

            if parsed.path.endswith("/efetch.fcgi"):
                self.assertEqual(query["email"], ["Binsheng.Gong@fda.hhs.gov"])
                return FakeResponse(
                    """
                    <TaxaSet>
                      <Taxon>
                        <TaxId>111</TaxId>
                        <ScientificName>Taxon Alpha</ScientificName>
                      </Taxon>
                    </TaxaSet>
                    """
                )

            raise AssertionError(f"Unexpected URL: {url}")

        taxonomy_map = build_taxonomy_map(
            contig_info,
            email="Binsheng.Gong@fda.hhs.gov",
            batch_size=50,
            opener=fake_opener,
        )

        self.assertEqual(
            taxonomy_map["ACC1.1"],
            {
                "ncbitaxon": "111",
                "ncbitaxonname": "Taxon Alpha",
                "organism": "Alpha virus",
            },
        )
        self.assertEqual(
            taxonomy_map["ACC2.1"],
            {
                "ncbitaxon": None,
                "ncbitaxonname": None,
                "organism": "Unknown virus",
            },
        )
        self.assertEqual(len(requested_urls), 2)


if __name__ == "__main__":
    unittest.main()
