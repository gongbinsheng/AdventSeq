#!/usr/bin/env python
"""Collect AdVentSeq per-sample pipeline outputs into consolidated tables.

The AdVentSeq pipeline writes each sample's results into a per-sample working
directory (``Pipeline.WD/$SID``).  This tool walks the pipeline output directory,
locates the relevant files for every sample listed in the SRA metadata sheet, and
assembles them into a set of tables:

* ``sample_info``                  - sample ID, title, fastp read counts (QC)
* host featureCounts matrices      - by reads and by read-pairs
* ``kraken2``                      - Kraken2 viral report matrix
* scan TaxonomyClassifier matrices - per-virus direct read/pair counts (RVDB scan)
* ``viraquant_scan_coverage``      - ViraQuant scan-mode coverage (long/tidy)
* targeted TaxonomyClassifier      - per-virus direct counts for each targeted list
* targeted ViraQuant coverage      - long table + a self-contained interactive HTML

Static tables are written both as individual ``.tsv`` files and bundled into a
single multi-sheet ``AdVentSeq_results.xlsx`` workbook.

File discovery is robust to two layouts:
  * per-sample folders  -- ``<results_dir>/<library_ID>/<library_ID>.*``
  * flattened per-type  -- ``<results_dir>/<type>/<library_ID>.*``
Files are matched by basename prefix ``<library_ID>.`` (the trailing dot prevents
``1_S1`` matching ``10_S10``).

The viral reference names used for the scan step and the targeted virus lists are
*not* hard-coded; they are supplied via ``--scan-ref`` and ``--targeted-refs`` and
matched against the alignment ``ref_genome`` embedded in the file names.
"""

import os
import re
import sys
import gzip
import json
import html
from argparse import ArgumentParser
from collections import OrderedDict, defaultdict
from pathlib import Path

import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # tqdm is a runtime dependency, but degrade gracefully.
    tqdm = None


# Curated ViraQuant measurement columns offered first in the interactive dropdown;
# any remaining numeric columns are appended after these.
VIRAQUANT_PRIMARY_FIELDS = [
    "mean_depth_ge_1",
    "mapped_reads",
    "mapped_per_bp",
    "mean_depth_all",
    "median_depth_ge_1",
    "frac_ge_1",
]
# Non-measurement identity columns in a ViraQuant TSV. "length" is a contig length,
# not a measurement, so it is excluded from the interactive measurement list too.
VIRAQUANT_ID_COLS = {"level", "virus", "seq_id", "n_contigs", "best_contig_id", "length"}

# Default per-measurement thresholds for the interactive table: a cell whose RAW value
# is below its measurement's threshold is rendered/exported empty. The rule is 1 for
# integer-valued measurements and 0.1 for float-valued ones; edit freely. Measurements
# not listed here fall back to that rule based on their observed value type. pass_*
# booleans are rendered as colour blocks and are not thresholded.
SUGGESTED_THRESHOLDS = {
    "mapped_reads": 1,
    "depth_at_50pct": 1,
    "depth_at_60pct": 1,
    "depth_at_70pct": 1,
    "depth_at_80pct": 1,
    "depth_at_90pct": 1,
    "depth_at_100pct": 1,
    "median_depth_ge_1": 1,
    "best_contig_depth_at_90pct": 1,
    "mean_depth_all": 0.1,
    "mean_depth_ge_1": 0.1,
    "frac_ge_1": 0.1,
    "breadth_cov_gt0": 0.1,
    "mapped_per_bp": 0.1,
    "best_contig_mapped_per_bp": 0.1,
    "best_contig_frac_ge_1": 0.1,
}

# Normalized-mode thresholds (value per million filtered reads). Normalized magnitudes
# are dataset-dependent, so these default to 0 (no filtering) — tune in the UI or the
# yml sidecar. Keyed by measurement name, same as SUGGESTED_THRESHOLDS.
SUGGESTED_THRESHOLDS_NORM = {}

WORKBOOK_NAME = "AdVentSeq_results.xlsx"
INTERACTIVE_HTML_NAME = "viraquant_targeted_interactive.html"

# Large/long tables: written gzip-compressed and kept out of the Excel workbook.
GZIP_TSV_TABLES = {"viraquant_scan_coverage"}
EXCEL_EXCLUDE_TABLES = {"viraquant_scan_coverage"}


def eprint(*args):
    msg = " ".join(str(a) for a in args)
    # Route through tqdm.write so messages don't corrupt an active progress bar.
    if tqdm is not None:
        tqdm.write(msg, file=sys.stderr)
    else:
        sys.stderr.write(msg + "\n")


def progress(items, desc):
    """Iterate samples with a TTY progress bar (or a plain start line otherwise)."""
    seq = list(items)
    if tqdm is not None and sys.stderr.isatty():
        return tqdm(seq, desc=desc, unit="sample", file=sys.stderr, leave=False)
    eprint("  %s (%d samples) ..." % (desc, len(seq)))
    return seq


def open_text(path):
    """Open a plain or gzip-compressed text file for reading."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


# ---------------------------------------------------------------------------
# Metadata + file discovery
# ---------------------------------------------------------------------------

def read_sample_metadata(sra_path, sheet):
    """Read library_ID + title (in sheet order) from the SRA metadata workbook."""
    df = pd.read_excel(sra_path, sheet_name=sheet, header=0)
    df = df.dropna(how="all")
    if "library_ID" not in df.columns:
        sys.exit('The "library_ID" column is not in sheet "%s".' % sheet)
    title_col = "title" if "title" in df.columns else None
    samples = OrderedDict()
    for _, row in df.iterrows():
        lib = row["library_ID"]
        if pd.isna(lib):
            continue
        lib = str(lib).strip()
        title = ""
        if title_col is not None and pd.notna(row[title_col]):
            title = str(row[title_col]).strip()
        samples[lib] = title
    if not samples:
        sys.exit("No samples found in the metadata sheet.")
    return samples


def index_result_files(results_dir):
    """Walk results_dir once and return a list of (basename, full_path)."""
    files = []
    for root, _dirs, names in os.walk(results_dir):
        for name in names:
            files.append((name, os.path.join(root, name)))
    return files


def files_for_sample(file_index, library_id):
    """Return files whose basename starts with '<library_id>.'."""
    prefix = library_id + "."
    return [(n, p) for (n, p) in file_index if n.startswith(prefix)]


# ---------------------------------------------------------------------------
# Mapper / reference detection from file names
# ---------------------------------------------------------------------------

# host-stage files:  <lib>.<host_mapper>.hg38.<...>
_HOST_RE = re.compile(r"^[^.]+\.(?P<host>[^.]+)\.hg38\.")
# virus-stage alignment outputs: <lib>.<host_mapper>.hg38.unmapped2host.<virus_mapper>.<ref>.<...>
_VIRUS_RE = re.compile(
    r"^[^.]+\.(?P<host>[^.]+)\.hg38\.unmapped2host\.(?P<virus>[^.]+)\.(?P<ref>[^.]+)\."
)
# Only these suffixes are genuine virus-alignment outputs. Other files share the
# "...unmapped2host.<token>.<token>..." shape but are NOT virus alignments, e.g.
# remove_read_pairs_mapped_to_host FASTQs (.unmapped2host.R1/R2/single/orphan.fastq.gz)
# and Kraken2 (.unmapped2host.Kraken2_viral.report); they must not contribute a
# spurious virus_mapper during auto-detection.
_VIRUS_OUTPUT_SUFFIXES = (
    ".unsorted.read_count.txt",
    ".ViraQuant.tsv", ".ViraQuant.tsv.gz",
    ".ViraQuant_scan.tsv", ".ViraQuant_scan.tsv.gz",
)


def _select_token(name, present, requested, what):
    """Auto-select a mapper when exactly one is present, else honour the request."""
    if requested:
        if requested not in present:
            eprint("Warning: requested %s '%s' not found (available: %s)."
                   % (what, requested, ", ".join(sorted(present)) or "none"))
        return requested
    if len(present) == 1:
        return next(iter(present))
    if len(present) == 0:
        return None
    sys.exit("Multiple %s values found: %s. Re-run with --%s to choose one."
             % (what, ", ".join(sorted(present)), what.replace("_", "-")))


def detect_mappers(file_index, samples, host_mapper, virus_mapper):
    """Scan all sample files to determine the host/virus mapper to use."""
    hosts, viruses = set(), set()
    for lib in samples:
        for name, _ in files_for_sample(file_index, lib):
            # Host mapper: the token before ".hg38." is the genuine host aligner in
            # every stage's file names.
            m = _HOST_RE.match(name)
            if m:
                hosts.add(m.group("host"))
            # Virus mapper: only trust genuine virus-alignment outputs.
            if name.endswith(_VIRUS_OUTPUT_SUFFIXES):
                mv = _VIRUS_RE.match(name)
                if mv:
                    viruses.add(mv.group("virus"))
    host = _select_token(None, hosts, host_mapper, "host_mapper")
    virus = _select_token(None, viruses, virus_mapper, "virus_mapper")
    return host, virus


# ---------------------------------------------------------------------------
# Per-file parsers
# ---------------------------------------------------------------------------

def parse_fastp(path):
    """Return total reads before/after filtering + passed-filter reads."""
    with open_text(path) as fh:
        data = json.load(fh)
    summary = data.get("summary", {})
    before = summary.get("before_filtering", {}).get("total_reads")
    after = summary.get("after_filtering", {}).get("total_reads")
    passed = data.get("filtering_result", {}).get("passed_filter_reads")
    return before, after, passed


def parse_featurecounts(path):
    """Return a (Geneid -> count) Series from a featureCounts TSV(.gz).

    The first line is a '# Program:' comment; the header follows with the count
    in the final column (the BAM file).
    """
    with open_text(path) as fh:
        first = fh.readline()
        if not first.startswith("#"):
            # No comment line; rewind by re-opening.
            fh.close()
            fh = open_text(path)
        df = pd.read_csv(fh, sep="\t")
    count_col = df.columns[-1]
    series = df.set_index("Geneid")[count_col]
    series.index = series.index.astype(str)
    return series


def parse_kraken2(path):
    """Parse a Kraken2 report into a DataFrame.

    Columns: pct, clade_reads, taxon_reads, rank_code, ncbi_taxid, name. The name
    is leading-space indented in the report to convey the rank depth; that
    indentation is preserved (only the trailing newline is removed).
    """
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pct, clade_reads, taxon_reads, rank_code, taxid = parts[:5]
            # Keep the leading whitespace: it encodes the taxon's rank depth.
            name = parts[5].rstrip()
            rows.append({
                "ncbi_taxid": taxid.strip(),
                "rank_code": rank_code.strip(),
                "name": name,
                "clade_reads": int(clade_reads),
                "taxon_reads": int(taxon_reads),
            })
    return pd.DataFrame(rows)


def parse_read_count(path):
    """Parse a TaxonomyClassifier read_count.txt.

    Returns (summary, reads, pairs) where summary maps the '#'-prefixed labels to
    their count, and reads/pairs map each per-virus key to its count.
    """
    summary = OrderedDict()
    reads = OrderedDict()
    pairs = OrderedDict()
    fallback = False
    with open_text(path) as fh:
        next(fh, None)  # header line ("<sample>\tCount Method\tCount Number")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if line.startswith("#Fallback"):
                fallback = True
                continue
            if len(parts) < 3:
                continue
            key, method, value = parts[0], parts[1], parts[2]
            if key.startswith("#"):
                summary[key.lstrip("#")] = value
                continue
            label = ("%s (fallback)" % key) if fallback else key
            if method == "reads":
                reads[label] = int(value)
            elif method == "pairs":
                pairs[label] = int(value)
    return summary, reads, pairs


def parse_viraquant(path):
    """Parse a ViraQuant TSV(.gz) into a DataFrame (all columns preserved)."""
    with open_text(path) as fh:
        df = pd.read_csv(fh, sep="\t", dtype={"seq_id": str})
    return df


# ---------------------------------------------------------------------------
# File locators (apply mapper / ref filters)
# ---------------------------------------------------------------------------

def find_host_file(sample_files, host_mapper, suffix):
    """featureCounts file: <lib>.<host>.hg38.featureCounts.*.<suffix>."""
    cands = []
    for name, path in sample_files:
        if not name.endswith(suffix):
            continue
        if host_mapper and (".%s.hg38." % host_mapper) not in name:
            continue
        cands.append((name, path))
    return cands[0][1] if cands else None


def find_kraken_file(sample_files, host_mapper):
    cands = []
    for name, path in sample_files:
        if not name.endswith(".Kraken2_viral.report"):
            continue
        if host_mapper and (".%s.hg38." % host_mapper) not in name:
            continue
        cands.append(path)
    return cands[0] if cands else None


def find_virus_file(sample_files, host_mapper, virus_mapper, ref, suffix,
                    exclude=None):
    """Locate a read_count/ViraQuant file for a given virus mapper + reference."""
    cands = []
    for name, path in sample_files:
        if not (name.endswith(suffix) or name.endswith(suffix + ".gz")):
            continue
        if exclude and exclude in name:
            continue
        if host_mapper and (".%s.hg38." % host_mapper) not in name:
            continue
        if virus_mapper and ref:
            if (".%s.%s." % (virus_mapper, ref)) not in name:
                continue
        elif ref and (".%s." % ref) not in name:
            continue
        cands.append(path)
    return cands[0] if cands else None


# ---------------------------------------------------------------------------
# Table builders
# ---------------------------------------------------------------------------

def build_sample_info(samples, file_index):
    rows = []
    for lib in progress(samples, "fastp QC / sample info"):
        title = samples[lib]
        sample_files = files_for_sample(file_index, lib)
        before = after = passed = None
        fastp_path = None
        for name, path in sample_files:
            if name.endswith(".fastp.json"):
                fastp_path = path
                break
        if fastp_path:
            try:
                before, after, passed = parse_fastp(fastp_path)
            except Exception as exc:  # noqa: BLE001
                eprint("Warning: failed to parse fastp for %s: %s" % (lib, exc))
        else:
            eprint("Warning: no fastp.json for sample %s" % lib)
        frac = (after / before) if (before and after is not None) else None
        rows.append({
            "library_ID": lib,
            "title": title,
            "fastp_before_total_reads": before,
            "fastp_after_total_reads": after,
            "passed_filter_reads": passed,
            "after_vs_before_fraction": frac,
        })
    return pd.DataFrame(rows)


def build_featurecounts_matrix(samples, file_index, host_mapper, count_kind):
    """count_kind is 'count_reads' or 'count_pairs'."""
    series_by_sample = OrderedDict()
    for lib in progress(samples, "host featureCounts (%s)" % count_kind):
        sample_files = files_for_sample(file_index, lib)
        path = find_host_file(sample_files, host_mapper,
                              "%s.tsv.gz" % count_kind)
        if path is None:
            path = find_host_file(sample_files, host_mapper, "%s.tsv" % count_kind)
        if path is None:
            eprint("Warning: no featureCounts %s for sample %s" % (count_kind, lib))
            continue
        try:
            series_by_sample[lib] = parse_featurecounts(path)
        except Exception as exc:  # noqa: BLE001
            eprint("Warning: failed featureCounts %s for %s: %s"
                   % (count_kind, lib, exc))
    if not series_by_sample:
        return None
    matrix = pd.DataFrame(series_by_sample)
    matrix.index.name = "Geneid"
    return matrix.reset_index()


def build_kraken2_matrix(samples, file_index, host_mapper):
    per_sample = OrderedDict()
    identity = OrderedDict()  # (taxid, rank) -> name
    for lib in progress(samples, "Kraken2"):
        sample_files = files_for_sample(file_index, lib)
        path = find_kraken_file(sample_files, host_mapper)
        if path is None:
            eprint("Warning: no Kraken2 report for sample %s" % lib)
            continue
        df = parse_kraken2(path)
        if df.empty:
            continue
        col = OrderedDict()
        for _, r in df.iterrows():
            key = (r["ncbi_taxid"], r["rank_code"])
            identity[key] = r["name"]
            col[key] = r["clade_reads"]
        per_sample[lib] = col
    if not per_sample:
        return None
    keys = list(identity.keys())
    rows = []
    for key in keys:
        row = {
            "ncbi_taxid": key[0],
            "rank_code": key[1],
            "name": identity[key],
        }
        for lib, col in per_sample.items():
            row[lib] = col.get(key, 0)
        rows.append(row)
    return pd.DataFrame(rows)


def build_taxonomy_matrices(samples, file_index, host_mapper, virus_mapper, ref,
                            taxonomy_map=None):
    """Return (reads_matrix, pairs_matrix, summary_df) for a virus reference."""
    reads_by_sample = OrderedDict()
    pairs_by_sample = OrderedDict()
    summary_rows = []
    found = False
    for lib in progress(samples, "TaxonomyClassifier [%s]" % ref):
        sample_files = files_for_sample(file_index, lib)
        path = find_virus_file(sample_files, host_mapper, virus_mapper, ref,
                               ".unsorted.read_count.txt")
        if path is None:
            eprint("Warning: no %s read_count for sample %s" % (ref, lib))
            continue
        found = True
        summary, reads, pairs = parse_read_count(path)
        reads_by_sample[lib] = reads
        pairs_by_sample[lib] = pairs
        srow = {"library_ID": lib}
        srow.update(summary)
        summary_rows.append(srow)
    if not found:
        return None, None, None

    def to_matrix(by_sample):
        keys = list(OrderedDict((k, None) for d in by_sample.values() for k in d))
        rows = []
        for key in keys:
            row = {"virus": key}
            if taxonomy_map is not None:
                row["name"] = taxonomy_map.get(str(key), "")
            for lib, d in by_sample.items():
                row[lib] = d.get(key, 0)
            rows.append(row)
        return pd.DataFrame(rows)

    return to_matrix(reads_by_sample), to_matrix(pairs_by_sample), \
        pd.DataFrame(summary_rows)


def build_viraquant_long(samples, file_index, host_mapper, virus_mapper, ref,
                         scan):
    """Concatenate ViraQuant TSVs across samples into a long table."""
    suffix = ".ViraQuant_scan.tsv" if scan else ".ViraQuant.tsv"
    exclude = None if scan else "_scan"
    desc = "ViraQuant %s [%s]" % ("scan" if scan else "targeted", ref)
    frames = []
    for lib in progress(samples, desc):
        sample_files = files_for_sample(file_index, lib)
        path = find_virus_file(sample_files, host_mapper, virus_mapper, ref,
                               suffix, exclude=exclude)
        if path is None:
            eprint("Warning: no %s ViraQuant%s for sample %s"
                   % (ref, "_scan" if scan else "", lib))
            continue
        df = parse_viraquant(path)
        df.insert(0, "sample_id", lib)
        frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Taxonomy map (optional) for readable scan names
# ---------------------------------------------------------------------------

def load_taxonomy_map(path):
    if not path:
        return None
    with open_text(path) as fh:
        raw = json.load(fh)
    out = {}
    for acc, rec in raw.items():
        if isinstance(rec, dict):
            taxid = rec.get("ncbitaxon")
            name = rec.get("ncbitaxonname") or rec.get("organism")
            if taxid is not None and name:
                out.setdefault(str(taxid), name)
    return out


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _sheet_name(base, used):
    """Excel sheet names: <=31 chars, unique."""
    name = base[:31]
    if name not in used:
        used.add(name)
        return name
    i = 1
    while True:
        suffix = "_%d" % i
        cand = name[:31 - len(suffix)] + suffix
        if cand not in used:
            used.add(cand)
            return cand
        i += 1


def write_outputs(tables, out_dir):
    """tables: list of (name, DataFrame). Writes TSVs + one xlsx workbook.

    Tables in GZIP_TSV_TABLES are written gzip-compressed; tables in
    EXCEL_EXCLUDE_TABLES are kept out of the workbook (they are large/long).
    """
    os.makedirs(out_dir, exist_ok=True)
    for name, df in tables:
        if name in GZIP_TSV_TABLES:
            df.to_csv(os.path.join(out_dir, name + ".tsv.gz"), sep="\t",
                      index=False, compression="gzip")
        else:
            df.to_csv(os.path.join(out_dir, name + ".tsv"), sep="\t", index=False)
    xlsx_path = os.path.join(out_dir, WORKBOOK_NAME)
    used = set()
    n_sheets = 0
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for name, df in tables:
            if name in EXCEL_EXCLUDE_TABLES:
                continue
            df.to_excel(writer, sheet_name=_sheet_name(name, used), index=False)
            n_sheets += 1
    eprint("Wrote %d tables to %s (%d sheets in %s)"
           % (len(tables), out_dir, n_sheets, WORKBOOK_NAME))


# ---------------------------------------------------------------------------
# Interactive HTML
# ---------------------------------------------------------------------------

def build_interactive_html(targeted_long, out_dir, filtered_reads=None,
                           sample_titles=None):
    """Build a self-contained HTML viewer for targeted ViraQuant coverage.

    targeted_long: dict ref_name -> long DataFrame (with 'sample_id', 'level',
    'virus', and measurement columns).
    filtered_reads: {sample_id: fastp after-filtering total reads} for CPM-style
    normalization in the viewer.
    sample_titles: {sample_id: title} shown next to the sample id in the table.
    """
    refs_payload = OrderedDict()
    measurement_cols = []
    meas_is_int = {}  # measurement -> all observed values are integral
    for ref, df in targeted_long.items():
        if df is None or df.empty:
            continue
        # Prefer virus-level rows; for any sample lacking them, fall back to contig.
        virus_df = df[df["level"] == "virus"] if "level" in df.columns else df
        if virus_df.empty:
            virus_df = df
        else:
            have = set(virus_df["sample_id"].unique())
            missing = set(df["sample_id"].unique()) - have
            if missing and "level" in df.columns:
                extra = df[(df["sample_id"].isin(missing)) & (df["level"] == "contig")]
                virus_df = pd.concat([virus_df, extra], ignore_index=True)
        numeric_cols = [c for c in virus_df.columns
                        if c not in VIRAQUANT_ID_COLS and c != "sample_id"
                        and pd.api.types.is_numeric_dtype(virus_df[c])]
        ordered = [c for c in VIRAQUANT_PRIMARY_FIELDS if c in numeric_cols]
        ordered += [c for c in numeric_cols if c not in ordered]
        for c in ordered:
            if c not in measurement_cols:
                measurement_cols.append(c)
            meas_is_int.setdefault(c, True)
        samples = list(dict.fromkeys(virus_df["sample_id"].tolist()))
        viruses = list(dict.fromkeys(virus_df["virus"].astype(str).tolist()))
        records = []
        for _, r in virus_df.iterrows():
            rec = {"sample_id": r["sample_id"], "virus": str(r["virus"])}
            for c in ordered:
                val = r[c]
                if pd.isna(val):
                    rec[c] = None
                elif float(val).is_integer():
                    rec[c] = int(val)
                else:
                    rec[c] = float(val)
                    meas_is_int[c] = False
            records.append(rec)
        refs_payload[ref] = {
            "samples": samples,
            "viruses": viruses,
            "records": records,
        }
    if not refs_payload:
        eprint("Warning: no targeted ViraQuant data; skipping interactive HTML.")
        return None

    # Per-measurement metadata: integer-ness, default raw/normalized thresholds, flag.
    meta = OrderedDict()
    for c in measurement_cols:
        is_pass = c.startswith("pass_")
        is_int = meas_is_int.get(c, True)
        thr_raw = SUGGESTED_THRESHOLDS.get(c, 1 if is_int else 0.1)
        thr_norm = SUGGESTED_THRESHOLDS_NORM.get(c, 0)
        meta[c] = {"is_integer": is_int, "is_pass": is_pass,
                   "threshold_raw": thr_raw, "threshold_norm": thr_norm}

    payload = {
        "refs": refs_payload,
        "measurements": measurement_cols,
        "meta": meta,
        "filtered_reads": filtered_reads or {},
        "sample_titles": sample_titles or {},
        "default_measurement": (VIRAQUANT_PRIMARY_FIELDS[0]
                                if VIRAQUANT_PRIMARY_FIELDS[0] in measurement_cols
                                else measurement_cols[0]),
    }
    data_json = json.dumps(payload)
    html_doc = _HTML_TEMPLATE.replace("__DATA__", data_json)
    out_path = os.path.join(out_dir, INTERACTIVE_HTML_NAME)
    with open(out_path, "w") as fh:
        fh.write(html_doc)
    eprint("Wrote interactive table: %s" % out_path)
    return out_path


_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AdVentSeq - Targeted ViraQuant coverage</title>
<style>
  body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
         margin: 1.5rem; color: #1f2937; }
  h1 { font-size: 1.3rem; }
  .controls { margin: 1rem 0; display: flex; gap: 1.2rem; flex-wrap: wrap;
              align-items: center; }
  .controls > div { display: flex; align-items: center; }
  label { font-weight: 600; margin-right: .4rem; }
  select, input[type=number] { font-size: 1rem; padding: .25rem .4rem; }
  input[type=number] { width: 6.5rem; }
  button { font-size: .9rem; padding: .35rem .7rem; cursor: pointer;
           border: 1px solid #9ca3af; border-radius: 4px; background: #f3f4f6; }
  button:hover { background: #e5e7eb; }
  table { border-collapse: collapse; margin-top: .5rem; font-variant-numeric: tabular-nums; }
  th, td { border: 1px solid #d1d5db; padding: .3rem .55rem; text-align: right;
           white-space: nowrap; }
  /* Two sticky left columns: sample id, then sample title. */
  th.scol, td.scol { text-align: left; position: sticky; background: #f9fafb;
                     font-weight: 600; }
  th.scol1, td.scol1 { left: 0; }
  th.scol2, td.scol2 { left: var(--col1w, 7rem); color: #4b5563; font-weight: 400; }
  thead th { background: #1e3a8a; color: #fff; position: sticky; top: 0; z-index: 2; }
  thead th.scol { background: #1e3a8a; color: #fff; z-index: 3; }
  tbody td.scol { z-index: 1; }
  .wrap { overflow: auto; max-height: 75vh; max-width: 100%; }
  .meta { color: #6b7280; font-size: .85rem; }
  td.zero { color: #9ca3af; }
  /* Boolean pass_* measurements are shown as colour blocks (no number). */
  td.passcell { padding: 0; }
  td.passtrue { background: #16a34a; }
  .legend { color: #6b7280; font-size: .85rem; margin: .25rem 0 0; min-height: 1.1rem; }
  .swatch { display: inline-block; width: .9rem; height: .9rem; vertical-align: middle;
            border: 1px solid #d1d5db; margin: 0 .25rem 0 .6rem; }
  .swatch.on { background: #16a34a; }
  .status { color: #374151; font-size: .85rem; margin: .35rem 0 0; }
</style>
</head>
<body>
<h1>Targeted ViraQuant coverage</h1>
<p class="meta">Rows are samples (id + title), columns are viruses. Cells below the active
threshold are shown empty. Each measurement has separate <b>raw</b> and <b>normalized</b>
thresholds that switch with the toggle. Normalized = value per million filtered (fastp
after-filtering) reads. Threshold edits persist in this browser; <b>Download table</b>
saves the TSVs to your browser's downloads folder.</p>
<div class="controls">
  <div id="refbox" style="display:none">
    <label for="refsel">Reference</label>
    <select id="refsel"></select>
  </div>
  <div>
    <label for="measuresel">Measurement</label>
    <select id="measuresel"></select>
  </div>
  <div>
    <label><input type="checkbox" id="normtoggle" checked> Normalized</label>
  </div>
  <div>
    <label for="thrinput">Threshold</label>
    <input type="number" id="thrinput" step="any" min="0">
    <button id="applybtn" title="Apply the threshold to the table and save thresholds to .yml">Apply</button>
  </div>
  <div>
    <button id="dlbtn" title="Download raw + threshold-filtered TSV of the current table">Download table (.tsv)</button>
  </div>
</div>
<div class="legend" id="legend"></div>
<div class="status" id="status"></div>
<div class="wrap"><table id="tbl"><thead></thead><tbody></tbody></table></div>
<script>
const DATA = __DATA__;
const LS_KEY = 'adventseq_viraquant_thresholds';

const refSel = document.getElementById('refsel');
const measureSel = document.getElementById('measuresel');
const refbox = document.getElementById('refbox');
const thrInput = document.getElementById('thrinput');
const normToggle = document.getElementById('normtoggle');

const refNames = Object.keys(DATA.refs);
refNames.forEach(r => {
  const o = document.createElement('option'); o.value = r; o.textContent = r;
  refSel.appendChild(o);
});
if (refNames.length > 1) { refbox.style.display = 'flex'; }

DATA.measurements.forEach(m => {
  const o = document.createElement('option'); o.value = m; o.textContent = m;
  if (m === DATA.default_measurement) o.selected = true;
  measureSel.appendChild(o);
});

// thresholds: flat map "<measure>.raw" / "<measure>.norm". start from embedded
// defaults; init() overlays sidecar yml + localStorage.
const thresholds = {};
DATA.measurements.forEach(m => {
  const mt = DATA.meta[m];
  if (mt && !mt.is_pass) {
    thresholds[m + '.raw'] = mt.threshold_raw;
    thresholds[m + '.norm'] = mt.threshold_norm;
  }
});

function metaOf(m) { return DATA.meta[m] || {}; }
function curMode() { return normToggle.checked ? 'norm' : 'raw'; }
function thrOf(measure) { return thresholds[measure + '.' + curMode()]; }

function fmt(v) {
  if (v === null || v === undefined) return '';
  if (Number.isInteger(v)) return v.toLocaleString();
  return (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(3) : v.toFixed(4);
}

function esc(s) {
  return String(s).replace(/[&<>"]/g, c =>
    ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]));
}

// The second sticky column (title) must offset by the rendered width of the first.
function fixStickyOffsets() {
  const c1 = document.querySelector('#tbl tbody td.scol1')
          || document.querySelector('#tbl thead th.scol1');
  if (c1) {
    document.getElementById('tbl').style.setProperty(
      '--col1w', c1.getBoundingClientRect().width + 'px');
  }
}

// normalized value = raw / (filtered_reads[sample] / 1e6); null if reads missing/zero.
function normalize(raw, sample) {
  const reads = DATA.filtered_reads[sample];
  if (!reads) return null;
  return raw / (reads / 1e6);
}

function persist() {
  try { localStorage.setItem(LS_KEY, JSON.stringify(thresholds)); } catch (e) {}
}
function loadLS() {
  try { return JSON.parse(localStorage.getItem(LS_KEY) || '{}'); } catch (e) { return {}; }
}

function downloadText(name, text, mime) {
  const blob = new Blob([text], {type: mime || 'text/plain'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
}

function safeName(s) { return String(s).replace(/[^A-Za-z0-9_.-]+/g, '_'); }

// Build a sample x virus TSV. filtered=false -> raw values, all cells. filtered=true ->
// the currently displayed values (raw or normalized) with the active threshold applied.
function tableToTsv(filtered) {
  const ref = refSel.value || refNames[0];
  const measure = measureSel.value;
  const isPass = !!metaOf(measure).is_pass;
  const normalized = normToggle.checked && !isPass;
  const thr = isPass ? 1 : thrOf(measure);
  const block = DATA.refs[ref];
  const samples = block.samples, viruses = block.viruses;
  const lookup = {};
  block.records.forEach(rec => { lookup[rec.virus + '\\u0000' + rec.sample_id] = rec[measure]; });
  const lines = ['sample_id\\ttitle\\t' + viruses.join('\\t')];
  samples.forEach(s => {
    const title = DATA.sample_titles[s] || '';
    const cells = viruses.map(v => {
      const raw = lookup[v + '\\u0000' + s];
      if (raw === null || raw === undefined) return '';
      if (!filtered) return String(raw);
      if (isPass) return raw >= 1 ? '1' : '';
      const val = normalized ? normalize(raw, s) : raw;
      if (val === null) return '';
      if (thr !== undefined && thr !== null && val < thr) return '';
      return String(val);
    });
    lines.push(s + '\\t' + title + '\\t' + cells.join('\\t'));
  });
  return lines.join('\\n') + '\\n';
}

function downloadTable() {
  const ref = safeName(refSel.value || refNames[0]);
  const measure = safeName(measureSel.value);
  downloadText('viraquant_' + ref + '_' + measure + '_raw.tsv',
               tableToTsv(false), 'text/tab-separated-values');
  downloadText('viraquant_' + ref + '_' + measure + '_filtered.tsv',
               tableToTsv(true), 'text/tab-separated-values');
}

// Sync the threshold input to the current measurement + view mode (raw vs normalized).
function syncControls() {
  const measure = measureSel.value;
  const isPass = !!metaOf(measure).is_pass;
  thrInput.disabled = isPass;
  normToggle.disabled = isPass;
  thrInput.value = isPass ? '' : thrOf(measure);
}

function render() {
  const ref = refSel.value || refNames[0];
  const measure = measureSel.value;
  const isPass = !!metaOf(measure).is_pass;
  const normalized = normToggle.checked && !isPass;
  // Threshold applies to the displayed (raw or normalized) value, using the matching key.
  const thr = isPass ? null : thrOf(measure);
  const block = DATA.refs[ref];
  const samples = block.samples;
  const viruses = block.viruses;
  const lookup = {};
  block.records.forEach(rec => { lookup[rec.virus + '\\u0000' + rec.sample_id] = rec[measure]; });

  document.getElementById('legend').innerHTML = isPass
    ? 'Boolean: <span class="swatch on"></span> true (1) <span class="swatch"></span> false (0)'
    : '';
  document.getElementById('status').textContent = isPass
    ? 'Showing pass/fail blocks.'
    : ('Values: ' + (normalized ? 'normalized (per million filtered reads)' : 'raw')
       + ' | threshold (on ' + (normalized ? 'normalized' : 'raw') + '): ' + thr);

  const thead = document.querySelector('#tbl thead');
  const tbody = document.querySelector('#tbl tbody');
  // Rows are samples (id + title); columns are viruses.
  let head = '<tr><th class="scol scol1">sample</th><th class="scol scol2">title</th>';
  viruses.forEach(v => { head += '<th>' + esc(v) + '</th>'; });
  head += '</tr>';
  thead.innerHTML = head;

  let body = '';
  samples.forEach(s => {
    const title = DATA.sample_titles[s] || '';
    body += '<tr><td class="scol scol1">' + esc(s) + '</td>'
          + '<td class="scol scol2">' + esc(title) + '</td>';
    viruses.forEach(v => {
      const raw = lookup[v + '\\u0000' + s];
      if (isPass) {
        const cls = (raw === 1) ? 'passcell passtrue' : 'passcell';
        const ttl = (raw === 1) ? 'true' : (raw === 0 ? 'false' : 'n/a');
        body += '<td class="' + cls + '" title="' + ttl + '"></td>';
        return;
      }
      // Compute the displayed value, then filter on it with the active threshold.
      let shown = raw;
      if (normalized) shown = (raw === null || raw === undefined) ? null : normalize(raw, s);
      if (shown === null || shown === undefined) {  // missing, or no read count
        body += '<td class="zero"></td>';
        return;
      }
      if (thr !== undefined && thr !== null && shown < thr) {
        body += '<td class="zero"></td>';
        return;
      }
      body += '<td>' + fmt(shown) + '</td>';
    });
    body += '</tr>';
  });
  tbody.innerHTML = body;
  fixStickyOffsets();
}

refSel.addEventListener('change', render);
measureSel.addEventListener('change', () => { syncControls(); render(); });
// Toggling raw/normalized swaps in that mode's threshold for the current measurement.
normToggle.addEventListener('change', () => { syncControls(); render(); });

// Apply: set the current measurement+mode threshold from the input, persist to
// localStorage, and re-render.
function applyThreshold() {
  const measure = measureSel.value;
  if (!metaOf(measure).is_pass) {
    const v = Number(thrInput.value);
    if (!isNaN(v)) thresholds[measure + '.' + curMode()] = v;
  }
  persist();
  render();
}
document.getElementById('applybtn').addEventListener('click', applyThreshold);
thrInput.addEventListener('keydown', (e) => { if (e.key === 'Enter') applyThreshold(); });
document.getElementById('dlbtn').addEventListener('click', downloadTable);

function init() {
  // localStorage overrides the embedded defaults (persists in this browser, file:// ok).
  Object.assign(thresholds, loadLS());
  syncControls();
  render();
}
init();
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = ArgumentParser(description="Collect AdVentSeq per-sample results into tables.")
    parser.add_argument("-s", "--sra", dest="sra", required=True,
                        help="SRA metadata workbook (.xlsx).")
    parser.add_argument("--sheet", dest="sheet", default="SRA_data",
                        help="Worksheet name [default: %(default)s].")
    parser.add_argument("-r", "--results-dir", dest="results_dir", required=True,
                        help="Pipeline output directory (the '-o' dir of the run).")
    parser.add_argument("-o", "--out-dir", dest="out_dir", default="./collected_tables",
                        help="Directory for the collected tables [default: %(default)s].")
    parser.add_argument("--scan-ref", dest="scan_ref", required=True,
                        help="Reference name of the full-BAM scan / direct-count "
                             "step (matched in file names), e.g. 'RVDBv29'.")
    parser.add_argument("--targeted-refs", dest="targeted_refs", required=True,
                        help="Comma-separated targeted virus-list reference names, "
                             "e.g. '7viruses,9viruses'.")
    parser.add_argument("--host-mapper", dest="host_mapper", default=None,
                        help="Host aligner name (auto-detected if only one present).")
    parser.add_argument("--virus-mapper", dest="virus_mapper", default=None,
                        help="Virus aligner name (auto-detected if only one present).")
    parser.add_argument("--taxonomy-map", dest="taxonomy_map", default=None,
                        help="Optional ncbi_taxonomy_map.json[.gz] to add readable "
                             "names to the scan-ref taxonomy tables.")
    o = parser.parse_args()

    results_dir = Path(o.results_dir).expanduser().resolve()
    if not results_dir.is_dir():
        sys.exit("Results directory does not exist: %s" % results_dir)

    samples = read_sample_metadata(o.sra, o.sheet)
    targeted_refs = [r.strip() for r in o.targeted_refs.split(",") if r.strip()]
    eprint("Indexing files under %s ..." % results_dir)
    file_index = index_result_files(str(results_dir))
    if not file_index:
        sys.exit("No files found under %s" % results_dir)
    eprint("Indexed %d files." % len(file_index))

    host_mapper, virus_mapper = detect_mappers(
        file_index, samples, o.host_mapper, o.virus_mapper)
    eprint("Samples: %d | host_mapper=%s | virus_mapper=%s | scan_ref=%s | targeted=%s"
           % (len(samples), host_mapper, virus_mapper, o.scan_ref,
              ",".join(targeted_refs)))
    eprint("Collecting results ...")

    tax_map = load_taxonomy_map(o.taxonomy_map)
    tables = []

    # 1. sample info / fastp QC
    sample_info = build_sample_info(samples, file_index)
    tables.append(("sample_info", sample_info))
    # Per-sample filtered (fastp after-filtering) total reads for CPM normalization.
    filtered_reads = {
        str(lib): (int(reads) if pd.notna(reads) else None)
        for lib, reads in zip(sample_info["library_ID"],
                              sample_info["fastp_after_total_reads"])
    }
    # Per-sample titles shown alongside the sample id in the interactive table.
    sample_titles = {
        str(lib): ("" if pd.isna(title) else str(title))
        for lib, title in zip(sample_info["library_ID"], sample_info["title"])
    }

    # 2. host gene quantification
    for kind, label in (("count_reads", "host_featureCounts_by_reads"),
                        ("count_pairs", "host_featureCounts_by_pairs")):
        mat = build_featurecounts_matrix(samples, file_index, host_mapper, kind)
        if mat is not None:
            tables.append((label, mat))

    # 3. Kraken2
    kraken = build_kraken2_matrix(samples, file_index, host_mapper)
    if kraken is not None:
        tables.append(("kraken2", kraken))

    # 4. scan TaxonomyClassifier direct counts (RVDB scan reference)
    reads_m, pairs_m, summary = build_taxonomy_matrices(
        samples, file_index, host_mapper, virus_mapper, o.scan_ref, tax_map)
    if reads_m is not None:
        tables.append(("scan_taxonomy_reads", reads_m))
        tables.append(("scan_taxonomy_pairs", pairs_m))
        tables.append(("scan_summary", summary))

    # 5. ViraQuant scan coverage (long)
    scan_cov = build_viraquant_long(
        samples, file_index, host_mapper, virus_mapper, o.scan_ref, scan=True)
    if scan_cov is not None:
        tables.append(("viraquant_scan_coverage", scan_cov))

    # 6 + 7. per targeted reference
    targeted_long = OrderedDict()
    for ref in targeted_refs:
        reads_m, pairs_m, summary = build_taxonomy_matrices(
            samples, file_index, host_mapper, virus_mapper, ref)
        if reads_m is not None:
            tables.append(("targeted_taxonomy_reads_%s" % ref, reads_m))
            tables.append(("targeted_taxonomy_pairs_%s" % ref, pairs_m))
            tables.append(("targeted_summary_%s" % ref, summary))
        cov = build_viraquant_long(
            samples, file_index, host_mapper, virus_mapper, ref, scan=False)
        if cov is not None:
            tables.append(("viraquant_targeted_coverage_%s" % ref, cov))
            targeted_long[ref] = cov

    if len(tables) <= 1:
        eprint("Warning: only sample_info was assembled; check --results-dir, "
               "mapper names, and reference names.")

    eprint("Writing %d tables to %s ..." % (len(tables), o.out_dir))
    write_outputs(tables, o.out_dir)
    build_interactive_html(targeted_long, o.out_dir, filtered_reads=filtered_reads,
                           sample_titles=sample_titles)
    eprint("Done.")


if __name__ == "__main__":
    main()
