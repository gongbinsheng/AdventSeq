# AdVentSeq

AdVentSeq is an HPC pipeline builder and viral-quantification toolkit for sequencing data. The
`Pipeline` class is provided as an importable Python package, while the analysis tools
(`ViraQuant.py`, `TaxonomyClassifier.py`, and the helpers) are standalone scripts.

## Package data flow

The diagram below shows each script in the package, its data files (inputs and outputs), and which
files flow between scripts. Rectangles are scripts/modules, cylinders are data files, and the
rounded node is the external NCBI service. Dashed arrows denote a code dependency (shared module) or
orchestration (generated job scripts running a tool) rather than a data file.

```mermaid
flowchart TD
    %% ---- External / shared ----
    NCBI(["NCBI E-utilities<br/>(nuccore + taxonomy)"])
    UTILS["contig_info_utils.py<br/>(shared loaders)"]

    %% ===== Reference data preparation =====
    subgraph REF["Reference preparation (RVDB &rarr; taxonomy map)"]
        RVDB[("RVDB SQLite<br/>U-RVDB&lt;release&gt;.sqlite.db[.gz]")]
        CONVERT["convert_rvdb_to_json.py"]
        CONTIG[("contig_info.json.gz")]
        BUILDMAP["build_ncbi_taxonomy_map.py"]
        PREVMAP[("previous<br/>ncbi_taxonomy_map.json.gz")]
        TAXMAP[("ncbi_taxonomy_map.json.gz")]

        RVDB --> CONVERT --> CONTIG
        CONTIG --> BUILDMAP --> TAXMAP
        PREVMAP -. reuse .-> BUILDMAP
    end
    NCBI --> BUILDMAP

    %% ===== Pipeline / alignment =====
    subgraph PIPE["Per-sample pipeline"]
        FASTQ[("FASTQ reads")]
        SETTINGS[("pipeline_settings.yml")]
        PIPELINE["AdVentSeq_Pipeline.py<br/>(Pipeline class)"]
        JOBS[("generated .sh job scripts")]
        SORTBAM[("sorted BAM (+ .bai)")]
        NAMEBAM[("name-sorted BAM")]

        FASTQ --> PIPELINE
        SETTINGS --> PIPELINE
        PIPELINE --> JOBS
        JOBS --> SORTBAM
        JOBS --> NAMEBAM
    end

    %% ===== Analysis tools =====
    subgraph ANALYSIS["Analysis & reporting"]
        TAXCLASS["TaxonomyClassifier.py"]
        SPLIT[("split BAMs + unmapped FASTQ<br/>+ &lt;sample&gt;.read_count.txt")]
        VIRAQUANT["ViraQuant.py"]
        VQOUT[("ViraQuant TSV + companion .yml")]
        SCAN[("ViraQuant_scan TSV[.gz]")]
        CONVSCAN["convert_ViraQuant_scan_to_ncbitaxon.py"]
        NCBITAX[("ncbitaxon-grouped TSV[.gz]")]
        KRAKEN[("Kraken2 inspect report")]
        ADDRVDB["add_rvdb_columns.py"]
        RVDBTSV[("inspect report + in_rvdb / rvdb_accessions")]

        NAMEBAM --> TAXCLASS --> SPLIT
        SORTBAM --> VIRAQUANT
        VIRAQUANT --> VQOUT
        VIRAQUANT --> SCAN --> CONVSCAN --> NCBITAX
        KRAKEN --> ADDRVDB --> RVDBTSV
    end

    %% ---- cross-stage data files ----
    CONTIG --> TAXCLASS
    TAXMAP --> TAXCLASS
    TAXMAP --> CONVSCAN
    TAXMAP --> ADDRVDB

    %% ---- orchestration / dependencies ----
    JOBS -. runs .-> VIRAQUANT
    JOBS -. runs .-> TAXCLASS
    UTILS -. imported by .-> CONVERT
    UTILS -. imported by .-> BUILDMAP
    UTILS -. imported by .-> TAXCLASS

    classDef script fill:#dbeafe,stroke:#1e40af,color:#1e3a8a;
    classDef data fill:#dcfce7,stroke:#166534,color:#14532d;
    classDef ext fill:#fef9c3,stroke:#854d0e,color:#713f12;
    class CONVERT,BUILDMAP,PIPELINE,TAXCLASS,VIRAQUANT,CONVSCAN,ADDRVDB,UTILS script;
    class RVDB,CONTIG,PREVMAP,TAXMAP,FASTQ,SETTINGS,JOBS,SORTBAM,NAMEBAM,SPLIT,VQOUT,SCAN,NCBITAX,KRAKEN,RVDBTSV data;
    class NCBI ext;
```

## Installation

AdVentSeq targets Python 3.12. Create a dedicated conda environment and install the package in
editable mode:

```bash
conda create -n AdVentSeq python=3.12 -y
conda activate AdVentSeq
pip install -e .
```

This installs the importable `AdVentSeq` package (the `Pipeline` class) along with its runtime
dependencies (`pysam`, `tqdm`, `pyyaml`). It also registers the analysis tools as console commands
(`ViraQuant`, `TaxonomyClassifier`, etc.) on your `PATH`, so they can be run from any directory
while the environment is active.

## Basic usage

Import the package:

```python
import AdVentSeq as AVS
from AdVentSeq import Pipeline, list_files_with_extensions, is_valid_DNA_sequence

# Configure class-level settings, then build a per-sample job script.
Pipeline.WD = "/path/to/work"
Pipeline.Data_folder = "/path/to/data"

p = Pipeline(
    sample_id="Sample01",
    file_list=["/path/to/data/Sample01_L001_R1.fastq.gz",
               "/path/to/data/Sample01_L001_R2.fastq.gz"],
    log_file="/path/to/logs/Sample01.log",
    ref_genome="hg38",
    ref_type="genome",
    feature="exon",
    library_source="TRANSCRIPTOMIC",
    library_layout="paired",
    platform="ILLUMINA",
    filetype="fastq",
)
p.set_qsub_parameters()
p.set_env_variables()
# ... call step methods, then serialize p.batch to a .sh script.
```

By default the `Pipeline` loads the bundled `pipeline_settings.yml.example` config. To use your own
environment/reference settings, call `Pipeline.load_pipeline_config("/path/to/your_settings.yml")`
before creating `Pipeline` objects. See [AdVentSeq/AdVentSeq_Pipeline.md](AdVentSeq/AdVentSeq_Pipeline.md)
for the full API.

The analysis tools are installed as console commands and can be run from any directory once the
environment is active:

```bash
conda activate AdVentSeq
ViraQuant --help
TaxonomyClassifier --help
```

## Documentation

### Package (importable)

- `AdVentSeq/AdVentSeq_Pipeline.py` (the `Pipeline` class) — see
  [AdVentSeq/AdVentSeq_Pipeline.md](AdVentSeq/AdVentSeq_Pipeline.md)

### Analysis tools (console commands)

- `ViraQuant` (`AdVentSeq/ViraQuant.py`) — see [AdVentSeq/ViraQuant.md](AdVentSeq/ViraQuant.md)
- `TaxonomyClassifier` (`AdVentSeq/TaxonomyClassifier.py`) — see [AdVentSeq/TaxonomyClassifier.md](AdVentSeq/TaxonomyClassifier.md)
- `build-ncbi-taxonomy-map` (`AdVentSeq/build_ncbi_taxonomy_map.py`) — see [AdVentSeq/build_ncbi_taxonomy_map.md](AdVentSeq/build_ncbi_taxonomy_map.md)
- `convert-ViraQuant-scan-to-ncbitaxon` (`AdVentSeq/convert_ViraQuant_scan_to_ncbitaxon.py`) — see
  [AdVentSeq/convert_ViraQuant_scan_to_ncbitaxon.md](AdVentSeq/convert_ViraQuant_scan_to_ncbitaxon.md)
- `convert-rvdb-to-json` (`AdVentSeq/convert_rvdb_to_json.py`)
- `add-rvdb-columns` (`AdVentSeq/add_rvdb_columns.py`)
