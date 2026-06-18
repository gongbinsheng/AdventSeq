# AdVentSeq

AdVentSeq is an HPC pipeline builder and viral-quantification toolkit for sequencing data. The
`Pipeline` class is provided as an importable Python package, while the analysis tools
(`ViraQuant.py`, `taxonomy_classifier.py`, and the helpers) are standalone scripts.

## Installation

AdVentSeq targets Python 3.12. Create a dedicated conda environment and install the package in
editable mode:

```bash
conda create -n AdVentSeq python=3.12 -y
conda activate AdVentSeq
pip install -e .
```

This installs the importable `AdVentSeq` package (the `Pipeline` class) along with its runtime
dependencies (`pysam`, `tqdm`, `pyyaml`), so the same environment can also run the standalone
scripts.

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

The standalone scripts are run directly within the activated environment, e.g.:

```bash
python ViraQuant.py --help
python taxonomy_classifier.py --help
```

## Documentation

### Package (importable)

- `AdVentSeq/AdVentSeq_Pipeline.py` (the `Pipeline` class) — see
  [AdVentSeq/AdVentSeq_Pipeline.md](AdVentSeq/AdVentSeq_Pipeline.md)

### Standalone scripts

- `ViraQuant.py` — see [ViraQuant.md](ViraQuant.md)
- `taxonomy_classifier.py` — see [taxonomy_classifier.md](taxonomy_classifier.md)
- `build_ncbi_taxonomy_map.py` — see [build_ncbi_taxonomy_map.md](build_ncbi_taxonomy_map.md)
- `convert_viraquant_scan_to_ncbitaxon.py` — see
  [convert_viraquant_scan_to_ncbitaxon.md](convert_viraquant_scan_to_ncbitaxon.md)
