# Installing AdventSeq

AdventSeq has two parts, and they are installed differently:

| Part | What it is | Install with |
| --- | --- | --- |
| The `AdventSeq` Python package | the `Pipeline` job-script builder and the analysis tools (`ViraQuant`, `TaxonomyClassifier`, `collect-results`, …) | conda **or** uv |
| The external tools | BWA, samtools, fastp, FastQC/MultiQC, Kraken2, featureCounts, GATK4 | conda only (Bioconda), one env per tool |

The generated pipeline jobs run each step inside its own conda env, including
the `AdventSeq` steps, which run in a conda env that must be named exactly
**`AdventSeq`**. So to run the full example pipeline you need conda
([Option A](#option-a-conda-recommended)). If you only want the Python package
and the analysis tools, for example to run `ViraQuant` or `collect-results` on
existing BAM files, uv ([Option B](#option-b-uv)) is enough.

## Supported systems

| System | Status |
| --- | --- |
| Linux (x86_64, aarch64) | Supported. |
| macOS (Intel and Apple Silicon) | Supported. |
| Windows 10/11 | Use **WSL2** (Windows Subsystem for Linux) and follow the Linux steps inside it. Bioconda and `pysam` do not provide native Windows packages. |

Requirements: Python ≥ 3.12 (installed for you by both options), `git`, and
internet access. The example pipeline also needs about 16 GB of RAM and
roughly 50 GB of free disk space for the reference data and indexes, plus room
for the results (see [README.md](README.md#running-the-example)).

### Windows: set up WSL2 first

In an Administrator PowerShell:

```powershell
wsl --install -d Ubuntu
```

Restart, open the **Ubuntu** app, and run every command below inside that
Ubuntu terminal. Keep the repository and data inside the Linux file system
(e.g. `~/AdventSeq`) rather than under `/mnt/c/`, which is much slower.

## Get the code

```bash
git clone https://github.com/gongbinsheng/AdventSeq.git
cd AdventSeq
```

## Option A: conda (recommended)

### 1. Install conda

If `conda --version` already works, skip this step. Otherwise install
[Miniforge](https://github.com/conda-forge/miniforge) (a minimal conda that
uses the conda-forge channel by default):

```bash
# Linux, macOS, and Windows (inside WSL2)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh" -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash    # use "init zsh" if your shell is zsh (macOS default)
```

Close and reopen the terminal so `conda` is on your `PATH`.

### 2. Create the `AdventSeq` env and install the package

The env **must** be named `AdventSeq` (the generated job scripts run
`conda activate AdventSeq`):

```bash
conda create -n AdventSeq -c conda-forge python=3.12 -y
conda activate AdventSeq
pip install -e .
```

`pip install -e .` installs the package in editable mode with its Python
dependencies (`pysam`, `pandas`, `openpyxl`, `pyyaml`, `tqdm`) and puts the
console commands on your `PATH` while the env is active.

### 3. Create the tool envs

Each pipeline step runs in its own env. Let AdventSeq write the install
commands for the envs it needs. This uses the example's env list,
[`examples/conda_tools.yml`](examples/conda_tools.yml):

```bash
conda activate AdventSeq
python -c "
from AdventSeq import Pipeline
Pipeline.set_conda_tools_config('examples/conda_tools.yml')
Pipeline.load_pipeline_config('examples/pipeline_settings.yml')
"
```

The first run reports the missing envs, writes
`examples/install_conda_envs.sh`, and exits. Review the script, then run it:

```bash
bash examples/install_conda_envs.sh
```

It creates the envs `MultiQC`, `fastp`, `BWA`, `Kraken2`, `samtools`,
`Subread` and `gatk` from the conda-forge and Bioconda channels. Run the
`python -c ...` command again. When every env exists it records the tool
versions at the end of `examples/pipeline_settings.yml` and finishes without
error.

> **Apple Silicon (M1/M2/M3/M4) Macs:** if conda reports that a Bioconda
> package is not available for `osx-arm64`, create that env with Intel
> packages (run through Rosetta 2) instead, e.g.
> `CONDA_SUBDIR=osx-64 conda create -y -n BWA -c conda-forge -c bioconda bwa`.

## Option B: uv

[uv](https://docs.astral.sh/uv/) installs the Python package and the analysis
tools only. It does not provide the external tools.

### 1. Install uv

```bash
# Linux, macOS, and Windows (inside WSL2)
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Reopen the terminal (or run `source $HOME/.local/bin/env`) so `uv` is on your `PATH`.

### 2. Install the package

From the repository root:

```bash
uv sync                      # creates .venv/ with Python 3.12 and all dependencies
source .venv/bin/activate    # puts the console commands on PATH
```

Instead of activating the venv, you can also prefix commands with `uv run`,
e.g. `uv run ViraQuant --help`.

### Running the full pipeline after a uv install

The pipeline jobs still need conda: install conda and the tool envs
([Option A](#option-a-conda-recommended), steps 1 and 3), and create the
`AdventSeq` env the jobs activate. uv can install the package into it:

```bash
conda create -n AdventSeq -c conda-forge python=3.12 -y
conda activate AdventSeq
uv pip install -e .          # installs into the active conda env
```

## Check the installation

With the env active (`conda activate AdventSeq` or `source .venv/bin/activate`):

```bash
python -c "import AdventSeq; print(AdventSeq.__version__)"
ViraQuant --help
TaxonomyClassifier --help
collect-results --help
build-ncbi-taxonomy-map --help
convert-ViraQuant-scan-to-ncbitaxon --help
convert-rvdb-to-json --help
add-rvdb-columns --help
```

Each command should print its version or help text without errors. Next, see
[Running the example](README.md#running-the-example) in the README.

## Uninstall

```bash
conda env remove -n AdventSeq          # Option A
rm -rf .venv                           # Option B
# tool envs created by install_conda_envs.sh:
for env in MultiQC fastp BWA Kraken2 samtools Subread gatk; do conda env remove -y -n "$env"; done
```
