#!/usr/bin/env bash
# Generate the per-sample job scripts for the example study, then run them.
#
# Before running: install AdventSeq (see INSTALL.md), activate the env
#   conda activate AdventSeq
# and prepare the FASTQ files and reference data (see "Running the example"
# in README.md). Run from anywhere:
#   bash examples/run_example.sh            # run the samples on this computer
#   bash examples/run_example.sh slurm      # submit the samples to Slurm
set -euo pipefail

EXAMPLES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCHEDULER="${1:-local}"
THREADS="${THREADS:-4}"
RESULT_DIR="$EXAMPLES_DIR/Results"

python "$EXAMPLES_DIR/create_pipeline_script_example.py" \
  --sra "$EXAMPLES_DIR/SRA_metadata_example.xlsx" \
  --data-folder "$EXAMPLES_DIR/data/fastq" \
  --config-yml "$EXAMPLES_DIR/pipeline_settings.yml" \
  --conda-tools-yml "$EXAMPLES_DIR/conda_tools.yml" \
  --out-folder "$RESULT_DIR" \
  --scheduler "$SCHEDULER" \
  --threads "$THREADS"

bash "$RESULT_DIR/submit.sh"
