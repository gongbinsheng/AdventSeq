#!/usr/bin/env bash
# Collect the per-sample outputs written by run_example.sh into tables.
#
# Run after all samples of run_example.sh have finished, with the AdventSeq
# env active:
#   bash examples/collect_results_example.sh
set -euo pipefail

EXAMPLES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT_DIR="$EXAMPLES_DIR/Results"          # the --out-folder used by run_example.sh
TABLES_DIR="$EXAMPLES_DIR/Collected_Tables"

collect-results \
  --sra "$EXAMPLES_DIR/SRA_metadata_example.xlsx" \
  --results-dir "$RESULT_DIR" \
  --out-dir "$TABLES_DIR" \
  --scan-ref RVDBv29 \
  --targeted-refs virus_panel \
  --taxonomy-map "$EXAMPLES_DIR/references/RVDB/v29.0/ncbi_taxonomy_map.json.gz" \
  --host-genes "$EXAMPLES_DIR/host_gene_families_example.yml"

echo "Tables written to $TABLES_DIR"
echo "Open $TABLES_DIR/viraquant_targeted_interactive.html in a browser for the interactive coverage table."
