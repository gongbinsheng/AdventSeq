conda activate AdventSeq

# Collect the per-sample pipeline outputs (produced by run_example.sh) into tables.
# RESULT_DIR must be the same directory passed to the pipeline via -o.
RESULT_DIR="Results"
TABLES_DIR="Collected_Tables"

collect-results \
  --sra SRA_metadata_example.xlsx \
  --results-dir ./"$RESULT_DIR" \
  --out-dir ./"$TABLES_DIR" \
  --scan-ref RVDBv29 \
  --targeted-refs 7viruses

# Open ./"$TABLES_DIR"/viraquant_targeted_interactive.html in a browser for the
# interactive coverage table (measurement drop-down).

conda deactivate
