

conda activate AdVentSeq
FASTQ_DATA_DIR="<folder of FASTQ data>"
RESULT_DIR="Results"
mkdir -p "$result_dir"
python create_pipeline_script_example.py \
  --sra SRA_metadata_example.xlsx \
  -d "$FASTQ_DATA_DIR" \
  --config-yml "pipeline_settings.yml" \
  -o ./"$RESULT_DIR"

conda deactivate
