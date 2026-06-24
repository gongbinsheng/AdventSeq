"""Build and submit per-sample AdVentSeq pipeline scripts for a whole study.

This is an *example* driver script. It reads a study metadata sheet (Excel),
and for every sample it generates one slurm batch script that runs the
AdVentSeq viral-discovery pipeline. It then writes a `submit.sh` that submits
all of those scripts to an HPC cluster with `sbatch`.

Pipeline stages built for each sample (single aligner: BWA-MEM):
    1. MultiQC            - QC report of the raw reads
    2. merge_lanes        - merge multi-lane FASTQ files per read
    3. fastp              - adapter / quality / polyG / polyX trimming
    4. BWA_MEM (host)     - align to the host genome (e.g. hg38)
    5. featureCounts      - gene-level quantification of host reads
    6. remove host reads  - keep only the read pairs that did NOT map to host
    7. Kraken2            - taxonomic screen of the remaining reads
    8. BWA_MEM (RVDB)     - align remaining reads to the RVDB viral database,
                            classify and quantify (TaxonomyClassifier + ViraQuant)
    9. BWA_MEM (7viruses) - align remaining reads to a small curated virus set,
                            classify and quantify (TaxonomyClassifier + ViraQuant)
   10. GATK4 Mutect2      - call variants in the curated viruses (e.g. Zika)

Expected metadata sheet (sheet name "SRA_data", same layout as an SRA export):
    Required columns: library_ID, library_source, library_layout, platform, filetype
    One or more file columns: filename, filename2, filename3, ...
      (each holding a data file name found under --data-folder)
    Optional column: library_selection  (used here only as an adapter example)

Usage:
    python create_pipeline_script.py \\
        -s study_metadata.xlsx \\
        -d /path/to/DATA \\
        -o /path/to/output \\
        --config-yml pipeline_settings.yml

    # then, on the HPC login node:
    bash /path/to/output/submit.sh
"""

import os
import sys
import shutil
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path

from AdVentSeq import Pipeline, list_files_with_extensions


# ============================================================================
# EDIT THESE PATHS / SETTINGS for your environment.
#
# Every absolute path below points to a reference resource on your cluster.
# Replace the "/path/to/..." placeholders with your own paths before running.
#
# Note: index / GTF / FASTA paths can alternatively be configured in the
# `reference_paths:` block of pipeline_settings.yml (keyed by ref_genome), in
# which case BWA_MEM / featureCounts resolve them automatically and you can drop
# the matching arguments below. The contig_info / taxonomy_map / virus_list
# resources are NOT covered by that schema, so they must always be passed here.
# ============================================================================

# Compute resources requested per slurm job.
THREADS = 16
MEMORY = "256G"  # large genome indices can be memory hungry; give jobs headroom.

# Host genome: reads mapping here are quantified, then removed before viral search.
HOST_GENOME = "hg38"
HOST_GTF = "/path/to/hg38/hg38.ncbiRefSeq.gtf"

# Kraken2 database root (used for the taxonomic screen of non-host reads).
KRAKEN2_DB_ROOT = "/path/to/Kraken2DB"

# RVDB: comprehensive viral reference database.
RVDB_VERSION = "v29"
RVDB_BWA_INDEX = f"/path/to/RVDB/{RVDB_VERSION}/BWAIndex/C-RVDB.fa"
RVDB_CONTIG_INFO = f"/path/to/RVDB/{RVDB_VERSION}/U-RVDB{RVDB_VERSION}.sqlite.db"
RVDB_TAXONOMY_MAP = f"/path/to/RVDB/{RVDB_VERSION}/ncbi_taxonomy_map.json.gz"

# Curated small set of viruses of interest (faster, and used for variant calling).
VIRUS_BWA_INDEX = "/path/to/7viruses/BWAIndex/viruses.fa"
VIRUS_CONTIG_INFO = "/path/to/7viruses/contig_info.json.gz"
VIRUS_GENOME_FASTA = "/path/to/7viruses/viruses.fa"
VIRUS_LIST = "/path/to/7viruses/7_virus_list.txt"

# Optional: custom adapter FASTA for libraries prepared with Nextera, etc.
# Leave as None to always use the standard Illumina adapters set below.
NEXTERA_FASTA = None  # e.g. "/path/to/Nextera_adapters.fasta"


def main():
    parser = ArgumentParser(description="Generate and submit AdVentSeq pipeline scripts for a study.")
    parser.add_argument('-s', '--sra', type=str, action='store', dest='sra_metadata', default=None, required=True,
                        help='Study metadata Excel sheet (SRA-style layout, sheet name "SRA_data").')
    parser.add_argument('-d', '--data-folder', type=str, action='store', dest='data_folder', default=None, required=True,
                        help='Path to DATA folder. File names must be unique even across sub-folders.')
    parser.add_argument('-o', '--out-folder', type=str, action='store', dest='wd', default=os.getcwd(), required=False,
                        help='Output folder. Default [current directory: %(default)s]')
    parser.add_argument('--config-yml', type=str, action='store', dest='config_yml', required=True,
                        help='Pipeline configuration file (e.g. pipeline_settings.yml).')

    o = parser.parse_args()

    # ------------------------------------------------------------------
    # Output folders. The working directory must already exist; the script
    # creates two sub-folders inside it:
    #   __slurm : one generated <sample>.slurm batch script per sample
    #   __log   : one <sample>.log file per sample (recreated fresh each run)
    # ------------------------------------------------------------------
    working_dir = Path(o.wd).absolute().resolve()  # normalise relative -> absolute path
    if not working_dir.is_dir():
        sys.exit("Working directory does not exist: [%s] %s\n" % (o.wd, working_dir.as_posix()))

    # Load the YAML config and set cluster-wide defaults shared by every sample.
    Pipeline.load_pipeline_config(o.config_yml)
    Pipeline.threadN = THREADS
    Pipeline.memory = MEMORY
    Pipeline.WD = working_dir.as_posix()

    batch_folder = working_dir.joinpath("__slurm")
    batch_folder.mkdir(parents=True, exist_ok=True)
    log_folder = working_dir.joinpath("__log")
    if log_folder.is_dir():
        shutil.rmtree(log_folder)
    log_folder.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Data folder. Index every .fastq/.fq/.gz/.bam file by its base name so we
    # can look up each sample's files by the name given in the metadata sheet.
    # ------------------------------------------------------------------
    data_folder = Path(o.data_folder).absolute().resolve()
    if not data_folder.is_dir():
        sys.exit("Data folder does not exist: [%s] %s\n" % (o.data_folder, data_folder.as_posix()))
    Pipeline.Data_folder = data_folder.as_posix()
    files = list_files_with_extensions(folder_path=data_folder.as_posix(), extensions=[".fastq", ".fq", ".gz", ".bam"])
    file_dict = {Path(filepath).name: Path(filepath).as_posix() for filepath in files}

    # Standard Illumina TruSeq adapters, used by fastp for every sample.
    # (A per-sample override is shown below via p.adapter_fasta.)
    Pipeline.adapter = {"R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                        "R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}

    # The submit script accumulates one `sbatch` line per sample that still
    # has work to do (samples already fully processed are skipped).
    submit = open(working_dir.joinpath("submit.sh"), 'w')
    num_jobs = 0

    # ------------------------------------------------------------------
    # Read the study metadata and validate that the required columns exist.
    # ------------------------------------------------------------------
    df = pd.read_excel(o.sra_metadata, sheet_name="SRA_data", header=0)
    df = df.dropna(how='all')  # drop fully empty rows
    for col_name in ["library_ID", "library_source", "library_layout", "platform", "filetype"]:
        if col_name not in df.columns:
            sys.exit('The "%s" column is NOT in the header.' % col_name)
    # File columns are named filename, filename2, filename3, ...
    filename_columns = [col for col in df.columns if col.startswith('filename')]
    if not filename_columns:
        sys.exit("No 'filename' columns are present in the header.")

    # ------------------------------------------------------------------
    # Build one pipeline / slurm script per sample (one row = one sample).
    # ------------------------------------------------------------------
    for index, row in df.iterrows():
        sample_id = row['library_ID']
        sys.stdout.write("%s\n" % sample_id)

        # Resolve this sample's data file names to absolute paths.
        file_list = [file_dict[row[col]] for col in filename_columns if pd.notna(row[col])]

        batch_file = batch_folder.joinpath(sample_id + '.slurm')
        log_file = log_folder.joinpath(sample_id + '.log')

        # Create the pipeline object for this sample. ref_genome is the host
        # genome the reads are first aligned (and quantified) against.
        p = Pipeline(sample_id=sample_id,
                     file_list=file_list,
                     log_file=log_file.as_posix(),
                     ref_genome=HOST_GENOME,
                     ref_type="genome",
                     feature="genome",
                     library_source=row["library_source"] if "library_source" in row else "UNKNOWN",
                     library_layout=row["library_layout"] if "library_layout" in row else "UNKNOWN",
                     platform=row["platform"] if "platform" in row else "UNKNOWN",
                     filetype=row["filetype"] if "filetype" in row else "UNKNOWN")

        # Set force=True to re-run steps that have already completed (the
        # pipeline otherwise skips any step whose output files already exist).
        force = False

        # Write the slurm header (#SBATCH lines) and set up the shell ENV.
        p.set_slurm_parameters()
        p.set_env_variables()

        # --- QC and read preparation --------------------------------------
        p.MultiQC()
        p.merge_lanes()

        # Optional per-sample adapter override (example): some library kits
        # (e.g. Nextera) need different adapters than the default TruSeq set.
        # adapter_fasta takes precedence over Pipeline.adapter inside fastp().
        if NEXTERA_FASTA and row.get("library_selection") == "other":
            p.adapter_fasta = NEXTERA_FASTA
        p.fastp()

        # --- Host alignment, quantification, and host-read removal ---------
        # Align all reads to the host genome.
        p.BWA_MEM(force=force)
        p.SAM2BAM(force=force)            # SAM -> unsorted BAM (used for read removal)
        p.BAM_stat(force=force)           # alignment summary statistics

        # Quantify host gene expression (featureCounts needs a coord-sorted BAM).
        p.sort_BAM(by_qname=False, force=force)
        p.featureCounts(gene_model="hg38_ncbiRefSeq", feature="exon", countReadPairs=True, gtf=HOST_GTF, force=force)
        p.featureCounts(gene_model="hg38_ncbiRefSeq", feature="exon", countReadPairs=False, gtf=HOST_GTF, force=force)

        # Drop read pairs that mapped to the host; only the rest go to viral search.
        p.remove_read_pairs_mapped_to_host(use_samtools=True, force=force)

        # Quick taxonomic screen of the remaining (non-host) reads.
        p.Kraken2(db_name="viral", kraken2_db_root=KRAKEN2_DB_ROOT, force=force)

        # --- Viral search 1: RVDB (comprehensive viral database) ----------
        p.BWA_MEM(ref_genome=f"RVDB{RVDB_VERSION}", bwa_index=RVDB_BWA_INDEX, force=force)
        p.SAM2BAM(force=force)
        p.TaxonomyClassifier(contig_info=RVDB_CONTIG_INFO,
                             taxonomy_map=RVDB_TAXONOMY_MAP,
                             virus_group_by="ncbitaxon",
                             force=force)
        p.sort_BAM(by_qname=False, force=force)
        p.ViraQuant(scan_by="mean_depth_ge_1>=1", force=force)

        # --- Viral search 2: curated 7-virus set --------------------------
        p.BWA_MEM(ref_genome="7viruses", bwa_index=VIRUS_BWA_INDEX, force=force)
        p.SAM2BAM(force=force)
        p.TaxonomyClassifier(contig_info=VIRUS_CONTIG_INFO, force=force)
        p.sort_BAM(by_qname=False, force=force)
        p.ViraQuant(virus_list=VIRUS_LIST, force=force)

        # --- Variant calling in the curated viruses (e.g. Zika) -----------
        p.sort_BAM(force=force)
        p.GATK4_Mutect2(ref_genome="7viruses",
                        genome_fasta=VIRUS_GENOME_FASTA,
                        no_filter=True,
                        split_multi_allelic=True,
                        force=force)

        # Write all accumulated step commands to this sample's slurm script.
        with open(batch_file, 'w') as batch:
            for step in p.batch:
                batch.write("".join(p.batch[step]))

        # Only submit samples that still have unfinished steps.
        if not p.is_completed:
            submit.write('sbatch "%s"\n' % batch_file.as_posix())
            num_jobs += 1

    submit.close()

    if num_jobs == 0:
        sys.stdout.write("No job to be run, ALL done!\n\n")
    else:
        sys.stdout.write(f"{num_jobs} jobs need to be run. please run:\nbash {working_dir}/submit.sh\n\n")


if __name__ == "__main__":
    main()
