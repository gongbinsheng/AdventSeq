"""Generate per-sample AdventSeq job scripts for the example study.

This example driver reads the study metadata sheet (SRA_metadata_example.xlsx)
and, for every sample, writes one bash job script that runs the AdventSeq
viral-discovery pipeline. It then writes `submit.sh`, which either runs all
job scripts one after another on this computer (--scheduler local, default) or
submits them to a Slurm cluster with `sbatch` (--scheduler slurm).

Pipeline stages built for each sample (single aligner: BWA-MEM):
    1. MultiQC            - FastQC + MultiQC report of the raw reads
    2. merge_lanes        - merge multi-lane FASTQ files per read
    3. fastp              - adapter / quality / polyG / polyX trimming
    4. BWA_MEM (host)     - align to the host genome (hg38)
    5. featureCounts      - gene-level quantification of host reads
    6. remove host reads  - keep only the read pairs that did NOT map to host
    7. Kraken2            - taxonomic screen of the remaining reads
    8. BWA_MEM (RVDB)     - align remaining reads to the RVDB viral database,
                            classify and quantify (TaxonomyClassifier + ViraQuant)
    9. BWA_MEM (virus_panel) - align remaining reads to a small virus panel,
                            classify and quantify (TaxonomyClassifier + ViraQuant)
   10. GATK4 Mutect2      - call variants in the panel viruses

Input data are NOT shipped with the repository. Before running this script,
prepare them as described in the "Running the example" section of README.md:
    examples/data/fastq/     - your FASTQ files, listed in the metadata sheet
    examples/references/     - host genome, Kraken2 DB, RVDB and virus panel

Usage (from the repository root, with the AdventSeq conda env active):
    python examples/create_pipeline_script_example.py
    bash examples/Results/submit.sh
"""

import os
import sys
import shutil
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from AdventSeq import Pipeline, list_files_with_extensions


EXAMPLES_DIR = Path(__file__).resolve().parent
REF_DIR = EXAMPLES_DIR / "references"

# ============================================================================
# Reference resources. The defaults match the layout produced by the
# preparation steps in README.md; edit them to use files stored elsewhere.
# All paths must be absolute (they are written into the job scripts).
# ============================================================================

# Host genome: reads mapping here are quantified, then removed before viral search.
HOST_GENOME = "hg38"
HOST_BWA_INDEX = REF_DIR / "hg38" / "BWAIndex" / "genome.fa"
HOST_GTF = REF_DIR / "hg38" / "hg38.ncbiRefSeq.gtf"
HOST_GENE_MODEL = "hg38_ncbiRefSeq"

# Kraken2 database root; the database itself is <root>/<name>.
KRAKEN2_DB_ROOT = REF_DIR / "Kraken2DB"
KRAKEN2_DB_NAME = "viral"

# RVDB: comprehensive viral reference database (https://rvdb.dbi.udel.edu).
RVDB_RELEASE = "v29.0"
RVDB_NAME = "RVDBv29"  # reference name used in output file names (collect-results --scan-ref)
RVDB_DIR = REF_DIR / "RVDB" / RVDB_RELEASE
RVDB_BWA_INDEX = RVDB_DIR / "BWAIndex" / "C-RVDB.fa"
RVDB_CONTIG_INFO = RVDB_DIR / "contig_info.json.gz"
RVDB_TAXONOMY_MAP = RVDB_DIR / "ncbi_taxonomy_map.json.gz"

# Small panel of viruses of interest (faster, and used for variant calling).
VIRUS_PANEL = "virus_panel"  # reference name used in output file names (collect-results --targeted-refs)
VIRUS_DIR = REF_DIR / VIRUS_PANEL
VIRUS_GENOME_FASTA = VIRUS_DIR / "viruses.fa"
VIRUS_BWA_INDEX = VIRUS_DIR / "BWAIndex" / "viruses.fa"
VIRUS_CONTIG_INFO = VIRUS_DIR / "contig_info.json"
VIRUS_LIST = VIRUS_DIR / "virus_list.txt"

BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")


def required_reference_files():
    """Every reference file the example needs, for a friendly up-front check."""
    files = [HOST_GTF, RVDB_CONTIG_INFO, RVDB_TAXONOMY_MAP,
             VIRUS_GENOME_FASTA, VIRUS_GENOME_FASTA.with_suffix(".fa.fai"),
             VIRUS_GENOME_FASTA.with_suffix(".dict"), VIRUS_CONTIG_INFO, VIRUS_LIST]
    for prefix in (HOST_BWA_INDEX, RVDB_BWA_INDEX, VIRUS_BWA_INDEX):
        files += [Path(f"{prefix}{suffix}") for suffix in BWA_INDEX_SUFFIXES]
    files.append(KRAKEN2_DB_ROOT / KRAKEN2_DB_NAME / "hash.k2d")
    return files


def find_conda_init_script():
    """Locate <conda base>/etc/profile.d/conda.sh from the active conda install."""
    candidates = []
    if os.environ.get("CONDA_EXE"):
        candidates.append(Path(os.environ["CONDA_EXE"]))
    if shutil.which("conda"):
        candidates.append(Path(shutil.which("conda")))
    for conda_exe in candidates:
        init_script = conda_exe.resolve().parents[1] / "etc" / "profile.d" / "conda.sh"
        if init_script.is_file():
            return init_script
    return None


def main():
    parser = ArgumentParser(description="Generate AdventSeq pipeline job scripts for the example study.")
    parser.add_argument('-s', '--sra', dest='sra_metadata', default=EXAMPLES_DIR / "SRA_metadata_example.xlsx",
                        help='Study metadata Excel sheet (SRA-style layout, sheet name "SRA_data"). '
                             'Default: %(default)s')
    parser.add_argument('-d', '--data-folder', dest='data_folder', default=EXAMPLES_DIR / "data" / "fastq",
                        help='Folder holding the FASTQ files. File names must be unique even across '
                             'sub-folders. Default: %(default)s')
    parser.add_argument('-o', '--out-folder', dest='wd', default=EXAMPLES_DIR / "Results",
                        help='Output folder (created if missing). Default: %(default)s')
    parser.add_argument('--config-yml', dest='config_yml', default=EXAMPLES_DIR / "pipeline_settings.yml",
                        help='Pipeline configuration file. Default: %(default)s')
    parser.add_argument('--conda-tools-yml', dest='conda_tools_yml', default=EXAMPLES_DIR / "conda_tools.yml",
                        help='Conda env definitions (step -> env). Default: %(default)s')
    parser.add_argument('--conda-init', dest='conda_init', default=None,
                        help='Conda init script sourced by the job scripts, usually '
                             '<conda base>/etc/profile.d/conda.sh. Auto-detected when conda is active.')
    parser.add_argument('--scheduler', choices=("local", "slurm"), default="local",
                        help='"local": submit.sh runs the samples one after another on this computer; '
                             '"slurm": submit.sh submits them with sbatch. Default: %(default)s')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Threads per job. Default: %(default)s')
    parser.add_argument('--skip-conda-check', action='store_true',
                        help='Do not verify the conda envs or record tool versions (only for inspecting '
                             'the generated scripts; the jobs still need the envs to run).')
    o = parser.parse_args()

    # ------------------------------------------------------------------
    # Check the inputs up front so a missing file gives a clear message.
    # ------------------------------------------------------------------
    missing = [p for p in required_reference_files() if not p.exists()]
    if missing:
        sys.exit("Missing reference files (see 'Running the example' in README.md):\n  "
                 + "\n  ".join(p.as_posix() for p in missing) + "\n")

    data_folder = Path(o.data_folder).absolute().resolve()
    if not data_folder.is_dir():
        sys.exit("Data folder does not exist: %s\n" % data_folder.as_posix())

    conda_init = Path(o.conda_init) if o.conda_init else find_conda_init_script()
    if conda_init is None:
        sys.exit("Could not find the conda init script. Activate conda first, or pass\n"
                 "  --conda-init <conda base>/etc/profile.d/conda.sh\n")

    # ------------------------------------------------------------------
    # Pipeline-wide settings shared by every sample.
    # ------------------------------------------------------------------
    Pipeline.set_conda_init_script(conda_init)
    Pipeline.set_conda_tools_config(o.conda_tools_yml)  # must come before load_pipeline_config
    Pipeline.load_pipeline_config(o.config_yml, check_conda_setup=not o.skip_conda_check)
    Pipeline.threadN = o.threads

    # Output folders: one job script per sample in __jobs, one log per sample in __log.
    working_dir = Path(o.wd).absolute().resolve()
    working_dir.mkdir(parents=True, exist_ok=True)
    Pipeline.WD = working_dir.as_posix()
    batch_folder = working_dir.joinpath("__jobs")
    batch_folder.mkdir(parents=True, exist_ok=True)
    log_folder = working_dir.joinpath("__log")
    if log_folder.is_dir():
        shutil.rmtree(log_folder)
    log_folder.mkdir(parents=True, exist_ok=True)

    # Index every .fastq/.fq/.gz/.bam file by its base name so each sample's
    # files can be looked up by the names given in the metadata sheet.
    Pipeline.Data_folder = data_folder.as_posix()
    files = list_files_with_extensions(folder_path=data_folder.as_posix(), extensions=[".fastq", ".fq", ".gz", ".bam"])
    file_dict = {Path(filepath).name: Path(filepath).as_posix() for filepath in files}

    # Standard Illumina TruSeq adapters, used by fastp for every sample.
    Pipeline.adapter = {"R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                        "R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}

    # ------------------------------------------------------------------
    # Read the study metadata and validate that the required columns exist.
    # ------------------------------------------------------------------
    df = pd.read_excel(o.sra_metadata, sheet_name="SRA_data", header=0)
    df = df.dropna(how='all')  # drop fully empty rows
    for col_name in ["library_ID", "library_source", "library_layout", "platform", "filetype"]:
        if col_name not in df.columns:
            sys.exit('The "%s" column is NOT in the header.' % col_name)
    # File columns are named filename, filename1, filename2, ...
    filename_columns = [col for col in df.columns if col.startswith('filename')]
    if not filename_columns:
        sys.exit("No 'filename' columns are present in the header.")

    missing_fastq = [row[col] for _, row in df.iterrows() for col in filename_columns
                     if pd.notna(row[col]) and row[col] not in file_dict]
    if missing_fastq:
        sys.exit("FASTQ files listed in %s were not found in %s:\n  %s\n"
                 % (Path(o.sra_metadata).name, data_folder.as_posix(), "\n  ".join(missing_fastq)))

    submit = open(working_dir.joinpath("submit.sh"), 'w')
    submit.write("#!/usr/bin/env bash\n")
    num_jobs = 0

    # ------------------------------------------------------------------
    # Build one job script per sample (one row = one sample).
    # ------------------------------------------------------------------
    for _, row in df.iterrows():
        sample_id = row['library_ID']
        sys.stdout.write("%s\n" % sample_id)

        file_list = [file_dict[row[col]] for col in filename_columns if pd.notna(row[col])]
        batch_file = batch_folder.joinpath(sample_id + ('.slurm' if o.scheduler == "slurm" else '.sh'))
        log_file = log_folder.joinpath(sample_id + '.log')

        # ref_genome is the host genome the reads are first aligned (and quantified) against.
        p = Pipeline(sample_id=sample_id,
                     file_list=file_list,
                     log_file=log_file.as_posix(),
                     ref_genome=HOST_GENOME,
                     ref_type="genome",
                     feature="genome",
                     library_source=row["library_source"],
                     library_layout=row["library_layout"],
                     platform=row["platform"],
                     filetype=row["filetype"])

        # Set force=True to re-run steps that have already completed (the
        # pipeline otherwise skips any step recorded in <sample>/_my_progress).
        force = False

        # Job header (#SBATCH lines are plain comments when run with bash) + shell ENV.
        p.set_slurm_parameters()
        p.set_env_variables()

        # --- QC and read preparation --------------------------------------
        p.MultiQC()
        p.merge_lanes()
        p.fastp()

        # --- Host alignment, quantification, and host-read removal ---------
        p.BWA_MEM(ref_genome=HOST_GENOME, bwa_index=HOST_BWA_INDEX.as_posix(), force=force)
        p.SAM2BAM(force=force)            # SAM -> unsorted BAM (used for read removal)
        p.BAM_stat(force=force)           # alignment summary statistics

        # Quantify host gene expression (featureCounts needs a coord-sorted BAM).
        p.sort_BAM(by_qname=False, force=force)
        for count_pairs in (True, False):
            p.featureCounts(gene_model=HOST_GENE_MODEL, feature="exon", countReadPairs=count_pairs,
                            gtf=HOST_GTF.as_posix(), force=force)

        # Drop read pairs that mapped to the host; only the rest go to viral search.
        p.remove_read_pairs_mapped_to_host(force=force)
        virus_start_step_id = p.get_current_step_id()

        # Quick taxonomic screen of the remaining (non-host) reads.
        p.Kraken2(db_name=KRAKEN2_DB_NAME, kraken2_db_root=KRAKEN2_DB_ROOT.as_posix(), force=force)

        # --- Viral search 1: RVDB (comprehensive viral database) ----------
        p.set_current_step_id(virus_start_step_id)
        p.BWA_MEM(ref_genome=RVDB_NAME, bwa_index=RVDB_BWA_INDEX.as_posix(), force=force)
        p.SAM2BAM(force=force)
        p.TaxonomyClassifier(contig_info=RVDB_CONTIG_INFO.as_posix(),
                             taxonomy_map=RVDB_TAXONOMY_MAP.as_posix(),
                             virus_group_by="ncbitaxon",
                             force=force)
        p.sort_BAM(by_qname=False, force=force)
        p.ViraQuant(scan_by="mean_depth_ge_1>=1", force=force)
        # Group the RVDB scan results to NCBI-taxon level (preferred by collect-results).
        p.convert_ViraQuant_scan_to_ncbitaxon(taxonomy_map=RVDB_TAXONOMY_MAP.as_posix(), force=force)

        # --- Viral search 2: virus panel ---------------------------------
        p.set_current_step_id(virus_start_step_id)
        p.BWA_MEM(ref_genome=VIRUS_PANEL, bwa_index=VIRUS_BWA_INDEX.as_posix(), force=force)
        p.SAM2BAM(force=force)
        p.TaxonomyClassifier(contig_info=VIRUS_CONTIG_INFO.as_posix(), force=force)
        p.sort_BAM(by_qname=False, force=force)
        p.ViraQuant(virus_list=VIRUS_LIST.as_posix(), force=force)

        # --- Variant calling in the panel viruses -------------------------
        p.GATK4_Mutect2(ref_genome=VIRUS_PANEL,
                        genome_fasta=VIRUS_GENOME_FASTA.as_posix(),
                        no_filter=True,
                        split_multi_allelic=True,
                        force=force)

        # Write all accumulated step commands to this sample's job script.
        with open(batch_file, 'w') as batch:
            for step in p.batch:
                batch.write("".join(p.batch[step]))

        # Only run samples that still have unfinished steps.
        if not p.is_completed:
            if o.scheduler == "slurm":
                submit.write('sbatch "%s"\n' % batch_file.as_posix())
            else:
                submit.write('echo "Running %s (log: %s)"\n' % (sample_id, log_file.as_posix()))
                submit.write('bash "%s" > "%s" 2>&1 || echo "  FAILED: %s, see the log above"\n'
                             % (batch_file.as_posix(), log_file.as_posix(), sample_id))
            num_jobs += 1

    submit.close()

    if num_jobs == 0:
        sys.stdout.write("No job to be run, ALL done!\n\n")
    else:
        sys.stdout.write(f"{num_jobs} jobs need to be run. Please run:\n  bash {working_dir}/submit.sh\n\n")


if __name__ == "__main__":
    main()
