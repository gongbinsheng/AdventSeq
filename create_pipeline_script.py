import os
import sys
import shutil
import pandas as pd
from argparse import ArgumentParser
from AdVentSeq import Pipeline, list_files_with_extensions
from HPC import HPC
from pathlib import Path

def main():
    parser = ArgumentParser()
    parser.add_argument('-s', '--sra', type=str, action='store', dest='sra_metadata', default=None, required=True,
                        help='SRA metadata sheet.')
    parser.add_argument('-d', '--data-folder', type=str, action='store', dest='data_folder', default=None, required=True,
                        help='Path to DATA folder. Data file names should be unique no matter if they are in different folders.')
    parser.add_argument('-o', '--out-folder', type=str, action='store', dest='wd', default=os.getcwd(), required=False,
                        help='output folder. Default [Current Directory: %(default)s]')
    parser.add_argument('--config-yml', type=str, action='store', dest='config_yml', required=True,
                        help='Pipeline configuration file')

    o = parser.parse_args()

    # Check working directory and create batch and log folders in it
    # Working directory must already exist.
    working_dir = Path(o.wd).absolute().resolve() # convert relative path to absolute path
    if not working_dir.is_dir():
        sys.exit("Working directory does not exist: [%s] %s\n" % (o.wd, working_dir.as_posix()))

    hpc = HPC()
    #hpc.add_hpc_nodes_to_be_skipped("ncshpc410","ncshpc411")

    Pipeline.load_pipeline_config(o.config_yml)
    Pipeline.threadN = 16
    Pipeline.memory = "256G" # when using 32G, star with kill itself at "inserting junctions into the genome indices"
    Pipeline.WD = working_dir.as_posix()
    batch_folder = working_dir.joinpath("__slurm")
    batch_folder.mkdir(parents=True, exist_ok=True)
    log_folder = working_dir.joinpath("__log")
    if log_folder.is_dir():
        shutil.rmtree(log_folder)
    log_folder.mkdir(parents=True, exist_ok=True)

    # Check data folder and look for .gz and .bam files. Data folder must exist.
    data_folder = Path(o.data_folder).absolute().resolve() # convert relative path to absolute path
    if not data_folder.is_dir():
        sys.exit("Data folder does not exist: [%s] %s\n" % (o.data_folder, data_folder.as_posix()))
    Pipeline.Data_folder = data_folder.as_posix()
    files = list_files_with_extensions(folder_path=data_folder.as_posix(), extensions=[".fastq", ".fq", ".gz", ".bam"])
    file_dict = {}
    for filepath in files:
        filepath = Path(filepath)
        file_dict[filepath.name] = filepath.as_posix()


    # Set adapter sequences
    #Pipeline.adapter_fasta = "adapters.fasta"
    # Read 1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA Read 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    Pipeline.adapter = {"R1":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "R2":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}
    # initiate a script. Add submission commands to it later
    submit = open(working_dir.joinpath("submit.sh"),'w')
    num_jobs = 0

    ### Read the matadata
    df = pd.read_excel(o.sra_metadata, sheet_name="SRA_data", header=0) # same format as SRA
    df = df.dropna(how='all') # Remove empty rows
    # Check if the header contains required columns
    for col_name in ["library_ID", "library_source", "library_layout", "platform", "filetype"]:
        if col_name not in df.columns:
            sys.exit('The "%s" column is NOT in the header.' % col_name)
    # Check for 'filename', 'filename2', 'filename3', etc.
    filename_columns = [col for col in df.columns if col.startswith('filename')]
    if not filename_columns:
        sys.exit("No 'filename' columns are present in the header.")

    ### Loop through each row
    for index, row in df.iterrows():
        # assign values to sample_ID and filenames
        sample_id = row['library_ID']
        sys.stdout.write("%s\n" % sample_id)
        file_list = [file_dict[row[col]] for col in filename_columns if pd.notna(row[col])]
        #file_list.sort()

        # assign filenames for batch and log files
        batch_file = batch_folder.joinpath(sample_id+'.slurm')
        log_file = log_folder.joinpath(sample_id+'.log')

        # create an instance of the Pipeline class
        #ref_genome = "virtual_genome"
        ref_genome = "hg38"
        p = Pipeline(sample_id=sample_id,
                     file_list=file_list,
                     log_file=log_file.as_posix(),
                     ref_genome=ref_genome,
                     ref_type="genome",
                     feature="genome",
                     library_source=row["library_source"] if "library_source" in row else "UNKNOWN",
                     library_layout = row["library_layout"] if "library_layout" in row else "UNKNOWN",
                     platform = row["platform"] if "platform" in row else "UNKNOWN",
                     filetype = row["filetype"] if "filetype" in row else "UNKNOWN")

        # run the pipeline step by step
        p.set_slurm_parameters()
        p.set_env_variables()

        p.MultiQC()
        p.merge_lanes()
        # trim adapter
        if row["library_selection"] == "other":
            p.adapter_fasta = "/galaxy001/bgong/Project_DeLab/Nextera_adapters.fasta"
            # note: the adapter_fasta should overwrite the Pipeline.adapter in the fastp() below
        p.fastp()

        ########################################################
        ### remove host reads then map to RVDB and 7 viruses ###
        ########################################################
        force = False
        host_genome = "hg38"
        #host_mappers = ("BWA_MEM","HISAT2","STAR")
        #host_mappers = ("BWA_MEM","STAR")
        host_mappers = ("BWA_MEM", )
        for host_mapper in host_mappers:
            if len(host_mappers) > 1:
                ### reset to initial values
                p.reset_pipeline(ref_genome=host_genome)

            ### map to host
            if host_mapper == "BWA_MEM":
                p.BWA_MEM(force=force)
            elif host_mapper == "HISAT2":
                p.HISAT2(force=force)
            elif host_mapper == "STAR":
                p.STAR(force=force)
            else:
                sys.exit(f"Wrong host mapper: {host_mapper}")
            # covert sam to unsorted bam for remove host read pairs
            p.SAM2BAM(force=force)
            p.BAM_stat(force=force)

            # quantification for host
            p.sort_BAM(by_qname=False, force=force) # sort by coordinate for featureCounts
            gtf_file = "/galaxy001/Resources/refGenomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2023-01-06/hg38.ncbiRefSeq.gtf"
            p.featureCounts(gene_model="hg38_ncbiRefSeq", feature="exon", countReadPairs=True, gtf=gtf_file, force=force)
            p.featureCounts(gene_model="hg38_ncbiRefSeq", feature="exon", countReadPairs=False, gtf=gtf_file, force=force)

            # remove host reads using unsorted bam
            p.remove_read_pairs_mapped_to_host(use_samtools=True, force=force)

            # get step id after removing host read-pairs
            virus_start_step_id = p.get_current_step_id()

            # Kraken2
            p.Kraken2(db_name="viral",kraken2_db_root="/galaxy001/Resources/Kraken2DB",force=force)

            # whether the unmapped reads can map to host using STAR
            #p.STAR(force=force)
            #p.SAM2BAM(force=force)

            ### map the rest read pairs to viruses
            #virus_mappers = ("BWA_MEM", "Bowtie2")
            virus_mappers = ("BWA_MEM",)
            for virus_mapper in virus_mappers:
                # map the rest reads to RVDB
                RVDB_version = "v29" # "v31"
                p.set_current_step_id(virus_start_step_id)  # reset step id for each aligner
                if virus_mapper == "Bowtie2":
                    p.Bowtie2(ref_genome=f"RVDB{RVDB_version}",
                              bowtie2_index=f"/galaxy001/Resources/RVDB/{RVDB_version}.0/Bowtie2/C-RVDB",
                              force=force)
                elif virus_mapper == "BWA_MEM":
                    p.BWA_MEM(ref_genome=f"RVDB{RVDB_version}",
                              bwa_index=f"/galaxy001/Resources/RVDB/{RVDB_version}.0/BWAIndex/C-RVDB.fa",
                              force=force)
                else:
                    sys.exit(f"Wrong virus mapper: {virus_mapper}")
                # # count reads for RVDB
                p.SAM2BAM(force=force)
                p.TaxonomyClassifier(contig_info=f"/galaxy001/Resources/RVDB/{RVDB_version}.0/U-RVDB{RVDB_version}.0.sqlite.db",
                                       taxonomy_map=f"/galaxy001/Resources/RVDB/{RVDB_version}.0/ncbi_taxonomy_map.json.gz",
                                       virus_group_by="ncbitaxon",
                                       force=force)
                p.sort_BAM(by_qname=False, force=force)

                p.ViraQuant(scan_by="mean_depth_ge_1>=1", force=force)

                # map the rest reads to viruses
                p.set_current_step_id(virus_start_step_id)  # reset step id for each aligner
                if virus_mapper == "Bowtie2":
                    p.Bowtie2(ref_genome="7viruses",
                              bowtie2_index="/galaxy001/bgong/Project_DeLab/Hg38plus7viruses/Bowtie2Index/viruses",
                              force=force)
                elif virus_mapper == "BWA_MEM":
                    p.BWA_MEM(ref_genome="7viruses",
                              bwa_index="/galaxy001/bgong/Project_DeLab/Hg38plus7viruses/BWAIndex/viruses.fa",
                              force=force)
                else:
                    sys.exit(f"Wrong virus mapper: {virus_mapper}")
                # count reads for viruses
                p.SAM2BAM(force=force)
                p.TaxonomyClassifier(contig_info="/galaxy001/bgong/Project_DeLab/Hg38plus7viruses/WholeGenomeFasta/contig_info.json.gz",
                                       force=force)
                p.sort_BAM(by_qname=False, force=force)
                p.ViraQuant(virus_list="/galaxy001/bgong/Project_DeLab/Virus_genome/7_virus_list.txt", force=force)


                # call mutation in viruses (especially for Zika)
                p.sort_BAM(force=force)
                p.GATK4_Mutect2(ref_genome="7viruses",
                                genome_fasta="/galaxy001/bgong/Project_DeLab/Hg38plus7viruses/WholeGenomeFasta/viruses.fa",
                                no_filter=True,
                                split_multi_allelic=True,
                                force=force)

        # write commands to batch file
        with open(batch_file, 'w') as batch:
            for step in p.batch:
                batch.write("".join(p.batch[step]))  # batch options


        # add submission command to a script
        if not p.is_completed:
            if len(hpc.hpc_nodes_to_be_skipped) > 0:
                hpc_nodes_to_be_skipped = ",".join(hpc.get_HPC_nodes_to_be_skipped())
                exclude_nodes_option = f"--exclude={hpc_nodes_to_be_skipped} "
            else:
                exclude_nodes_option = ""
            submit.write('sbatch %s"%s"\n' % (exclude_nodes_option, batch_file.as_posix()))
            num_jobs += 1

    if num_jobs == 0:
        sys.stdout.write("No job to be run, ALL done!\n\n")
    else:
        sys.stdout.write(f"{num_jobs} jobs need to be run. please run:\nbash {working_dir}/submit.sh\n\n")


if __name__ == "__main__":
    main()
