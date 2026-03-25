import sys
import os
import pysam
import gzip
import json
import shutil
import sqlite3
import tempfile
from collections import defaultdict
from argparse import ArgumentParser
from pathlib import Path

#from get_info_from_Entrez import contig_json


def normalize_record(columns, row):
    record = dict(zip(columns, row, strict=True))
    accession = record.pop("accs")

    seqlen = record.get("seqlen")
    if isinstance(seqlen, str):
        try:
            record["seqlen"] = int(seqlen)
        except ValueError:
            pass

    return accession, record


def load_contig_info_from_json(contig_info_path):
    suffixes = contig_info_path.suffixes
    if suffixes[-2:] == [".json", ".gz"]:
        with gzip.open(contig_info_path, "rt", encoding="utf-8") as handle:
            return json.load(handle)

    with contig_info_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_contig_info_from_sqlite_db(db_path):
    contig_info = {}
    with sqlite3.connect(db_path) as connection:
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM rvdb")
        columns = [description[0] for description in cursor.description]

        for row in cursor:
            accession, record = normalize_record(columns, row)
            contig_info[accession] = record

    return contig_info


def load_contig_info_from_sqlite(contig_info_path):
    suffixes = contig_info_path.suffixes
    if suffixes[-1:] == [".gz"]:
        with tempfile.NamedTemporaryFile(suffix=".sqlite.db", delete=False) as tmp_handle:
            tmp_path = Path(tmp_handle.name)

        try:
            with gzip.open(contig_info_path, "rb") as source, tmp_path.open("wb") as destination:
                shutil.copyfileobj(source, destination)
            return load_contig_info_from_sqlite_db(tmp_path)
        finally:
            tmp_path.unlink(missing_ok=True)

    return load_contig_info_from_sqlite_db(contig_info_path)


def load_contig_info(contig_info_filepath):
    contig_info_path = Path(contig_info_filepath)
    suffixes = contig_info_path.suffixes

    if suffixes[-2:] == [".json", ".gz"] or suffixes[-1:] == [".json"]:
        return load_contig_info_from_json(contig_info_path)

    if suffixes[-3:] == [".sqlite", ".db", ".gz"] or suffixes[-2:] == [".db", ".gz"]:
        return load_contig_info_from_sqlite(contig_info_path)

    if suffixes[-2:] == [".sqlite", ".db"] or suffixes[-1:] == [".db"]:
        return load_contig_info_from_sqlite(contig_info_path)

    sys.exit(
        f"Unsupported file type for {contig_info_filepath}. "
        "Supported formats: .json, .json.gz, .db, .db.gz, .sqlite.db, .sqlite.db.gz"
    )


def main():
    parser = ArgumentParser()
    parser.add_argument('--bam', type=str, action='store', dest='in_sam', default=None, required=True,
                        help='BAM file sorted by query name')
    parser.add_argument('--combined_genome', type=str, action='store', dest='combined_genome', default=None, required=True,
                        help='name of the combined genome, used for identify contig json file')
    parser.add_argument('--contig_info', type=str, action='store', dest='contig_info', default=None,
                        required=True, help='contig information in JSON or SQLite format')

    o = parser.parse_args()

    if o.contig_info:
        contig_json_filepath = o.contig_info
    else:
        sys.exit(f"contig_info file is not provided for {o.combined_genome}")

    contig_info = load_contig_info(contig_json_filepath)

    sample_id = os.path.splitext(os.path.basename(o.in_sam))[0]

    primary_mapped = defaultdict(lambda: {"R1":[], "R2":[]})
    discordant_mapped = defaultdict(lambda: {"R1":[], "R2":[]})
    unmapped = defaultdict(lambda: {"R1":[], "R2":[]})


    def get_entrez_id(contig_id):
        # acc|GENBANK|HM118302.1|HIV-1 isolate BREPM3072_06 from Brazil envelope glycoprotein (env) gene, partial cds|Human immunodeficiency virus 1|VRL|25-JUL-2016
        # chr1
        # chr1_KI270706v1_random
        strs = contig_id.strip().split("|")
        if len(strs) > 2:
            return strs[2]
        else:
            return strs[0]


    def write_primary_mapped_read_pairs():
        reads_written = set()
        for id, read_pair in primary_mapped.items():
            if len(read_pair["R1"]) == 1 and len(read_pair["R2"]) == 1:
                R1 = read_pair["R1"][0]
                R2 = read_pair["R2"][0]
                if R1.reference_name.startswith("chr"):
                    human_bam.write(R1)
                    human_bam.write(R2)
                    read_pair_count["human"] += 1
                else:
                    viruses_bam.write(R1)
                    viruses_bam.write(R2)
                    read_pair_count["viruses"] += 1
                    virus_count[contig_info[get_entrez_id(R1.reference_name)]['organism']]["pair"] += 1
                read_pair_count["primary"] += 1
                read_pair_count["mapped"] += 1
                read_pair_count["total"] += 1
                reads_written.add(id)
        for id in reads_written:
            primary_mapped.pop(id)


    def write_discordant_read_pairs():
        reads_written = set()
        for id, read_pair in discordant_mapped.items():
            if len(read_pair["R1"]) == 1 and len(read_pair["R2"]) == 1:
                R1 = read_pair["R1"][0]
                R2 = read_pair["R2"][0]
                discordant_bam.write(R1)
                discordant_bam.write(R2)
                read_pair_count["discordant"] += 1
                read_pair_count["mapped"] += 1
                read_pair_count["total"] += 1
                reads_written.add(id)
        for id in reads_written:
            discordant_mapped.pop(id)


    def write_unmapped_read_pairs():
        reads_written = set()
        for id, read_pair in unmapped.items():
            if len(read_pair["R1"]) == 1 and len(read_pair["R2"]) == 1:
                R1 = read_pair["R1"][0]
                R2 = read_pair["R2"][0]
                unmapped_R1.write(f'@{R1.query_name}\n{R1.query_sequence}\n+\n{R1.qual}\n')
                unmapped_R2.write(f'@{R2.query_name}\n{R2.query_sequence}\n+\n{R2.qual}\n')
                read_pair_count["unmapped"] += 1
                read_pair_count["total"] += 1
                reads_written.add(id)
        for id in reads_written:
            unmapped.pop(id)


    def print_progress():
        sys.stdout.write(f"Total: {read_count['total']}/{read_pair_count['total']}\t"
                         f"QCfail: {read_count['qcfail']}\t"
                         f"Mapped: {read_count['mapped']}/{read_pair_count['mapped']}\t"
                         f"Unmapped: {read_count['unmapped']}/{read_pair_count['unmapped']}\t"
                         f"Secondary: {read_count['secondary']}\t"
                         f"Supplementary: {read_count['supplementary']}\t"
                         f"Primary: {read_pair_count['primary']}\t"
                         f"Discordant: {read_pair_count['discordant']}\t"
                         f"Human: {read_count['human']}/{read_pair_count['human']}\t"
                         f"Viruses: {read_count['viruses']}/{read_pair_count['viruses']}\n")


    def final_process():
        message = ""
        if len(primary_mapped) > 0:
            message += "    primary mapped dictionary is not empty\n"
        if len(discordant_mapped) > 0:
            message += "    discordant mapped dictionary is not empty\n"
        if len(unmapped) > 0:
            message += "    unmapped dictionary is not empty\n"
        if len(message) > 0:
            with pysam.AlignmentFile("%s.other.bam" % sample_id, "wb", header=in_sam.header) as other:
                for read_pair in primary_mapped.values():
                    for which_end in read_pair:
                        for read in read_pair[which_end]:
                            other.write(read)
                for read_pair in discordant_mapped.values():
                    for which_end in read_pair:
                        for read in read_pair[which_end]:
                            other.write(read)
                for read_pair in unmapped.values():
                    for which_end in read_pair:
                        for read in read_pair[which_end]:
                            other.write(read)
            sys.exit("error: \n" + message + "\n")
        else:
            sys.stdout.write("\ncomplete.\n")


    read_count = defaultdict(lambda: 0)
    read_pair_count = defaultdict(lambda:0)
    virus_count = defaultdict(lambda:defaultdict(lambda: 0))
    with (pysam.AlignmentFile("%s" % o.in_sam) as in_sam,
          pysam.AlignmentFile("%s.human.bam" % sample_id, "wb", header=in_sam.header) as human_bam,
          pysam.AlignmentFile("%s.viruses.bam" % sample_id, "wb", header=in_sam.header) as viruses_bam,
          pysam.AlignmentFile("%s.secondary.bam" % sample_id, "wb", header=in_sam.header) as secondary_bam,
          pysam.AlignmentFile("%s.supplementary.bam" % sample_id, "wb", header=in_sam.header) as supplementary_bam,
          pysam.AlignmentFile("%s.discordant.bam" % sample_id, "wb", header=in_sam.header) as discordant_bam,
          pysam.AlignmentFile("%s.qcfail.bam" % sample_id, "wb", header=in_sam.header) as qcfail_bam,
          gzip.open("%s.Unmapped.R1.fastq.gz" % sample_id, "wt") as unmapped_R1,
          gzip.open("%s.Unmapped.R2.fastq.gz" % sample_id, "wt") as unmapped_R2):
        for read in in_sam:
            read_count["total"] += 1

            if read.is_duplicate:
                sys.exit("Duplicate reads found")

            if read.is_qcfail:
                qcfail_bam.write(read)
                read_count["qcfail"] += 1
                continue

            #== check if the read is secondary mapped
            if read.is_secondary:
                # Secondary alignments occur when a read maps to multiple locations
                # in the genome with nearly the same quality.
                secondary_bam.write(read)
                read_count["secondary"] += 1
                continue
            #== check if the read is supplementary mapped.
            if read.is_supplementary:
                # Supplementary alignments are used to represent the alignment of a read that maps in multiple segments,
                # often due to structural variations such as large insertions, deletions, or inversions.
                # In such cases, a single read may be split into multiple parts,
                # each aligning to different regions of the genome.
                supplementary_bam.write(read)
                read_count["supplementary"] += 1
                continue

            ### count reads ###
            if read.is_mapped:
                read_count["mapped"] += 1
                contig_name = get_entrez_id(read.reference_name)
                if contig_name.startswith("chr"):
                    read_count["human"] += 1
                else:
                    read_count["viruses"] += 1
                    virus_count[contig_info[contig_name]['organism']]["read"] += 1
            else:
                read_count["unmapped"] += 1

            ### count read pairs ###
            #== orphans
            if not read.is_paired:
                sys.exit(f"Orphaned reads: {read.query_name}\n")

            #== identify which end of read pair
            if read.is_read1:
                which_end = "R1"
            elif read.is_read2:
                which_end = "R2"
            else:
                sys.exit(f"Unknown read end: {read.query_name}\n")

            #== check if both reads are mapped
            if read.is_mapped and read.mate_is_mapped:
                #== primary mapped reads
                # In the SAM/BAM format, primary, secondary, and supplementary alignments are
                # mutually exclusive categories used to describe the status of a read's alignment.

                #== write the ready read pairs
                if read.query_name not in primary_mapped:
                    write_primary_mapped_read_pairs()
                if read.query_name not in discordant_mapped:
                    write_discordant_read_pairs()

                #== check the read and add to one of the dictionaries
                contig_name = get_entrez_id(read.reference_name)
                next_contig_name = get_entrez_id(read.next_reference_name)
                if contig_name.startswith("chr") and next_contig_name.startswith("chr"):
                    #== for human chromosomes
                    # chr1
                    # chr1_KI270706v1_random
                    if contig_name.split('_')[0] == next_contig_name.split('_')[0]:
                        #== check if both reads mapped to the same chromosome
                        # The main alignment for the read, typically covering the entire read length.
                        primary_mapped[read.query_name][which_end].append(read)
                    # chrUn_GL000218v1
                    elif contig_name.startswith('chrUn') or next_contig_name.startswith('chrUn'):
                        # == chrUn is an exception
                        primary_mapped[read.query_name][which_end].append(read)
                    else:
                        # When read1 and read2 from a paired-end sequencing run are mapped to different chromosomes,
                        # this is referred to as discordant or anomalous mapping.
                        # These reads are often called discordant read pairs or inter-chromosomal read pairs.
                        discordant_mapped[read.query_name][which_end].append(read)
                elif (not contig_name.startswith("chr")) and (not next_contig_name.startswith("chr")):
                    #== for viruses
                    if contig_info[contig_name]['organism'] == contig_info[next_contig_name]['organism']:
                        #== determine if both reads mapped to REO1
                        primary_mapped[read.query_name][which_end].append(read)
                    else:
                        #== differnt viruses
                        discordant_mapped[read.query_name][which_end].append(read)
                else:
                    #== one is human and one is virus
                    discordant_mapped[read.query_name][which_end].append(read)
            #== either read is unmapped
            else:
                if read.query_name not in unmapped:
                    write_unmapped_read_pairs()
                unmapped[read.query_name][which_end].append(read)

            #== write to files every 1,000,000 reads
            if read_pair_count["total"] % 1000000 == 0:
                write_primary_mapped_read_pairs()
                write_discordant_read_pairs()
                write_unmapped_read_pairs()
                print_progress()

        # write the rest reads
        write_primary_mapped_read_pairs()
        write_discordant_read_pairs()
        write_unmapped_read_pairs()

        final_process()
        print_progress()

    with open("%s.read_count.txt" % sample_id, "w") as outfile:
        viruses = list(virus_count.keys())
        viruses.sort()

        # header1 = "\t".join(("#Sample",
        #                      "Total", "Total",
        #                      "QC Failed",
        #                      "Mapped", "Mapped", "Unmapped", "Unmapped",
        #                      "Secondary", "Supplementary",
        #                      "Primary", "Discordant",
        #                      "Human", "Human", "Viruses", "Viruses"))
        # header2 = "\t".join(("#Count Method",
        #                      "reads", "pairs",
        #                      "reads",
        #                      "reads", "pairs", "reads", "pairs",
        #                      "reads", "reads",
        #                      "pairs", "pairs",
        #                      "reads", "pairs", "reads", "pairs"))
        # header1 = header1 + "\t" + "\t".join([v for v in viruses for _ in range(2)])
        # header2 = header2 + "\t" + "\t".join(["reads", "pairs"] * len(viruses))
        # outfile.write(header1 + "\n" + header2 + "\n")

        # outfile.write(f"{sample_id}\t"
        #               f"{read_count['total']}\t{read_pair_count['total']}\t"
        #               f"{read_count['qcfail']}\t"
        #               f"{read_count['mapped']}\t{read_pair_count['mapped']}\t"
        #               f"{read_count['unmapped']}\t{read_pair_count['unmapped']}\t"
        #               f"{read_count['secondary']}\t{read_count['supplementary']}\t"
        #               f"{read_pair_count['primary']}\t{read_pair_count['discordant']}\t"
        #               f"{read_count['human']}\t{read_pair_count['human']}\t"
        #               f"{read_count['viruses']}\t{read_pair_count['viruses']}")
        # for virus in viruses:
        #     outfile.write(f"\t{virus_count[virus]['read']}\t{virus_count[virus]['pair']}")
        # outfile.write("\n")

        header = "\t".join((sample_id, "Count Method", "Count Number"))
        outfile.write(header + "\n")
        outfile.write(f"#Total\treads\t{read_count['total']}\n"
                      f"#Total\tpairs\t{read_pair_count['total']}\n"
                      f"#QC Failed\treads\t{read_count['qcfail']}\n"
                      f"#Mapped\treads\t{read_count['mapped']}\n"
                      f"#Mapped\tpairs\t{read_pair_count['mapped']}\n"
                      f"#Unmapped\treads\t{read_count['unmapped']}\n"
                      f"#Unmapped\tpairs\t{read_pair_count['unmapped']}\n"
                      f"#Secondary\treads\t{read_count['secondary']}\n"
                      f"#Supplementary\treads\t{read_count['supplementary']}\n"
                      f"#Primary\tpairs\t{read_pair_count['primary']}\n"
                      f"#Discordant\tpairs\t{read_pair_count['discordant']}\n"
                      f"#Human\treads\t{read_count['human']}\n"
                      f"#Human\tpairs\t{read_pair_count['human']}\n"
                      f"#Virus\treads\t{read_count['viruses']}\n"
                      f"#Virus\tpairs\t{read_pair_count['viruses']}\n")
        lines = ""
        for virus in viruses:
            lines += f"{virus}\treads\t{virus_count[virus]['read']}\n"
            lines += f"{virus}\tpairs\t{virus_count[virus]['pair']}\n"
        outfile.write(lines)

if __name__ == "__main__":
    main()
