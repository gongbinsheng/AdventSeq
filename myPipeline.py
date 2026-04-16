import re
import sys
import subprocess
import inspect
from collections import defaultdict, OrderedDict
from pathlib import Path

import yaml

def list_files_with_extensions(folder_path, extensions):
    file_list = []
    for item in Path(folder_path).iterdir():
        if item.is_file() and item.suffix in extensions:
            file_list.append(item.as_posix())
        elif item.is_dir():
            file_list.extend(list_files_with_extensions(item, extensions))
    return file_list

class OrderedDefaultDict(defaultdict):
    def __init__(self, default_factory=None, *args, **kwargs):
        if default_factory and not callable(default_factory):
            raise TypeError("first argument must be callable")
        super(OrderedDefaultDict, self).__init__(default_factory, *args, **kwargs)
        self._keys = []

    def __setitem__(self, key, value):
        if key not in self:
            self._keys.append(key)
        super(OrderedDefaultDict, self).__setitem__(key, value)

    def __delitem__(self, key):
        super(OrderedDefaultDict, self).__delitem__(key)
        self._keys.remove(key)

    def __iter__(self):
        return iter(self._keys)

    def keys(self):
        return list(self._keys)

    def values(self):
        return [self[key] for key in self._keys]

    def items(self):
        return [(key, self[key]) for key in self._keys]

def is_valid_DNA_sequence(s: object) -> bool:
    if not isinstance(s, str):
        sys.exit("Adapter is not a string.\n")
    valid_characters = set("ATCG")
    return all(char in valid_characters for char in s)

def parse_gtf(file_path):

    with open(file_path, 'r') as gtf_file:
        for line in gtf_file:
            if line.startswith('#'):
                continue  # Skip comment lines
            fields = line.strip().split('\t')
            # parse line
            try:
                chrom, source, feature_type, start1, end1, score, strand, frame, attributes = fields[:9]
            except Exception as e:
                sys.exit(f"Could not parse line:\n{line}\n\n{e}\n")
            # parse attributes
            try:
                attributes_dict = {}
                attributes_list = attributes.strip().split(';')
                for pair in attributes_list:
                    if pair.strip():  # Check if pair is not empty
                        key, value = pair.strip().split(' ', 1)  # Split by first space to separate key and value
                        attributes_dict[key] = value.strip('"')
            except Exception as e:
                sys.exit(f"Could not parse attributes:\n{line}\n\n{e}\n")
            # save information into a dictionary
            gtf_record = {
                'chrom': chrom,
                'source': source,
                'feature_type': feature_type,
                'start1': int(start1),
                'end1': int(end1),
                'score': score,
                'strand': strand,
                'frame': frame,
                'attributes': attributes_dict
            }
            yield gtf_record

DEFAULT_ENV_CONFIG = Path(__file__).with_name("pipeline_settings.yml")


def _normalize_path(path):
    return Path(path).expanduser().as_posix()

class Pipeline:
    __HPC_nodes_to_skipped = set()
    HPC_nodes_to_use = set()
    threadN = 4
    memory = None
    WD = None
    Data_folder = None
    adapter = None
    adapter_fasta = None

    config_path = DEFAULT_ENV_CONFIG
    config = {}
    envs4steps = {}
    configured_reference_paths = {}

    @classmethod
    def _resolve_config_path(cls, config_path=None):
        return Path(config_path or cls.config_path)


    @classmethod
    def _extract_envs4steps(cls, config, config_path):
        step_envs = config.get("envs4steps", config)
        if not isinstance(step_envs, dict):
            sys.exit(
                f"envs4steps must be a YAML mapping.\n"
                f"    {config_path}\n"
            )

        envs4steps = {}
        for step_name, env_name in step_envs.items():
            if not isinstance(step_name, str) or not isinstance(env_name, str):
                sys.exit(
                    f"Each envs4steps entry must map a step name to a conda env string:\n"
                    f"    {config_path}\n"
                )
            envs4steps[step_name] = env_name
        return envs4steps


    @classmethod
    def _extract_reference_paths(cls, config):
        reference_paths = config.get("reference_paths", {})
        if reference_paths is None:
            return {}
        if not isinstance(reference_paths, dict):
            sys.exit("reference_paths must be a YAML mapping.\n")

        bundled_reference_paths = {}
        for bundle_name, bundle_value in reference_paths.items():
            if not isinstance(bundle_name, str):
                sys.exit("Each reference_paths key must be a string.\n")
            if not isinstance(bundle_value, dict):
                sys.exit(
                    "Each reference_paths entry must map a ref_genome to a YAML mapping of paths.\n"
                )
            normalized_bundle = {}
            for path_name, path_value in bundle_value.items():
                if not isinstance(path_name, str) or not isinstance(path_value, str):
                    sys.exit(
                        "Each reference_paths bundle entry must map a path name to a string path.\n"
                    )
                normalized_bundle[path_name] = path_value
            bundled_reference_paths[bundle_name] = normalized_bundle

        return bundled_reference_paths


    @classmethod
    def load_pipeline_config(cls, config_path=None):
        env_config_path = cls._resolve_config_path(config_path)

        if not env_config_path.exists():
            example_path = env_config_path.with_suffix(env_config_path.suffix + ".example")
            sys.exit(
                f"Cannot find conda env config file:\n"
                f"  {env_config_path}\n"
                f"Create it from the example file:\n"
                f"  {example_path}\n"
            )

        try:
            with open(env_config_path, "r", encoding="utf-8") as handle:
                config = yaml.safe_load(handle) or {}
        except yaml.YAMLError as e:
            sys.exit(f"Could not parse YAML file:\n  {env_config_path}\n{e}\n")
        except OSError as e:
            sys.exit(f"Could not read env config file:\n  {env_config_path}\n{e}\n")

        if not isinstance(config, dict):
            sys.exit(f"Env config must be a YAML mapping:\n  {env_config_path}\n")

        cls.config_path = env_config_path
        cls.config = config
        cls.envs4steps = cls._extract_envs4steps(config, env_config_path)
        cls.configured_reference_paths = cls._extract_reference_paths(config)
        return config


    @classmethod
    def ensure_pipeline_config_loaded(cls):
        if cls.config:
            return cls.config
        return cls.load_pipeline_config()

    def add_nodes_to_be_skipped(self,more_nodes_to_be_skipped):
        self.__HPC_nodes_to_skipped = self.__HPC_nodes_to_skipped.union(more_nodes_to_be_skipped)


    def __init__(self, sample_id, file_list, log_file, ref_genome, ref_type, feature, library_source, library_layout, platform, filetype):
        self.ensure_pipeline_config_loaded()
        if self.WD is None:
            sys.stderr.write("Working directory (WD) is not set. Set it to current directory.")
            self.WD = Path.cwd()
        else:
            self.WD = Path(self.WD)
        if self.Data_folder is None:
            sys.exit("Data folder is not given.")

        self.batch = OrderedDefaultDict(lambda:[])
        self.__current_step_id = ""
        self.__reset = 0

        self.sample_id = sample_id
        self.file_list = file_list
        self.log_file = log_file

        self.ref_genome = ref_genome
        self.ref_type = ref_type # genome, transcriptome
        self.feature = feature # transcript, exon, CDS
        self.library_source = library_source.upper() # GENOMIC, TRANSCRIPTOMIC
        self.library_layout = self.__get_library_layout(library_layout) # paired, single
        self.check_adapter()
        self.platform = platform.upper() # ILLUMINA, PACBIO_SMRT, OXFORD_NANOPORE
        self.filetype = filetype.lower() # fastq, bam
        self.genome_fasta = None
        self.bwa_index = None
        self.bowtie2_index = None
        self.minimap2_index = None
        self.hisat2_index = None
        self.star_index = None
        self.gtf = None
        self.__load_reference_paths_from_config()

        self.required_conda_envs = set()
        self.read_progress()
        self.is_completed = True
        try:
            self.raw_data_type = self.__get_data_type()
        except ValueError as e:
            sys.exit(f"Error: {e}")
        if not self.raw_data_type.lower().startswith(self.filetype.lower()):
            sys.exit("File type in SRA metadata does not match the actual file type in data folder.\n")
        self.mycat = {"FASTQ.GZ":"zcat",
                      "FASTQ":"cat",
                      "BAM":"samtools view"}

        self.FASTQ_filename_pattern = re.compile(fr'{self.sample_id}[_.-]*(?P<part1>.*?)[_.-]+(r|read)[_.-]*[0-9]+[_.-]*(?P<part2>.*?)\.(fastq|fq).*', re.IGNORECASE)
        self.fastq_data = self.__get_fastq_data()
        self.__FASTQ_ENV = []

        self.RG_is_added = False
        self.host_read_pairs_removed = False
        self.FASTQ_ready_for_contigs = set()


    def __get_data_type(self):
        extensions = set()
        gzext = ""
        for file in self.file_list:
            suffixes = Path(file).suffixes
            if suffixes[-1].lower() == ".gz":
                gzext = ".GZ"
                suffixes = suffixes[:-1]
            extensions.add(suffixes[-1].lower()) # Extract the file extension

        if len(extensions) == 1:
            ext = extensions.pop()
            if ext in [".fastq", ".fq"]:
                return 'FASTQ' + gzext
            elif ext == '.bam':
                return 'BAM'
            else:
                raise ValueError('Unsupported file type: %s\n' % ext)
        else:
            raise ValueError('Multiple file types detected: %s\n' % ", ".join(list(extensions)))


    def __get_library_layout(self, library_layout):
        if library_layout is None:
            return None
        else:
            if library_layout.lower() not in ("single","paired"):
                sys.exit("Unknown library layout: %s" % library_layout)
            else:
                return library_layout.lower()


    def __get_fastq_data(self):
        if self.filetype.upper() != "FASTQ":
            return None
        lanes = set()
        if self.library_layout.lower() == "single":
            fastq_data = OrderedDict()
            for file in self.file_list:
                parts = self.FASTQ_filename_pattern.match(Path(file).name)
                if parts is not None:
                    part = "_".join((parts.group('part1').strip("[_.-]"),
                                     parts.group('part2').strip("[_.-]"))).strip("_")
                else:
                    sys.exit(f"File name:\n\t{file}\ndoes not match the pattern:\n\t{self.FASTQ_filename_pattern.pattern}\n")
                if part not in lanes:
                    lanes.add(part)
                    fastq_data[part] = file
                else:
                    sys.exit(f"Duplicated lane name [ {part} ] for sample [ {self.sample_id} ]\n")
        elif self.library_layout.lower() == "paired":
            fastq_data = defaultdict(lambda:{"R1":None, "R2":None})
            if len(self.file_list) % 2 != 0:
                sys.exit("Odd number of files found for paired-end data.\n")
            i = 0
            for file in self.file_list:
                end = i % 2 + 1
                parts = self.FASTQ_filename_pattern.match(Path(file).name)
                if parts is not None:
                    part = "_".join((parts.group('part1').strip("[_.-]"),
                                     parts.group('part2').strip("[_.-]"))).strip("_")
                else:
                    sys.exit(f"File name:\n\t{file}\ndoes not match the pattern:\n\t{self.FASTQ_filename_pattern.pattern}\n")
                if fastq_data[part]["R%d" % end] is None:
                    fastq_data[part]["R%d" % end] = (file)
                    i += 1
                else:
                    sys.exit(f"Duplicated lane name [ {part} ] for sample [ {self.sample_id} ]\n")
            for part in fastq_data:
                if fastq_data[part]["R1"] is None or fastq_data[part]["R2"] is None:
                    sys.exit(f"The part/lane name does not match between R1 and R2 for sample: {self.sample_id}\n")
            return fastq_data


    def check_conda_envs(self):
        missing_envs = []
        try:
            # Get a list of all Conda environments
            available_conda_envs = subprocess.check_output(["conda", "env", "list"]).decode("utf-8")
            # Iterate through the provided environment names
            for env in self.required_conda_envs:
                if env not in available_conda_envs:
                    missing_envs.append(env)
            return missing_envs
        except Exception as e:
            sys.exit(f"An error occurred: {e}")


    def __validate_existing_file(self, path, label):
        normalized_path = _normalize_path(path)
        if not Path(normalized_path).is_file():
            sys.exit(f"Invalid {label}: file not found.\n  {normalized_path}\n")
        return normalized_path


    def __validate_existing_directory(self, path, label):
        normalized_path = _normalize_path(path)
        if not Path(normalized_path).is_dir():
            sys.exit(f"Invalid {label}: directory not found.\n  {normalized_path}\n")
        return normalized_path


    def __validate_index_prefix(self, path, label, suffix_groups):
        normalized_path = _normalize_path(path)
        for suffixes in suffix_groups:
            if all(Path(f"{normalized_path}{suffix}").is_file() for suffix in suffixes):
                return normalized_path
        sys.exit(f"Invalid {label}: index files not found for prefix.\n  {normalized_path}\n")


    def __set_reference_path(self, attr_name, path, validator):
        normalized_path = validator(path)
        setattr(self, attr_name, normalized_path)


    def __get_required_reference_path(self, attr_name, label, setter_name, arg_name):
        path = getattr(self, attr_name)
        if path is None:
            sys.exit(
                f"{label} is not set.\n"
                f"Use {setter_name}() before running this step, or pass {arg_name}=...\n"
            )
        return path


    def __get_ref_genome(self):
        if self.ref_genome is None:
            sys.exit(
                "ref_genome is not set.\n"
                "Use set_ref_genome(), or pass ref_genome=...\n"
            )
        return self.ref_genome


    def __load_reference_paths_from_config(self):
        setter_by_attr = {
            "genome_fasta": self.set_genome_fasta,
            "bwa_index": self.set_bwa_index,
            "bowtie2_index": self.set_bowtie2_index,
            "minimap2_index": self.set_minimap2_index,
            "hisat2_index": self.set_hisat2_index,
            "star_index": self.set_star_index,
            "gtf": self.set_gtf,
        }

        if self.ref_genome is None:
            return

        selected_reference_paths = self.configured_reference_paths.get(self.ref_genome)
        if selected_reference_paths is None:
            return

        for attr_name, path in selected_reference_paths.items():
            setter = setter_by_attr.get(attr_name)
            if setter is None:
                valid_keys = ", ".join(sorted(setter_by_attr))
                sys.exit(
                    f"Unknown reference_paths key in YAML for {self.ref_genome}: {attr_name}\n"
                    f"Valid keys are: {valid_keys}\n"
                )
            setter(path)


    def __clear_reference_paths(self):
        self.genome_fasta = None
        self.bwa_index = None
        self.bowtie2_index = None
        self.minimap2_index = None
        self.hisat2_index = None
        self.star_index = None
        self.gtf = None


    def set_ref_genome(self, ref_genome):
        if not isinstance(ref_genome, str) or not ref_genome.strip():
            sys.exit("ref_genome must be a non-empty string.\n")
        if self.ref_genome == ref_genome:
            return
        self.ref_genome = ref_genome
        self.__clear_reference_paths()
        self.__load_reference_paths_from_config()

    def set_genome_fasta(self, path):
        self.__set_reference_path(
            "genome_fasta",
            path,
            lambda value: self.__validate_existing_file(value, "genome FASTA"),
        )


    def set_bwa_index(self, path):
        self.__set_reference_path(
            "bwa_index",
            path,
            lambda value: self.__validate_index_prefix(
                value,
                "BWA index",
                [(".amb", ".ann", ".bwt", ".pac", ".sa")],
            ),
        )


    def set_bowtie2_index(self, path):
        self.__set_reference_path(
            "bowtie2_index",
            path,
            lambda value: self.__validate_index_prefix(
                value,
                "Bowtie2 index",
                [
                    (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
                    (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"),
                ],
            ),
        )


    def set_minimap2_index(self, path):
        self.__set_reference_path(
            "minimap2_index",
            path,
            lambda value: self.__validate_existing_file(value, "minimap2 index"),
        )


    def set_hisat2_index(self, path):
        self.__set_reference_path(
            "hisat2_index",
            path,
            lambda value: self.__validate_index_prefix(
                value,
                "HISAT2 index",
                [
                    tuple(f".{i}.ht2" for i in range(1, 9)),
                    tuple(f".{i}.ht2l" for i in range(1, 9)),
                ],
            ),
        )


    def set_star_index(self, path):
        self.__set_reference_path(
            "star_index",
            path,
            lambda value: self.__validate_existing_directory(value, "STAR index"),
        )


    def set_gtf(self, path):
        self.__set_reference_path(
            "gtf",
            path,
            lambda value: self.__validate_existing_file(value, "GTF"),
        )


    def __get_genome_fasta(self):
        return self.__get_required_reference_path(
            "genome_fasta",
            "Genome FASTA",
            "set_genome_fasta",
            "genome_fasta",
        )


    def __get_bwa_index(self):
        return self.__get_required_reference_path(
            "bwa_index",
            "BWA index",
            "set_bwa_index",
            "bwa_index",
        )


    def __get_bowtie2_index(self):
        return self.__get_required_reference_path(
            "bowtie2_index",
            "Bowtie2 index",
            "set_bowtie2_index",
            "bowtie2_index",
        )


    def __get_minimap2_index(self):
        return self.__get_required_reference_path(
            "minimap2_index",
            "minimap2 index",
            "set_minimap2_index",
            "minimap2_index",
        )


    def __get_hisat2_index(self):
        return self.__get_required_reference_path(
            "hisat2_index",
            "HISAT2 index",
            "set_hisat2_index",
            "hisat2_index",
        )


    def __get_star_index(self):
        return self.__get_required_reference_path(
            "star_index",
            "STAR index",
            "set_star_index",
            "star_index",
        )


    def __get_gtf(self):
        return self.__get_required_reference_path(
            "gtf",
            "GTF",
            "set_gtf",
            "gtf",
        )


    def check_adapter(self):
        if self.adapter is not None:
            if self.library_layout == "single":
                if not is_valid_DNA_sequence(self.adapter):
                    sys.exit("The adapter is not a valid DNA sequence.")
            elif self.library_layout == "paired":
                if not isinstance(self.adapter, dict):
                    sys.exit("The PE adapter is not a dictionary.")
                if set(self.adapter) != {"R1", "R2"}:
                    sys.exit('The PE adapter should only contains "R1" and "R2" as dictionary keys.')
                for key in self.adapter:
                    if not is_valid_DNA_sequence(self.adapter[key]):
                        sys.exit("The adapter[%s] is not a valid DNA sequence." % key)
        if self.adapter_fasta is not None:
            valid_nucleotides = set('ATCG')
            if not Path(self.adapter_fasta).is_file():
                sys.exit("Can't find the adapter fasta file: \n  %s\n" % self.adapter_fasta)
            with open(self.adapter_fasta, 'r') as file:
                is_header = False
                seq_id = ""
                for line in file:
                    line = line.strip()
                    if not line:
                        continue  # Skip empty lines
                    if line.startswith('>'):
                        is_header = True
                        seq_id = line
                        continue
                    elif is_header:
                        is_header = False
                        if not all(base in valid_nucleotides for base in line.upper()):
                            sys.exit(f"Invalid sequence line in adapter fasta file:\n  {self.adapter_fasta}\n    {seq_id}\n    {line}\n")
                    else:
                        sys.exit(f"Invalid format in adapter fasta file:\n  {self.adapter_fasta}\n    {line}")


    def set_current_step_id(self, step_id):
        self.__current_step_id = step_id


    def get_current_step_id(self):
        return self.__current_step_id


    def __set_FASTQs(self, *args):
        if len(args) == 1:
            self.__FASTQ_ENV = [args[0]]
        elif len(args) == 2:
            self.__FASTQ_ENV = args
        else:
            sys.exit("Error: __set_FASTQs() accepts either one (single-end) or two (paired-end) arguments.\n")


    def reset_pipeline(self, **kwargs):
        batch = []
        self.__reset += 1
        #fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.sample_id + "|" + f"Round{self.__reset}"
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        # reset SID
        batch.append('SID="%s"\n' % self.sample_id)
        # reset ref_genome
        if "ref_genome" in kwargs:
            self.set_ref_genome(kwargs["ref_genome"])
            batch.append('ref_genome="%s"\n' % kwargs["ref_genome"])
        # reset FASTQ or FASTQ_R1 and FASTQ_R2
        if len(self.__FASTQ_ENV) == 1:
            batch.append('FASTQ="%s"\n' % self.__FASTQ_ENV[0])
        elif len(self.__FASTQ_ENV) == 2:
            batch.append('FASTQ_R1="%s"\n' % self.__FASTQ_ENV[0])
            batch.append('FASTQ_R2="%s"\n' % self.__FASTQ_ENV[1])
        else:
            sys.exit("Error: FASTQ_ENV should contain either one (single-end) or two (paired-end) elements.\n")
        self.set_current_step_id(self.sample_id)  # reset current step id to sample_id
        self.batch[step_id].extend(batch)

    def set_qsub_parameters(self):
        batch = []
        batch.append('#!/bin/bash\n')
        batch.append('#\n')
        batch.append('#$ -N "BG_%s"\n' % self.sample_id)
        batch.append('#$ -S /bin/bash\n')
        batch.append('#$ -cwd\n')
        batch.append('#$ -j y\n')
        if len(self.HPC_nodes_to_use) > 0:
            batch.append('#$ -l h=(%s)\n' % "|".join(self.HPC_nodes_to_use))
        else:
            if len(self.__HPC_nodes_to_skipped) > 0:
                batch.append('#$ -l h=!(%s)\n' % "|".join(self.__HPC_nodes_to_skipped))
        if self.memory is not None:
            batch.append('#$ -l h_vmem=%s\n' % self.memory)
        batch.append('#$ -pe smp %s\n' % self.threadN)
        batch.append('#$ -R y\n')
        batch.append('#$ -o "%s"\n\n' % self.log_file)
        self.batch["HPC"].extend(batch)


    def set_slurm_parameters(self):
        batch = []
        batch.append('#!/bin/bash\n')
        batch.append('#\n')
        batch.append('#SBATCH --job-name="BG_%s"\n' % self.sample_id)
        batch.append('#SBATCH --output="%s"\n' % self.log_file)
        batch.append('#SBATCH --error="%s"\n' % self.log_file)
        if self.memory is not None:
            batch.append('##SBATCH --mem=%s\n' % self.memory)
        batch.append('#SBATCH --cpus-per-task=%s\n' % self.threadN)
        batch.append('\n')
        self.batch["HPC"].extend(batch)


    def set_env_variables(self):
        batch = []
        step_id = "ENV"
        batch.append('### %s ###\n' % "Set ENV variables")
        batch.append('# Conda initiation\n')
        batch.append('source /account001/bgong/init_conda.sh\n\n\n')  # Conda initiation
        batch.append('\n')
        batch.append('SID="%s"\n' % self.sample_id)
        self.set_current_step_id(self.sample_id) # init the current step id
        batch.append('WD="%s/$SID"\n' % self.WD)
        batch.append('[ -d "$WD" ] || mkdir -p "$WD"\n')
        batch.append('cd "$WD"\n')
        batch.append('my_progress="_my_progress"\n')
        batch.append('\n\n')
        self.batch[step_id].extend(batch)


    def __update_env_variable(self, env_var, value):
        for i in range(len(self.batch["ENV"])):
            if self.batch["ENV"][i].startswith(env_var+"="):
                self.batch["ENV"][i] = '%s="%s"\n' % (env_var, value)
                return


    def read_progress(self):
        self.progress = {}
        my_progress = self.WD / self.sample_id / "_my_progress"
        if my_progress.is_file():
            with open(my_progress, 'r') as f:
                for line in f:
                    strs = line.strip().split("\t")
                    if len(strs) == 2:
                        self.progress[strs[0]] = strs[1]
            # write progress back
            with open(my_progress, 'w') as f:
                for step, timestamp in self.progress.items():
                    f.write(f"{step}\t{timestamp}\n")


    def __comment_lines(self, lines):
        commented_lines = []
        for line in lines:
            if (
                    line[0] == "#"
                    or line.startswith("mapper=")
                    or line.startswith("BAM=")
                    or line.startswith("SAM=")
                    or line.startswith("sorted_BAM=")
                    or line.startswith("FASTQ")
                    or line.startswith("VCF")
                    or line.startswith("BED")
                    or line.startswith("ref_genome=")
                    or line.startswith("SID=")
                    or line.startswith("OLD")
            ):
                commented_lines.append(line)
            else:
                commented_lines.append("# " + line)
        return commented_lines


    def MultiQC(self, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        if self.raw_data_type.upper() == "BAM":
            batch.append('### %s was skipped ###\n\n' % step_id)
            batch.append("# %s can only run with FASTQ files.\n#\t%s\n" % (step_id, "\n#\t".join(self.file_list)))
            return
        # set up env variables
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        batch.append('[ -d "%s" ] && rm -fr "%s"\n' % (fn_name, fn_name))
        batch.append('mkdir "%s"\n' % fn_name)
        # FastQC
        batch.append('fastqc \\\n')
        batch.append('    --threads %d -o "%s" \\\n' % (self.threadN, fn_name))
        for file in self.file_list:
            batch.append('    "%s" \\\n' % file)
        batch.append('    2>&1\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - fastqc failed."; exit 1; fi\n' % step_id)
        # MultiQC
        batch.append('multiqc --force --outdir "%s" --filename "${SID}.MultiQC.html" "%s"\n' % (fn_name, fn_name))
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - multiqc failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def merge_lanes(self, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        #batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        #self.required_conda_envs.add(self.envs4steps[fn_name]) # add env name to required env list
        # commands
        if self.library_layout.lower() == "single":
            if len(self.fastq_data) == 1:
                lane = next(iter(self.fastq_data.values()))
                batch.append('FASTQ="%s"\n' % lane)
                self.__set_FASTQs(lane)
            else:
                batch.append('FASTQ="${SID}.merge_lanes.fastq.gz"\n')
                self.__set_FASTQs("${SID}.merge_lanes.fastq.gz")
                batch.append('%s <' % self.mycat[self.raw_data_type])
                for part in self.fastq_data:
                    batch.append(' "%s"' % self.fastq_data[part])
                batch.append(' > "$FASTQ"\n\n')
        elif self.library_layout.lower() == "paired":
            if len(self.fastq_data) == 1:
                lane = next(iter(self.fastq_data.values()))
                batch.append('FASTQ_R1="%s"\n' % lane['R1'])
                batch.append('FASTQ_R2="%s"\n' % lane['R2'])
                self.__set_FASTQs(lane['R1'], lane['R2'])
            else:
                batch.append('FASTQ_R1="${SID}.merge_lanes.R1.fastq.gz"\n')
                batch.append('FASTQ_R2="${SID}.merge_lanes.R2.fastq.gz"\n')
                self.__set_FASTQs("${SID}.merge_lanes.R1.fastq.gz", "${SID}.merge_lanes.R2.fastq.gz")
                # R1
                batch.append('%s <' % self.mycat[self.raw_data_type])
                for part in self.fastq_data:
                    batch.append(' "%s"' % self.fastq_data[part]["R1"])
                batch.append(' | gzip -c > "$FASTQ_R1"\n')
                batch.append('if [ $? -ne 0 ]; then echo "Error: %s - gzip failed."; exit 1; fi\n\n' % step_id)
                # R2
                batch.append('%s <' % self.mycat[self.raw_data_type])
                for part in self.fastq_data:
                    batch.append(' "%s"' % self.fastq_data[part]["R2"])
                batch.append(' | gzip -c > "$FASTQ_R2"\n')
                batch.append('if [ $? -ne 0 ]; then echo "Error: %s - gzip failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        #batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is just merging lanes, no change to data, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def fastp(self, dedup=False, trim_adapter=True, trim_polyG=True, trim_polyX=True, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name]) # add env name to required env list
        # commands
        tags = []
        if dedup:
            tags.append("dedup")
        if trim_adapter:
            tags.append("trim_adapter")
        tag = "_".join(tags)
        batch.append('[ -d "%s" ] && rm -fr "%s"\n' % (fn_name, fn_name))
        batch.append('mkdir "%s"\n\n' % fn_name)

        batch.append('fastp \\\n')
        if self.library_layout.lower() == "single":
            batch.insert(-1, 'oldFASTQ="$FASTQ"\n')
            batch.insert(-1, 'FASTQ="${SID}.fastp_%s.fastq.gz"\n' % tag)
            self.__set_FASTQs("${SID}.fastp_%s.fastq.gz" % tag)
            batch.append('    -i="$oldFASTQ" \\\n')
            batch.append('    -o="$FASTQ" \\\n')
            if trim_adapter:
                if self.adapter_fasta is not None:
                    batch.append('    --adapter_fasta="%s" \\\n' % self.adapter_fasta)
                elif self.adapter is not None:
                    self.batch["ENV"].insert(-1, 'ADAPTER=%s\n' % self.adapter)
                    batch.append('    --adapter_sequence=$ADAPTER \\\n')
            else:
                batch.append('    --disable_adapter_trimming \\\n')
        elif self.library_layout.lower() == "paired":
            batch.insert(-1, 'oldFASTQ_R1="$FASTQ_R1"\n')
            batch.insert(-1, 'oldFASTQ_R2="$FASTQ_R2"\n')
            batch.insert(-1, 'FASTQ_R1="${SID}.fastp_%s.R1.fastq.gz"\n' % tag)
            batch.insert(-1, 'FASTQ_R2="${SID}.fastp_%s.R2.fastq.gz"\n' % tag)
            self.__set_FASTQs("${SID}.fastp_%s.R1.fastq.gz" % tag, "${SID}.fastp_%s.R2.fastq.gz" % tag)
            batch.append('    --in1="$oldFASTQ_R1" \\\n')
            batch.append('    --in2="$oldFASTQ_R2" \\\n')
            batch.append('    --out1="$FASTQ_R1" \\\n')
            batch.append('    --out2="$FASTQ_R2" \\\n')
            if trim_adapter:
                if self.adapter_fasta is not None:
                    batch.append('    --adapter_fasta="%s" \\\n' % self.adapter_fasta)
                elif self.adapter is not None:
                    self.batch["ENV"].insert(-1, 'ADAPTER_R1=%s\n' % self.adapter["R1"])
                    self.batch["ENV"].insert(-1, 'ADAPTER_R2=%s\n' % self.adapter["R2"])
                    batch.append('    --adapter_sequence=$ADAPTER_R1 \\\n')
                    batch.append('    --adapter_sequence_r2=$ADAPTER_R2 \\\n')
                else:
                    self.batch["ENV"].insert(-1, "# Adapter will be auto detected by fastp.\n")
                    batch.append('    --detect_adapter_for_pe \\\n')
            else:
                batch.append('    --disable_adapter_trimming \\\n')
        #batch.append('    --disable_quality_filtering \\\n')
        #batch.append('    --disable_length_filtering \\\n')
        if trim_polyG:
            batch.append('    --trim_poly_g \\\n')
        if trim_polyX:
            batch.append('    --trim_poly_x \\\n')
        if dedup:
            batch.append('    --dedup \\\n')
        batch.append('    --json="%s/${SID}.fastp.json" \\\n' % fn_name)
        batch.append('    --html="%s/${SID}.fastp.html" \\\n' % fn_name)
        batch.append('    --thread=%d\n\n' % self.threadN)
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # for the current practice, fastp should be always run on raw data, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def sort_BAM(self, force=False, by_qname=False):
        batch = []
        self.sorted_by_qname = by_qname
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        if by_qname:
            step_id += "_by_qname"
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('sorted_BAM="${BAM%%.*.*}.sorted%s.bam"\n' % ("_by_qname" if by_qname else ""))
        batch.append('if [[ "sorted_BAM" == "$BAM" ]]; then\n')
        batch.append('  echo "BAM file was already sorted in the same way."\n')
        batch.append('else\n')
        batch.append('  samtools sort%s -@ %d -o "${sorted_BAM}" "${BAM}"\n' % (" -n" if by_qname else "", self.threadN))
        batch.append('fi\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - samtools sort failed."; exit 1; fi\n' % step_id)
        if not by_qname:
            batch.append('samtools index --bai "${sorted_BAM}"\n')
            batch.append('if [ $? -ne 0 ]; then echo "Error: %s - samtools index failed."; exit 1; fi\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # sort the BAM file, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def SAM2BAM(self, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('BAM="${SAM%.*}.unsorted.bam"\n')
        batch.append('samtools view -b -o "${BAM}" "${SAM}"\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - SAM to BAM failed."; exit 1; fi\n\n' % step_id)
        batch.append('rm -f "${SAM}"\n')
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # convert SAM to BAM, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def BAM_stat(self, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('echo "# samtools view" > "${BAM}.stats"\n')
        batch.append('echo -n -e "primary mapped pairs\\t" >> "${BAM}.stats"\n')
        batch.append('samtools view -c -f 65 -F 2316  -@ %d "$BAM" >> "${BAM}.stats"\n' % self.threadN)
        # -f 65 = include 1 (paired) + 64 (read1)
        # -F 2316 = remove 256 (secondary) + 2048 (supplementary) + 4 (read-unmapped) + 8 (mate-unmapped)
        batch.append('echo -n -e "primary mapped properly paired\\t" >> "${BAM}.stats"\n')
        batch.append('samtools view -c -f 67 -F 2304  -@ %d "$BAM" >> "${BAM}.stats"\n' % self.threadN)
        # -f 65 = 1 (paired) + 2 (properly paired) + 64 (read1)
        # -F 2304 = remove 256 (secondary) + 2048 (supplementary)
        # no need to remove 4 + 8, because 2 already guarantees both reads mapped
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - samtools view failed."; exit 1; fi\n\n' % step_id)

        batch.append(f'echo -e "\\n{"="*20}\\n" >> "${{BAM}}.stats"\n')
        batch.append('echo "# samtools flagstat" >> "${BAM}.stats"\n')
        batch.append('samtools flagstat -@ %d -O tsv "$BAM" >> "${BAM}.stats"\n' % self.threadN)
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - samtools stats failed."; exit 1; fi\n\n' % step_id)

        batch.append(f'echo -e "\\n{"="*20}\\n" >> "${{BAM}}.stats"\n')
        batch.append('echo "# samtools stats" >> "${BAM}.stats"\n')
        batch.append('samtools stats -@ %d "$BAM" >> "${BAM}.stats"\n' % self.threadN)
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - samtools stats failed."; exit 1; fi\n\n' % step_id)

        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is just generating stats, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def BAM_not_in_BED(self, region, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + "|" + region
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name]) # add env name to required env list
        # BAM
        if region == "exon":
            batch.append('BED="$EXOME_BED"\n')
        elif region == "genome":
            batch.append('BED="$GENOME_BED"\n')
        batch.append('BAM_not_in_BED="${BAM%%.bam}.not_in_%s.bam"\n' % region)
        batch.append(' bedtools intersect \\\n')
        batch.append('    -v \\\n') # Only report those entries in A that have no overlap in B.
        batch.append('    -abam "$BAM" \\\n')
        batch.append('    -b "$BED" \\\n')
        batch.append(' > "$BAM_not_in_BED"\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - bedtools intersect failed."; exit 1; fi\n\n' % step_id)
        # BED
        batch.append('bedtools genomecov \\\n')
        batch.append('    -ibam "$BAM_not_in_BED" \\\n')
        batch.append('    -bg \\\n')
        batch.append('    -pc \\\n')
        batch.append(' > "${BAM_not_in_BED}.bed"\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - bedtools genomecov failed."; exit 1; fi\n\n' % step_id)
        # coverage filter
        for cov in (100, 1000, 5000, 10000):
            batch.append('cat "${BAM_not_in_BED}.bed"')
            batch.append(" | awk '$4 > %d'" % (cov - 1,))
            batch.append(' | bedtools sort')
            batch.append(' | bedtools merge')
            batch.append(' > "${BAM_not_in_BED}.cov%d.bed"\n' % cov)
            batch.append('if [ $? -ne 0 ]; then echo "Error: %s - cov%d failed."; exit 1; fi\n\n' % (step_id,cov))
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def BWA_MEM(self, ref_genome=None, bwa_index=None, SAM=None, force=False):
        batch = []
        self.mapper = "bwa_mem"
        if bwa_index is not None:
            if ref_genome is None:
                sys.exit(f"Please set ref_genome for specified BWAIndex:\n    {bwa_index}\n")
            effective_ref_genome = ref_genome
            effective_bwa_index = self.__validate_index_prefix(
                bwa_index,
                "BWA index",
                [(".amb", ".ann", ".bwt", ".pac", ".sa")],
            )
        else:
            if ref_genome is not None:
                self.set_ref_genome(ref_genome)
            effective_ref_genome = self.__get_ref_genome()
            effective_bwa_index = self.__get_bwa_index()

        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + "_" + effective_ref_genome
        self.set_current_step_id(step_id)
        if self.raw_data_type == "BAM":
            batch.append('### %s was skipped ###\n\n' % step_id)
            batch.append("# %s can only run with FASTQ files.\n#\t%s\n" % (step_id, "\n#\t".join(self.file_list)))
            return
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('mapper="%s"\n' % self.mapper)
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands (timing, index, SAM)
        batch.append('ref_genome="%s"\n' % effective_ref_genome)
        batch.append('BWAIndex="%s"\n' % effective_bwa_index)
        if SAM is None:
            batch.append('SAM="${SID}.${mapper}.${ref_genome}.sam"\n')
        else:
            batch.append('SAM="%s"\n' % SAM)
        batch.append("start=$(date +%s)\n")  # timer
        batch.append('bwa mem \\\n')
        batch.append('    -R "@RG\\tID:$SID\\tSM:$SID\\tLB:$SID\\tPL:%s" \\\n' % (self.platform,))
        batch.append('    -t %d \\\n' % self.threadN)
        batch.append('    "$BWAIndex" \\\n')
        if self.library_layout.lower() == "single":
            batch.append('    "$FASTQ" \\\n')
        elif self.library_layout.lower() == "paired":
            batch.append('    "$FASTQ_R1" \\\n')
            batch.append('    "$FASTQ_R2" \\\n')
        batch.append(' > "${SAM}"\n\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append("end=$(date +%s)\n")  # timer
        batch.append("runtime=$((end - start))\n")  # timer
        batch.append(
            f'echo -e "$SID\\t$(hostname)\\t$mapper\\t$ref_genome\\t{self.threadN}\\t$runtime\\tseconds\\t$(date +\'%Y-%m-%d %H:%M:%S\')" >> ../runtime.txt\n')  # timer
        self.RG_is_added = True
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        self.set_current_step_id(step_id) # mapping is a critical step, set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def Bowtie2(self, ref_genome=None, bowtie2_index=None, SAM=None, end2end=True, force=False):
        batch = []
        self.mapper = "bowtie2"
        if bowtie2_index is not None:
            if ref_genome is None:
                sys.exit(f"Please set ref_genome for specified Bowtie2Index:\n    {bowtie2_index}\n")
            effective_ref_genome = ref_genome
            effective_bowtie2_index = self.__validate_index_prefix(
                bowtie2_index,
                "Bowtie2 index",
                [
                    (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
                    (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"),
                ],
            )
        else:
            if ref_genome is not None:
                self.set_ref_genome(ref_genome)
            effective_ref_genome = self.__get_ref_genome()
            effective_bowtie2_index = self.__get_bowtie2_index()

        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + "_" + effective_ref_genome
        if self.raw_data_type == "BAM":
            batch.append('### %s was skipped ###\n\n' % step_id)
            batch.append("# %s can only run with FASTQ files.\n#\t%s\n" % (step_id, "\n#\t".join(self.file_list)))
            return
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('mapper="%s"\n' % self.mapper)
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands (timing, index, SAM)
        batch.append('ref_genome="%s"\n' % effective_ref_genome)
        batch.append('Bowtie2Index="%s"\n' % effective_bowtie2_index)
        if SAM is None:
            batch.append('SAM="${SID}.${mapper}.${ref_genome}.sam"\n')
        else:
            batch.append('SAM="%s"\n' % SAM)
        batch.append("start=$(date +%s)\n")  # timer
        batch.append('bowtie2 \\\n')
        batch.append('    --rg-id "$SID" \\\n')
        batch.append('    --rg "SM:$SID" \\\n')
        batch.append('    --rg "LB:$SID" \\\n')
        batch.append('    --rg "PL:%s" \\\n' % (self.platform,))
        batch.append('    -x "$Bowtie2Index" \\\n')
        batch.append('    --threads %d \\\n' % self.threadN)
        if self.library_layout.lower() == "single":
            batch.append('    -U "$FASTQ" \\\n')
        elif self.library_layout.lower() == "paired":
            batch.append('    -1 "$FASTQ_R1" \\\n')
            batch.append('    -2 "$FASTQ_R2" \\\n')
        #batch.append('    -q \\\n')
        batch.append('    --phred33 \\\n')
        if end2end:
            batch.append('    --end-to-end  \\\n')
        batch.append('    -S "${SAM}"\n\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append("end=$(date +%s)\n")  # timer
        batch.append("runtime=$((end - start))\n")  # timer
        batch.append(
            f'echo -e "$SID\\t$(hostname)\\t$mapper\\t$ref_genome\\t{self.threadN}\\t$runtime\\tseconds\\t$(date +\'%Y-%m-%d %H:%M:%S\')" >> ../runtime.txt\n')  # timer
        self.RG_is_added = True
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        self.set_current_step_id(step_id) # mapping is a critical step, set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def minimap2(self, ref_genome=None, minimap2_index=None, SAM=None, force=False):
        batch = []
        self.mapper = "minimap2"
        if minimap2_index is not None:
            if ref_genome is None:
                sys.exit(f"Please set ref_genome for specified minimap2Index:\n    {minimap2_index}\n")
            effective_ref_genome = ref_genome
            effective_minimap2_index = self.__validate_existing_file(minimap2_index, "minimap2 index")
        else:
            if ref_genome is not None:
                self.set_ref_genome(ref_genome)
            effective_ref_genome = self.__get_ref_genome()
            effective_minimap2_index = self.__get_minimap2_index()

        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + "_" + effective_ref_genome
        if self.raw_data_type == "BAM":
            batch.append('### %s was skipped ###\n\n' % step_id)
            batch.append("# %s can only run with FASTQ files.\n#\t%s\n" % (step_id, "\n#\t".join(self.file_list)))
            return
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('mapper="%s"\n' % self.mapper)
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands (timing, index, SAM)
        batch.append('ref_genome="%s"\n' % effective_ref_genome)
        batch.append('minimap2Index="%s"\n' % effective_minimap2_index)
        if SAM is None:
            batch.append('SAM="${SID}.${mapper}.${ref_genome}.sam"\n')
        else:
            batch.append('SAM="%s"\n' % SAM)
        batch.append("start=$(date +%s)\n")  # timer
        batch.append('minimap2 \\\n')
        batch.append('    -a \\\n') # output in the SAM format
        batch.append('    -R "@RG\\tID:$SID\\tSM:$SID\\tLB:$SID\\tPL:%s" \\\n' % (self.platform,))
        batch.append('    -t %d \\\n' % self.threadN)
        batch.append('    -x sr \\\n') # short reads against a reference
        batch.append('    "${minimap2Index}" \\\n')
        if self.library_layout.lower() == "single":
            batch.append('    "$FASTQ" \\\n')
        elif self.library_layout.lower() == "paired":
            batch.append('    "$FASTQ_R1" \\\n')
            batch.append('    "$FASTQ_R2" \\\n')
        batch.append(' > "${SAM}"\n\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n' % step_id)
        batch.append(
            'echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append("end=$(date +%s)\n")  # timer
        batch.append("runtime=$((end - start))\n")  # timer
        batch.append(
            f'echo -e "$SID\\t$(hostname)\\t$mapper\\t$ref_genome\\t{self.threadN}\\t$runtime\\tseconds\\t$(date +\'%Y-%m-%d %H:%M:%S\')" >> ../runtime.txt\n')  # timer
        self.RG_is_added = True
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        self.set_current_step_id(step_id) # mapping is a critical step, set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def HISAT2(self, ref_genome=None, hisat2_index=None, SAM=None, force=False):
        batch = []
        self.mapper = "HISAT2"
        if hisat2_index is not None:
            if ref_genome is None:
                sys.exit(f"Please set ref_genome for specified HISAT2Index:\n    {hisat2_index}\n")
            effective_ref_genome = ref_genome
            effective_hisat2_index = self.__validate_index_prefix(
                hisat2_index,
                "HISAT2 index",
                [
                    tuple(f".{i}.ht2" for i in range(1, 9)),
                    tuple(f".{i}.ht2l" for i in range(1, 9)),
                ],
            )
        else:
            if ref_genome is not None:
                self.set_ref_genome(ref_genome)
            effective_ref_genome = self.__get_ref_genome()
            effective_hisat2_index = self.__get_hisat2_index()

        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + "_" + effective_ref_genome
        if self.raw_data_type == "BAM":
            batch.append('### %s was skipped ###\n\n' % step_id)
            batch.append(
                "# %s can only run with FASTQ files.\n#\t%s\n" % (step_id, "\n#\t".join(self.file_list)))
            return
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('mapper="%s"\n' % self.mapper)
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands (timing, index, SAM)
        batch.append('ref_genome="%s"\n' % effective_ref_genome)
        batch.append('HISAT2Index="%s"\n' % effective_hisat2_index)
        if SAM is None:
            batch.append('SAM="${SID}.${mapper}.${ref_genome}.sam"\n')
        else:
            batch.append('SAM="%s"\n' % SAM)
        batch.append("start=$(date +%s)\n")  # timer
        batch.append('hisat2 \\\n')
        batch.append('    --rg-id "$SID" --rg "SM:$SID" --rg "LB:$SID" --rg "PL:%s" \\\n' % (self.platform,))
        batch.append('    -p %d \\\n' % self.threadN)
        batch.append('    -x "$HISAT2Index" \\\n')
        if self.library_layout.lower() == "single":
            batch.append('    -U "$FASTQ" \\\n')
        elif self.library_layout.lower() == "paired":
            batch.append('    -1 "$FASTQ_R1" \\\n')
            batch.append('    -2 "$FASTQ_R2" \\\n')
        batch.append('    -S "${SAM}"\n\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n' % step_id)
        batch.append(
            'echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append("end=$(date +%s)\n")  # timer
        batch.append("runtime=$((end - start))\n")  # timer
        batch.append(
            f'echo -e "$SID\\t$(hostname)\\t$mapper\\t$ref_genome\\t{self.threadN}\\t$runtime\\tseconds\\t$(date +\'%Y-%m-%d %H:%M:%S\')" >> ../runtime.txt\n')  # timer
        self.RG_is_added = True
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        self.set_current_step_id(step_id) # mapping is a critical step, set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def STAR(self, star_index=None, ref_genome=None, gene_model=None, gtf = None, quantification=False, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        self.mapper = "STAR"
        if star_index is not None:
            if ref_genome is None:
                sys.exit("ERROR: ref_genome must be provided")
            effective_ref_genome = ref_genome
            effective_star_index = self.__validate_existing_directory(star_index, "STAR index")
        else:
            if ref_genome is not None:
                self.set_ref_genome(ref_genome)
            effective_ref_genome = self.__get_ref_genome()
            effective_star_index = self.__get_star_index()
        step_id = step_id + "_" + effective_ref_genome

        if gtf is None:
            #gtf = self.__get_gtf()
            effective_gtf = None
            gene_model = None # set gene_model to None as the default GTF will be used
        else:
            effective_gtf = self.__validate_existing_file(gtf, "GTF")
            if gene_model is None:
                gene_model = Path(effective_gtf).name.split(".")[0]
        if gene_model is not None:
            step_id = step_id + "_" + gene_model

        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('mapper="%s"\n' % self.mapper)
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('ref_genome="%s"\n' % effective_ref_genome)
        batch.append('STARIndex="%s"\n' % effective_star_index)
        if gene_model is not None:
            batch.append('gene_model="%s"\n' % gene_model)
        if effective_gtf:
            batch.append("GTF=%s\n" % effective_gtf)
        batch.append('rm -fr "${SID}.${mapper}.${ref_genome}%s.*"\n\n' % ("_{gene_model}" if gene_model else ""))
        # batch.append('FASTQ_unmapped_R1="${SID}.${mapper}.${ref_genome}.unmapped.R1.fastq.gz"\n')
        # batch.append('FASTQ_unmapped_R2="${SID}.${mapper}.${ref_genome}.unmapped.R2.fastq.gz"\n\n')
        batch.append('SAM="${SID}.${mapper}.${ref_genome}%s.sam"\n' % ("_${gene_model}" if gene_model else ""))
        batch.append('STAR_outFileNamePrefix="${SAM%.sam}"\n')
        batch.append("start=$(date +%s)\n")  # timer
        batch.append('STAR \\\n')
        batch.append('    --runMode alignReads \\\n')
        if quantification:
            batch.append('    --quantMode TranscriptomeSAM GeneCounts \\\n')
        batch.append('    --twopassMode Basic \\\n')
        batch.append('    --runThreadN %d \\\n' % (self.threadN,))
        batch.append('    --genomeDir "$STARIndex" \\\n')
        if self.library_layout.lower() == "single":
            batch.append('    --readFilesIn "$FASTQ" \\\n')
        elif self.library_layout.lower() == "paired":
            batch.append('    --readFilesIn "$FASTQ_R1" "$FASTQ_R2" \\\n')
        batch.append('    --readFilesCommand zcat \\\n')
        if effective_gtf:
            batch.append('    --sjdbGTFfile "$GTF" \\\n')
        batch.append('    --outFileNamePrefix "${STAR_outFileNamePrefix}." \\\n')
        batch.append('    --outSAMattrRGline "ID:$SID\\tSM:$SID\\tLB:$SID\\tPL:%s" \\\n' % (self.platform,))
        batch.append('    --outSAMtype SAM \\\n')
        batch.append('    --outSAMunmapped Within KeepPairs \\\n')
        batch.append('    --outReadsUnmapped None\n\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n' % step_id)
        batch.append('mv "${STAR_outFileNamePrefix}.Aligned.out.sam" "${SAM}"\n')
        batch.append(
            'echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append("end=$(date +%s)\n")  # timer
        batch.append("runtime=$((end - start))\n")  # timer
        batch.append(
            f'echo -e "$SID\\t$(hostname)\\t$mapper\\t$ref_genome\\t{self.threadN}\\t$runtime\\tseconds\\t$(date +\'%Y-%m-%d %H:%M:%S\')" >> ../runtime.txt\n')  # timer
        self.RG_is_added = True
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        self.set_current_step_id(step_id) # mapping is a critical step, set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def taxomomic_classifier(self, contig_info, taxonomy_map=None, virus_group_by=None,force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('python /account001/bgong/workspace/AdVentSeq/taxonomic_classifier.py \\\n')
        batch.append('    --bam "${BAM}" \\\n')
        batch.append('    --combined_genome %s \\\n' % self.ref_genome)
        if taxonomy_map:
            batch.append('    --taxonomy_map "%s" \\\n' % taxonomy_map)
        if virus_group_by:
            batch.append('    --virus_group_by "%s" \\\n' % virus_group_by)
        batch.append('    --contig_info "%s"\n' % contig_info)
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def remove_read_pairs_mapped_to_host(self, use_samtools=True, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        # sort SAM file by query name before this step
        #self.sort_SAM(force=force, by_qname=True)
        # continue this step
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        if use_samtools:
            batch.append('samtools view -b -f 12 -F 2304 -@ %d "${BAM}" \\\n' % self.threadN)
            # -f 12: Ensures that at least one read in the pair is unmapped (READ1_UNMAPPED or READ2_UNMAPPED)
            # -F 256: Excludes secondary alignments
            # -F 2048: Excludes supplementary alignments
            # -F 2048 + 256 = 2304, so we use -F 2304 to exclude both secondary and supplementary alignments
            batch.append(' | samtools sort -n -@ %d \\\n' % self.threadN)
            batch.append(' | samtools fastq -n -@ %d \\\n' % self.threadN)
            batch.append('    -1 "${BAM%.*.*}.unmapped2host.R1.fastq.gz" \\\n')
            batch.append('    -2 "${BAM%.*.*}.unmapped2host.R2.fastq.gz" \\\n')
            batch.append('    -0 "${BAM%.*.*}.unmapped2host.single.fastq.gz" \\\n')
            batch.append('    -s "${BAM%.*.*}.unmapped2host.orphan.fastq.gz"\n')
            """
            batch.append('samtools view -b -f 12 -F 256 -@ %d \\\n' % self.threadN)
            # -f 12: Ensures that at least one read in the pair is unmapped (READ1_UNMAPPED or READ2_UNMAPPED)
            # -F 256: Excludes secondary alignments
            batch.append('    "${BAM}" \\\n')
            batch.append('    > "${BAM%.*}.unmapped2host.bam"\n')
            batch.append('if [ $? -ne 0 ]; then echo "Error: %s | samtools remove host reads failed."; exit 1; fi\n\n' % step_id)

            batch.append('samtools sort -n -@ %d \\\n' % self.threadN)
            batch.append('    -o "${BAM%.*}.unmapped2host.sorted_by_qname.bam" \\\n')
            batch.append('    "${BAM%.*}.unmapped2host.bam"\n')
            batch.append('if [ $? -ne 0 ]; then echo "Error: %s | samtools sort by qname failed."; exit 1; fi\n\n' % step_id)

            batch.append('samtools fastq -n -@ %d \\\n' % self.threadN)
            batch.append('    -1 "${BAM%.*}.unmapped2host.R1.fastq.gz" \\\n')
            batch.append('    -2 "${BAM%.*}.unmapped2host.R2.fastq.gz" \\\n')
            batch.append('    -0 "${BAM%.*}.unmapped2host.single.fastq.gz" \\\n')
            batch.append('    -s "${BAM%.*}.unmapped2host.orphan.fastq.gz" \\\n')
            batch.append('    "${BAM%.*}.unmapped2host.sorted_by_qname.bam"\n')
            batch.append('if [ $? -ne 0 ]; then echo "Error: %s | samtools convert to fastq failed."; exit 1; fi\n\n' % step_id)
            """
        else:
            batch.append('python /account/bgong/workspace/Project_DeLab/remove_host_read_pairs.py \\\n')
            batch.append('    --bam "${BAM}" \\\n')
            batch.append('    --prefix "${BAM%.*.*}.unmapped2host"\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n')
        self.host_read_pairs_removed = True
        batch.append('SID="${BAM%.*.*}.unmapped2host"\n')
        batch.append('FASTQ_R1="${SID}.R1.fastq.gz"\n')
        batch.append('FASTQ_R2="${SID}.R2.fastq.gz"\n')
        batch.append('\n')
        # test if this step has already been completed
        self.set_current_step_id(step_id) # this step set SID and FASTQ_R1/2, must set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def Kraken2(self, db_name="standard", kraken2_db_root="/galaxy001/Resources/Kraken2DB", force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + "_" + db_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append(f'Kraken2DB_ROOT="{kraken2_db_root}"\n')
        batch.append(f'kraken2db="{db_name}"\n')
        batch.append('Kraken2Report="${SID}.Kraken2_${kraken2db}.report"\n')
        batch.append('if [ ! -d "${Kraken2DB_ROOT}/${kraken2db}" ]; then echo "Error: kraken2db (%s) not exists."; exit 1; fi\n\n' % db_name)
        batch.append('kraken2 \\\n')
        batch.append('    --db "${Kraken2DB_ROOT}/${kraken2db}" \\\n')
        batch.append('    --threads %d \\\n' % self.threadN)
        batch.append('    --unclassified-out "${SID}.Kraken2_${kraken2db}.unclassified#.fastq" \\\n')
        batch.append('    --classified-out "${SID}.Kraken2_${kraken2db}.classified#.fastq" \\\n')
        batch.append('    --output "${SID}.Kraken2_${kraken2db}.out" \\\n')
        batch.append('    --report "${Kraken2Report}" \\\n')
        batch.append('    --paired \\\n')
        batch.append('    --use-names \\\n')
        batch.append('    --gzip-compressed \\\n')
        batch.append('    "${FASTQ_R1}" "${FASTQ_R2}"\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('find . -name "*.fastq" -exec gzip {} \\;\n')
        batch.append('gzip "${SID}.Kraken2_${kraken2db}.out"\n')
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def featureCounts(self, gene_model=None, feature="exon", countReadPairs = True, saf = None, gtf = None, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        if countReadPairs:
            step_id += "_count_pairs"
        else:
            step_id += "_count_reads"

        annotation_format = "GTF"
        if saf and gtf:
            sys.exit("Please specify either GTF or SAF for featureCounts, but not both.\n")
        if saf:
            gene_model = gene_model or Path(saf).name.split(".")[0]
            annotation_format = "SAF"
        elif gtf:
            gene_model = gene_model or Path(gtf).name.split(".")[0]
        else:
            gtf = self.__get_gtf()
            gene_model = gene_model or Path(gtf).name.split(".")[0]

        step_id = f"{step_id}|{annotation_format}_{gene_model}_{feature}"
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # set annotation
        if annotation_format == "SAF":
            batch.append('GENE_MODEL="%s"\n' % saf)
        else:
            batch.append('GENE_MODEL="%s"\n' % gtf)
        # commands
        batch.append('rm -f *.tmp\n\n')

        batch.append('CountRes="${SID}.${mapper}.${ref_genome}.%s.%s.%s.tsv"'
                                  '\n' % (fn_name,
                                          "_".join([gene_model, feature]),
                                          "count_pairs" if countReadPairs else "count_reads"))
        batch.append('featureCounts \\\n')
        batch.append('    -F %s \\\n' % annotation_format)
        if feature in ("exon","transcript","CDS"):
            batch.append('    -t %s \\\n' % feature)
        else:
            sys.exit("Wrong feature: %s\n" % feature)
        self.__update_env_variable("feature", feature)
        # batch.append('    -g gene_id \\\n') # gene_id is the default meta-feature
        if self.library_layout.lower() == "paired":
            batch.append('    -p \\\n') # libraries are assumed to contain paired-end reads
            if countReadPairs:
                batch.append('    --countReadPairs \\\n')  #  fragments (or templates) will be counted
        batch.append('    --primary \\\n')  # Count primary alignments only.
        batch.append('    -B \\\n') #  Only count read pairs that have both ends aligned
        batch.append('    -C \\\n')  #  Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome  but on different strands
        if self.sorted_by_qname:
            batch.append('    --donotsort \\\n') #  Do not sort reads in BAM/SAM input
        batch.append('    -T %d \\\n' % self.threadN)
        batch.append('    -a "${GENE_MODEL}" \\\n')
        batch.append('    -o "${CountRes}" \\\n')
        batch.append('    "${sorted_BAM}"\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('gzip -f "${CountRes}"\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - gzip failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def depth_by_pos(self, min_MAPQ=1, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('samtools depth \\\n')
        batch.append('    -q %d  \\\n' % min_MAPQ)
        batch.append('    "${sorted_BAM}" \\\n')
        batch.append(' | gzip > "${BAM}_depth.txt.gz"\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def coverage_from_BigWig(self, bed_file, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('python /account/bgong/workspace/Project_DeLab/bg2table.py \\\n')
        batch.append('    --cov "${SID}.${mapper}.${ref_genome}.coverageBed.gz" \\\n')
        batch.append('    --out_prefix "${SID}.${mapper}.${ref_genome}.boundary200"\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def BAM2BigWig(self, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name]) # add env name to required env list
        # commands
        batch.append('assigned_reads=$(awk \'NR==2 {print $2}\''
                                  ' "${CountRes}.summary")\n')
        batch.append('scale_factor=$(awk "BEGIN { printf \\\"%.10f\\\", 1000000 / ${assigned_reads} }")\n')
        batch.append('echo "scale_factor=$scale_factor"\n')
        batch.append('bamCoverage \\\n')
        batch.append('   --bam "${sorted_BAM}" \\\n')
        batch.append('   --outFileName "${SID}.${mapper}.${ref_genome}.bw" \\\n')
        batch.append('   --outFileFormat bigwig \\\n')
        batch.append('   --scaleFactor ${scale_factor} \\\n')
        batch.append('   --binSize 1 \\\n')
        batch.append('   --smoothLength 1 \\\n')
        batch.append('   --numberOfProcessors %d \\\n' % self.threadN)
        batch.append('   --normalizeUsing None \\\n')
        batch.append('   --minMappingQuality 1 \\\n')
        batch.append('   --skipNonCoveredRegions\n')

        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def unmapped_to_fastq(self):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name]) # add env name to required env list
        # commands
        batch.append('FASTQ_unmapped_R1="${SID}.${mapper}.${ref_genome}_unmapped.R1.fastq.gz"\n')
        batch.append('FASTQ_unmapped_R2="${SID}.${mapper}.${ref_genome}_unmapped.R2.fastq.gz"\n')
        batch.append('FASTQ_unmapped_S="${SID}.${mapper}.${ref_genome}_unmapped.singleton.fastq.gz"\n')
        batch.append('samtools fastq \\\n')
        batch.append('    -s "$FASTQ_unmapped_S" \\\n')
        batch.append('    -1 "$FASTQ_unmapped_R1" \\\n')
        batch.append('    -2 "$FASTQ_unmapped_R2" \\\n')
        batch.append('    -f 0x4 \\\n')  # 0x4	UNMAP	segment unmapped
        batch.append('    -@ %d \\\n' % self.threadN)
        batch.append('    "$BAM"\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - samtools fastq failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # to be decided whether to set current step id
        if step_id in self.progress:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def __BAM_to_proper_paired_FASTQ_of_a_contig(self, contig, contig_name, force=False):
        batch = []
        if contig_name in self.FASTQ_ready_for_contigs:
            return
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + f"[{contig_name}]"
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append(f'FASTQ_CONTIG_R1="${{sorted_BAM%.*.*}}.{contig_name}.R1.fastq.gz"\n')
        batch.append(f'FASTQ_CONTIG_R2="${{sorted_BAM%.*.*}}.{contig_name}.R2.fastq.gz"\n')
        #batch.append(f'BAM_CONTIG="${{sorted_BAM%.*.*}}.{contig_name}.bam"\n')
        batch.append('samtools view \\\n')
        batch.append('    -b \\\n')
        batch.append('    -f 3 \\\n') # include: paired (1) + porper pair (2)
        batch.append('    -F 3852 \\\n') # remove: read unmapped (4) + mate unmapped (8) + secondary (256) + QC failed (512) + PCR duplicates (1024) +  supplementary (2048)
        batch.append('    ${sorted_BAM} \\\n')
        batch.append(f'    {contig} | \\\n')
        batch.append('samtools collate -O -u - | \\\n') # -O Output to stdout, -u uncompressed
        batch.append('samtools fastq \\\n')
        batch.append('    -1 ${FASTQ_CONTIG_R1} \\\n')
        batch.append('    -2 ${FASTQ_CONTIG_R2} \\\n')
        batch.append('    -0 /dev/null \\\n')
        batch.append('    -s /dev/null \\\n')
        batch.append('    -n -\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # do not set step_id, this step is an internal function and is called by SPAdes and others.
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        # a little different from other functions. this internal function may be called multiple times.
        # add contig_name to set for testing to prevent duplicate script lines
        self.batch[step_id].extend(batch)
        self.FASTQ_ready_for_contigs.add(contig_name)


    def SPAdes(self, contig, mode="rnaviral", contig_name=None, force=False):
        batch = []
        if mode not in ("rna","rnaviral"):
            mode = "general_careful"
        if not contig_name:
            contig_name = str(contig)
        # call function to prepare FASTQ files for the contig    
        self.__BAM_to_proper_paired_FASTQ_of_a_contig(contig=contig, contig_name=contig_name, force=force)
        # SPAdes
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name + f"[{contig_name}]_{mode}"
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name]) # add env name to required env list
        # PE
        batch.append(f'SPAdes_Output="SPAdes/{contig_name}/{mode}"\n')
        batch.append('[ -d "${SPAdes_Output}" ] && rm -fr "${SPAdes_Output}"\n')
        batch.append('mkdir -p "${SPAdes_Output}"\n')
        batch.append('spades.py \\\n')
        batch.append('    -o "${SPAdes_Output}" \\\n')
        if mode == "rnaviral":
            batch.append('    --rnaviral \\\n')
        elif mode == "rna":
            batch.append('    --rna \\\n')
        else:
            batch.append('    --careful \\\n')
        batch.append('    -1 "$FASTQ_CONTIG_R1" \\\n')
        batch.append('    -2 "$FASTQ_CONTIG_R2" \\\n')
        batch.append('    --threads %d \\\n' % self.threadN)
        batch.append('    --memory 64\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s - PE failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def GATK4_Mutect2(self, ref_genome=None, genome_fasta=None, mode="DNA", no_filter=False, split_multi_allelic=False,force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        self.caller = "GATK4_Mutect2"
        if genome_fasta is not None:
            if ref_genome is None:
                sys.exit(f"Please set ref_genome for specified genome_fasta:\n    {genome_fasta}\n")
            effective_ref_genome = ref_genome
            effective_genome_fasta = self.__validate_existing_file(genome_fasta, "genome FASTA")
        else:
            if ref_genome is not None:
                self.set_ref_genome(ref_genome)
            effective_ref_genome = self.__get_ref_genome()
            effective_genome_fasta = self.__get_genome_fasta()

        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('caller="%s"\n' % self.caller)
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        ### commands
        batch.append('ref_genome="%s"\n' % effective_ref_genome)
        batch.append('genome_fasta="%s"\n' % effective_genome_fasta)
        batch.append('VCF="${sorted_BAM%.*}.vcf"\n')
        # call variants
        batch.append('gatk Mutect2 \\\n')
        batch.append('    -R "${genome_fasta}" \\\n')
        batch.append('    -I "${sorted_BAM}" \\\n')
        if mode == "RNA":
            batch.append('    --max-mnp-distance 0 \\\n') # Prevents MNPs from being grouped, ensential for RNA-seq
            batch.append('    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\\n') # Prevents loss of useful RNA-seq reads that are properly mapped but don’t conform to DNA-seq pairing expectations
        batch.append('    --output "${VCF}"\n')
        batch.append(' > "${SAM}"\n\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s|caller failed."; exit 1; fi\n' % step_id)
        # filter
        batch.append('gatk FilterMutectCalls \\\n')
        batch.append('    -R "${genome_fasta}" \\\n')
        batch.append('    -V "${VCF}" \\\n')
        batch.append('    --output "${VCF%.*}.filtered.vcf"\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s|filter failed."; exit 1; fi\n' % step_id)
        # VCF to table
        batch.append('gatk VariantsToTable \\\n')
        batch.append('    -V "${VCF%.*}.filtered.vcf" \\\n')
        batch.append('    -F CHROM -F POS -F REF -F ALT -F FILTER -F DP -GF DP -GF AF \\\n')
        if split_multi_allelic:
            batch.append('    --split-multi-allelic \\\n')
        if no_filter:
            batch.append('    --show-filtered \\\n')
            batch.append('    -O "${VCF%.*}.nofilter.tsv"\n')
        else:
            batch.append('    -O "${VCF%.*}.filtered.tsv"\n')
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s|to table failed."; exit 1; fi\n' % step_id)

        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def ViraQuant(self, virus_list=None, top_n=10, scan_by=None, force=False):
        batch = []
        surfix = ""
        if scan_by:
            surfix = "_scan"
        elif virus_list:
            surfix = ""
        else:
            surfix = f"_top{top_n}"
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + f"{fn_name}{surfix}"
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append('python /account001/bgong/workspace/AdVentSeq/ViraQuant.py \\\n')
        batch.append('    --bam "${sorted_BAM}" \\\n')
        if scan_by:
            batch.append(f'    --scan-by "{scan_by}" \\\n')
        elif virus_list:
            batch.append(f'    --viruses "{virus_list}" \\\n')
        else:
            batch.append(f'    --top-n {top_n} \\\n')
        batch.append(f'    --out "${{sorted_BAM%.*.*}}.ViraQuant{surfix}.tsv"\n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        # this is an endpoint, no need to set current step id
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)


    def __example(self, force=False):
        batch = []
        fn_name = inspect.currentframe().f_code.co_name  # get the name of the function
        step_id = self.__current_step_id + "|" + fn_name
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('### %s ###\n' % step_id)
        batch.append(f'{"#" * (len(step_id) + 8)}\n')
        batch.append('conda activate %s\n\n' % self.envs4steps[fn_name])
        self.required_conda_envs.add(self.envs4steps[fn_name])  # add env name to required env list
        # commands
        batch.append(' \\\n')
        batch.append('    \\\n')
        batch.append('    \\\n')
        batch.append('    \\\n')
        batch.append('    \\\n')
        batch.append('    \\\n')

        batch.append('    \n')
        # check if commands were completed successfully
        batch.append('if [ $? -ne 0 ]; then echo "Error: %s failed."; exit 1; fi\n\n' % step_id)
        batch.append('echo -e "%s\\t$(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')" >> "$my_progress"\n\n' % step_id)
        batch.append('conda deactivate\n\n\n')
        # test if this step has already been completed
        self.set_current_step_id(step_id)
        if step_id in self.progress and not force:
            batch = self.__comment_lines(batch)
        else:
            self.is_completed = False
        self.batch[step_id].extend(batch)
