# -*- coding: utf-8 -*-
# @Time    : 2021/8/26 11:35
# @Author  : Zhongyi Hua
# @File    : utilities.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
import sys
from typing import Dict

from .global_params import SUB_STEP_GATK, SUB_STEP_SAMTOOLS
from .utilities import split_interval


def pp_out_dir(_output_dir: Path, _sub_dir=''):
    if not _output_dir.exists():
        _output_dir.mkdir()
    else:
        if _output_dir.is_file():
            print("The output directory seems a file?")
            sys.exit(1)
    # make some output directory
    try:
        (_output_dir/ _sub_dir).mkdir()
    except FileExistsError:
        pass


@dataclass
class SampleIO:
    """
    Attributes and Value: file_type and Path
    
    Possible file_type:
    fq1,fq2,cleanfq1,cleanfq2,bam,gvcf,vcf,
    
    **IMPORTANCE:** All SampleIO instance should be created in `pp_meta*` func
    """
    name: str = None


@dataclass
class StepIO:
    """
    Although some attributes seems duplicated with the ConfigPipe, I believe it is necessary

    outsh: Path to write the script
    in_lst: [{'{sample}': sample1, '{file_type}': sample1_file_path1, 'file_type2': ..., <Maybe other keys>},
             {'sample': sample2, '_name: sample2_file_path1},
             ...]
    out_lst: [{}]
    """
    step_name: str # Required
    software: str = None
    ref: str = None
    outsh: Path = None
    outdir: Path = None
    logdir: Path = None
    outdct: Dict = field(default_factory=dict) # For most time, the key should be the sample name


class ConfigPipe:
    def __init__(self, config:dict , outdir:Path , logdir:Path, project: str, units: list, task:str):
        """

        Args:
            config: yaml.load results
            units[list]: e.g. ['all'] 

        Attrs:
            ConfigPipe.project = 'project name'
            ConfigPipe.task = 'RNA/WGS'
            ConfigPipe.units = ['QC', 'align', 'call']
            ConfigPipe.otherInfo = {step1: {Info1: ..., Info2:...}, ...} # Key-values other than 'software', 'bin', and 'params' in config file
            ConfigPipe.software = {step1: software1, step2: software2}
            ConfigPipe.softwareBin = {software1: bin1, software2: bin2}
            ConfigPipe.softwarePara = {software1: {para1: ..., para2: ...},
                                       software2: {para1: ..., para2: ...}}
            ConfigPipe.outdir = Path(The root outdir directory)
            ConfigPipe.logdir = Path(The root logdir directory)
            ConfigPipe.tmp = Path(The root tmp directory)
            ConfigPipe.ref =  {software1: Ref1, software2: Ref2, ...}
        """
        self.project: str = project
        self.task: str = task
        self.units: deque = self.sort_unit(task, units)
        self.software: dict = self._parse_software(config, self.units)
        self.softwareBin: dict = self._parse_bin(config)
        self.softwarePara: dict = self._parse_para(config, self.units)
        self.otherInfo: dict = self._parse_info(config, self.units)
        self.outdir: Path = outdir
        self.logdir: Path = logdir
        self.tmp: Path = self._pp_tmp(config, outdir)
        self.ref: Path = self._pp_ref(config, project, outdir)
    
    @staticmethod
    def sort_unit(task, units):
        """_summary_

        Args:
            task (_type_): _description_
            units (_type_): _description_
        Return:
            Ordered pipeline unit
        """
        if task == 'WGS':
            start_dict = {'QC': 0, 'align': 1, 'call': 2}
            if units == ['all']:
                units = ['QC', 'align', 'call']
        else: # RNA
            start_dict = {'QC': 0, 'align': 1, 'quantify': 2}
            if units == ['all']:
                units = ['QC', 'align', 'quantify']
        units = sorted(units, key=lambda x: start_dict[x])
        return units

    @staticmethod
    def _parse_software(_config, _units):
        units_dict = defaultdict()
        for _step in _units:
            units_dict[_step] = _config[_step]['software']
        return units_dict
    
    @staticmethod
    def _parse_para(_config, _units):
        """Returns:
        Legal key for now: threads, HaplotypeCaller, CombineGVCFs, GenotypeGVCFs, OtherParams
        {
        software1: {threads: 2, otherPara: []} 
        GATK:      {threads: 2, HaplotypeCaller: [], ...}
        }
        """
        para_dict = defaultdict()
        for _step in _units:
            soft_name = _config[_step]['software']
            soft_dct = defaultdict()
            if _config[_step].get('params'):
                if otherparams := _config[_step]['params'].get(soft_name):
                    soft_dct['OtherParams'] = [str(_) for _ in otherparams] # Change number to str
                if threads := _config[_step]['params'].get("threads", '1'):
                    soft_dct['threads'] = threads
            para_dict[soft_name] = soft_dct
        return para_dict

    @staticmethod
    def _parse_bin(_config):
        soft_dict = defaultdict()
        for _para_dct in _config.values():
            if bin_dct := _para_dct.get('bin'):
                bin_dct = {_software: Path(_binpath) for _software, _binpath in bin_dct.items()}
                soft_dict.update(bin_dct)
        return soft_dict
    
    @staticmethod
    def _parse_info(_config, units):
        """Parse key-values other than 'software', 'bin', and 'params'.
        e.g. intervals for 'call' step.

        Args:
            _config (_type_): _description_
            units (_type_): _description_

        Returns:
            _type_: _description_
        """
        info_dict = defaultdict()
        for _step in units:
            info_dict[_step] = {_key: _value
                                for _key, _value in _config[_step].items()
                                if _key not in {'software', 'bin', 'params'}}
        return info_dict

    @staticmethod
    def _pp_tmp(_config, _outdir):
        # set temporary directory
        if _config.get('resources').get('temporary') is None:
            _tmp_path = (_outdir / 'tmp')
        else:
            _tmp_path = Path(_config['resources']['temporary'])
        pp_out_dir(_tmp_path)
        return _tmp_path

    @staticmethod
    def _pp_ref(_config, _project, _outdir: Path):
        """
        Only ref file path for Index step is generated in ConfigPipe. Because no matter
        which steps you choose, it is necessary to index. Therefore, Index always step0. 
        """
        ref_dct = defaultdict()

        # cmds path
        ref_dct['outsh'] = _outdir / ( _project + '.step0.index.sh')
        # Raw Input
        ref_dct['raw_ref'] = _config['resources'].get('reference')

        # Output
        indexdir = _outdir / '00_index'
        ref_index_path = indexdir / Path(ref_dct['raw_ref']).name

        if _tmp_path := _config['resources'].get('gff'):
            ref_dct['raw_gff'] = _tmp_path
            ref_dct['gff'] = indexdir / Path(ref_dct['raw_gff']).name
        if _tmp_path := _config['resources'].get('gtf'):
            ref_dct['raw_gtf'] = _tmp_path
            ref_dct['gtf'] = indexdir / Path(ref_dct['raw_gtf']).name
        else:
            # Just set an empty value to avoid raising error in index_dct
            ref_dct['gtf'] = ''
        # Output
        pp_out_dir(indexdir)
        prefix = ref_index_path.stem
        # All possible index output path
        index_dct = {'hisat2': indexdir / prefix,
                     'STAR': indexdir, # STAR use genome dir directly
                     'RSEM': ref_index_path,
                     'htseq-count': indexdir / Path(ref_dct['gtf']).name,
                     'bwa': ref_index_path,
                     'bowtie': indexdir / prefix,
                     'bowtie2': indexdir / prefix,
                     'GATK': ref_index_path,
                     'bcftools': ref_index_path,
                     'samtools': ref_index_path
                     }
        #for _step in units:
        #    _software = _config[_step]['software']
        #    ref_dct[_software] = index_dct.get(_software)
        ref_dct['ref'] = ref_index_path
        ref_dct['refindex'] = indexdir / (Path(ref_dct['raw_ref']).name + '.fai') # .fai file for fasta file. It will be created by IndexPipe
        ref_dct.update(index_dct)
        return ref_dct


class MetaPipe:
    """
    All methods under MetaPipe (except for internal methods) should
    return a StepIO instance. 
    
    Notice: 
    1. All samples' outputs returned here refer to the str used in cmds,
     not the file name (e.g., STAR and RSEM accepeted a prefix to generate
      multiple output files). Therefore, some software need a convertion function.

    2. This class only return the final output file path. The interminate file will
    be handled in each PP_* class.
    """

    def __init__(self, meta: Path, pipe_config: ConfigPipe):
        self.meta = meta
        self.pipe_config = pipe_config
        pp_out_dir(pipe_config.outdir)
        # Notice! pp_meta always runs before all other
        # functions to generate samplelist
        # The only possible value for sample_dct is the RG
        self.sample_dct: dict = {}


    def _pp_meta1(self, step='QC'):
        """For type 1 meta file. Only for `QC` or `align`
        :return ["{sample1}": SampleIO, "{sample2}": SampleIO, ...]
        """
        sample_dct = defaultdict()
        with open(self.meta, 'r') as f_in:
            for _line in f_in.read().splitlines():
                if not _line.startswith('#'):
                    _sample_ins = SampleIO()
                    _sample_ins.name, _sample_ins.RG, _fq1, _fq2 = _line.split()
                    setattr(_sample_ins, 'fq1' if step == 'QC' else 'cleanfq1', Path(_fq1))
                    setattr(_sample_ins, 'fq2' if step == 'QC' else 'cleanfq2', Path(_fq2))
                    sample_dct[_sample_ins.name] = _sample_ins
                    # Notice! pp_meta always runs before all other functions to generate sample dict! 
                    self.sample_dct[_sample_ins.name] = _sample_ins.RG
        return sample_dct


    def _pp_meta2(self, file_type='bam'):
        """For type2 meta file

        Args:
            _file (Path): sample meta info path
            _name (str, optional): _description_. Defaults to 'bam'/'vcf'/'gvcf'/....

        Returns:
            [{'sample': sample1, _name: sample1_file_path1},
             {'sample': sample2, _name: sample2_file_path1},
             ...]
        """
        sample_dct = defaultdict()
        with open(self.meta, 'r') as f_in:
            for _ in f_in.read().splitlines():
                _sample_ins = SampleIO()
                _sample_ins.name, _file_path =  _.split('\t')
                setattr(_sample_ins, file_type, Path(_file_path))
                sample_dct[_sample_ins.name] = _sample_ins
                # Notice! pp_meta always runs before all other functions to generate sample list! 
                # For type2 meta, the RG is null.
                self.sample_dct[_sample_ins.name] = ''
        return sample_dct
    

    def pp_qc(self, order=1):
        step_io = StepIO('QC')
        step_io.outsh = self.pipe_config.outdir / (self.pipe_config.project + f'.step{order}.QC.sh')
        step_io.outdir = self.pipe_config.outdir / f'{order:02d}_cleandata'
        step_io.logdir = self.pipe_config.logdir / f'{order:02d}_QC'
        pp_out_dir(step_io.outdir)
        pp_out_dir(step_io.logdir)
        for _sample_name in self.sample_dct:
            _sample_ins = SampleIO(_sample_name)
            _sample_ins.RG = self.sample_dct.get(_sample_name)
            _sample_ins.cleanfq1 = step_io.outdir / f"{_sample_ins.name}_1_clean.fq.gz"
            _sample_ins.cleanfq2 = step_io.outdir / f"{_sample_ins.name}_2_clean.fq.gz"
            step_io.outdct[_sample_name] = _sample_ins
        return [step_io]

    def pp_align(self, order=1):
        step_io = StepIO('align')
        # Output
        step_io.outsh = self.pipe_config.outdir / (self.pipe_config.project + f'.step{order}.align.sh')
        step_io.outdir = self.pipe_config.outdir / f'{order:02d}_mapping'
        step_io.logdir = self.pipe_config.logdir / f'{order:02d}_mapping'
        pp_out_dir(step_io.outdir)
        pp_out_dir(step_io.logdir)
        # Generate different output path for different softwares
        match self.pipe_config.software['align']:
            # From now, RNA-seq
            case 'STAR':
                _suffix = ''
                # The above _suffix is used for alignment output. However, the actual output file in the BAM format
                # ends with 'toTranscriptome.out.bam'. Therefore, it has to set the software attribute, so that the
                # PPquan class could handle the change in suffix when the STAR aligner was used.
                step_io.software = 'STAR'
            case 'hisat2':
                _suffix = '.bam'
            # From now, WGS
            case 'bowtie' | 'bowtie2' | 'bwa':
                _suffix = '.sorted.bam'
        for _sample in self.sample_dct:
            _sample_ins = SampleIO(_sample)
            _sample_ins.bam = step_io.outdir / (_sample_ins.name + _suffix)
            step_io.outdct[_sample] = _sample_ins
        return [step_io]

    def pp_quan(self, order=1):
        step_io = StepIO('quantify')
        step_io2 = StepIO('MergeQuantify') # For merging results from samples
        step_io.outsh = self.pipe_config.outdir / (self.pipe_config.project + f'.step{order}.quan.sh')
        step_io.outdir = self.pipe_config.outdir / f'{order:02d}_expression'
        step_io.logdir = self.pipe_config.logdir / f'{order:02d}_expression'
        pp_out_dir(step_io.logdir)
        pp_out_dir(step_io.outdir)
        match self.pipe_config.software['quantify']:
            case 'RSEM':
                _suffix = ''
                _suffix_merge = 'sh'
                step_io.software = step_io2.software = 'RSEM'
            case 'htseq-count':
                _suffix = '.tsv'
                _suffix_merge = 'py'
                step_io.software = step_io2.software = 'htseq-count'
        # step1. For each sample
        for _sample in self.sample_dct:
            _sample_ins = SampleIO(_sample)
            _sample_ins.tsv = step_io.outdir / (_sample_ins.name + _suffix)
            step_io.outdct[_sample] = _sample_ins
        # step2. Merge samples
        step_io2.outsh = self.pipe_config.outdir / (self.pipe_config.project + f'.step{order+1}.merge.{_suffix_merge}')
        step_io2.outdir = self.pipe_config.outdir
        step_io2.logdir = self.pipe_config.logdir
        step_io2.outdct['Alltsv'] = self.pipe_config.outdir / 'count.tsv'
        return [step_io, step_io2]

    def pp_call(self, order=1):
        step_io_lst = []
        match self.pipe_config.software['call']:
            case 'GATK':
                for sub_step in SUB_STEP_GATK:
                    step_io = StepIO(sub_step)
                    # For output directory
                    if sub_step in ['gvcf', 'genotype']:
                        step_io.outdir = self.pipe_config.outdir / f"{order:02d}_{sub_step}"
                        step_io.logdir = self.pipe_config.logdir / f"{order:02d}_{sub_step}"
                    if sub_step.startswith('merge'):
                        step_io.outdir = self.pipe_config.outdir
                        step_io.logdir = self.pipe_config.logdir
                    pp_out_dir(step_io.outdir)
                    pp_out_dir(step_io.logdir)
                    # For generate script
                    step_io.outsh = self.pipe_config.outdir / (self.pipe_config.project + f'.step{order}.{sub_step}.sh')
                    if sub_step == 'gvcf':
                        for _sample_name in self.sample_dct:
                            _sample_ins = SampleIO(_sample_name)
                            _sample_ins.gvcf = step_io.outdir / f"{_sample_name}.g.vcf.gz"
                            step_io.outdct[_sample_name] = _sample_ins
                    if sub_step == 'merge_sample':
                        step_io.outdct['gvcfDB'] = step_io.outdir / f"{order:02d}_genomeDB.g.vcf.gz"
                    if sub_step == 'genotype':
                        intervals = self.pipe_config.otherInfo['call'].get('intervals')
                        if intervals is None:
                            with open(self.pipe_config.ref['refindex'], 'r', encoding='utf-8') as f_in:
                                intervals = [_.split()[0] for _ in f_in.read().splitlines()]
                        for _interval in intervals:
                            _sample_ins = SampleIO(_interval)
                            _sample_ins.vcf = step_io.outdir / f"{_interval}.vcf.gz"
                            step_io.outdct[_interval] = _sample_ins
                    if sub_step == 'merge_seq':
                        step_io.outdct['vcf'] = step_io.outdir / f"{order:02d}_result.vcf.gz"
                    step_io_lst.append(step_io)
                    order += 1
            case 'bcftools':
                for sub_step in SUB_STEP_SAMTOOLS:
                    step_io = StepIO(sub_step)
                    step_io.outsh = self.pipe_config.outdir / (self.pipe_config.project + f'.step{order}.{sub_step}.sh')
                    if sub_step == 'genotype':
                        step_io.outdir = self.pipe_config.outdir / f"{order:02d}_{sub_step}"
                        step_io.logdir = self.pipe_config.logdir / f"{order:02d}_{sub_step}"
                        intervals = self.pipe_config.otherInfo['call'].get('intervals')
                        if intervals is None:
                            with open(self.pipe_config.ref['refindex'], 'r', encoding='utf-8') as f_in:
                                intervals = [_.split()[0] for _ in f_in.read().splitlines()]
                        intervals = split_interval(self.pipe_config.ref['refindex'], intervals, 50_000_000)
                        # each interval generate a result file
                        for _part, _region in enumerate(intervals):
                            _sample_ins = SampleIO(_part)
                            _sample_ins.vcf = step_io.outdir / f"part{_part:03d}.vcf.gz"
                            _sample_ins.region = _region
                            step_io.outdct[_part] = _sample_ins
                    if sub_step == 'merge_seq':
                        step_io.outdir = self.pipe_config.outdir
                        step_io.logdir = self.pipe_config.logdir                        
                        step_io.outdct['vcf'] = step_io.outdir / f"{order:02d}_result.vcf.gz"
                    pp_out_dir(step_io.outdir)
                    pp_out_dir(step_io.logdir)
                    step_io_lst.append(step_io)
                    order += 1
        return step_io_lst


    def pp_in(self, step='QC'):
        step_io = StepIO('meta')
        if step in ['QC', 'align']:
            sample_dct = self._pp_meta1(step)
        else:
        # TODO This may raise ERROR for now, because _pp_meta2 accept file type as
        # arguments for now, but pp_in was designed to accept step name.
            sample_dct = self._pp_meta2()
        step_io.outdct = sample_dct
        return step_io


    def pp_out(self, unit, order):
        """Prepare output files

        Args:
            unit (_type_): _description_
            order (_type_): _description_

        Raises:
            ValueError: The unit is not in present unit list
        """
        match unit:
            case 'QC':
                return self.pp_qc(order)
            case 'align':
                return self.pp_align(order)
            case 'call':
                return self.pp_call(order)
            case 'quantify':
                return self.pp_quan(order)
            case _:
                raise ValueError(f'{unit} is illegal')


class PP:
    def __init__(self, config: ConfigPipe, input_io: StepIO,  output_io: StepIO):
        """_summary_

        Args:
            config (ConfigPipe):
                ConfigPipe.project = 'project name'
                ConfigPipe.task = 'RNA/WGS'
                ConfigPipe.units = ['QC', 'align', 'call']
                ConfigPipe.otherInfo = {step1: {Info1: ..., Info2:...}, ...} # Key-values other than 'software', 'bin', and 'params' in config file
                ConfigPipe.software = {step1: software1, step2: software2}
                ConfigPipe.softwareBin = {software1: bin1, software2: bin2}
                ConfigPipe.softwarePara = {software1: {para1: ..., para2: ...},
                                        software2: {para1: ..., para2: ...}}
                ConfigPipe.outdir = Path(The root outdir directory)
                ConfigPipe.logdir = Path(The root logdir directory)
                ConfigPipe.tmp = Path(The root tmp directory)
                ConfigPipe.ref =  {software1: Ref1, software2: Ref2, ...}
            meta (ConfigIO):

        """
        self.pipe_config = config
        self.inIO = input_io
        self.outIO = output_io
