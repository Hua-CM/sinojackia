# -*- coding: utf-8 -*-
# @Time : 2021/9/28 1:27
# @Author : Zhongyi Hua
# @FileName: index.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com
import logging
import subprocess as sp
from Utilities.Pipe import ConfigPipe


logging.basicConfig(
    level='INFO'
)
logger = logging.getLogger('')


class IndexPipe:
    def __init__(self, pipe_config: ConfigPipe):
        """_summary_

        Args:
            pipe_config (ConfigPipe):
                pipe_config.ref: All possible index file output path
        """
        self.pipe_config = pipe_config

    def _copy_file(self):
        try:
            self.pipe_config.ref['ref'].symlink_to(self.pipe_config.ref['raw_ref'])
            if gff_path := self.pipe_config.ref.get('gff'):
                gff_path.symlink_to(self.pipe_config.ref['raw_gff'])
            if gff_path := self.pipe_config.ref.get('gtf'):
                gff_path.symlink_to(self.pipe_config.ref['raw_gtf'])
        except FileExistsError:
            pass

    def _index_fasta(self):
        logger.info('Index your reference genome. It may take a few minutes')
        _bin = self.pipe_config.softwareBin['samtools']
        _out = self.pipe_config.ref['ref']
        sp.run([_bin, 'faidx', _out], check=False)

    def _write_cmds(self, cmds):
        with open(self.pipe_config.ref['outsh'], 'a', encoding='utf-8') as fp:
            fp.write(cmds+'\n')

    def _pp_star_index(self):
        _bin = self.pipe_config.softwareBin['STAR']
        _out = self.pipe_config.ref['STAR']
        _params = ''
        if params_dct := self.pipe_config.softwarePara.get('STAR'):
            _thread = params_dct.get('threads', 1)
            _params = f"--runThreadN {_thread}"
        _index_cmd = (f"{_bin} {_params} "
                      f"--runMode genomeGenerate --sjdbGTFtagExonParentTranscript Parent "
                      f"--genomeDir {_out} "
                      f"--genomeFastaFiles {self.pipe_config.ref['ref']} "
                      f"--sjdbGTFfile {self.pipe_config.ref['gff']}")
        self._write_cmds(_index_cmd)

    def _pp_hisat_index(self):
        _bin = self.pipe_config.softwareBin['hisat2-build']
        _out = self.pipe_config.ref['hisat2']
        _params = ""
        if params_dct := self.pipe_config.softwarePara.get('hisat2'):
            _thread = params_dct.get('threads', 1)
            _params = f"--threads {_thread}"
        _index_cmd = f"{_bin} {_params} {self.pipe_config.ref['ref']} {_out}"
        self._write_cmds(_index_cmd)

    def _pp_rsem_index(self):
        _bin = self.pipe_config.softwareBin['RSEM'] / 'rsem-prepare-reference'
        _out = self.pipe_config.ref['RSEM']
        _index_cmd = f"{_bin} --gff3 {self.pipe_config.ref['gff']} {self.pipe_config.ref['ref']} {_out}"
        self._write_cmds(_index_cmd)
        
    def _pp_bwa_index(self):
        _bin = self.pipe_config.softwareBin['bwa']
        _out = self.pipe_config.ref['bwa']
        _index_cmd = f"{_bin} index {_out}"
        self._write_cmds(_index_cmd)
    
    def _pp_bowtie_index(self):
        _bin = self.pipe_config.softwareBin['bowtie-build']
        _in = self.pipe_config.ref['ref']
        _out = self.pipe_config.ref['bowtie']
        _index_cmd = f"{_bin} {_in} {_out}"
        self._write_cmds(_index_cmd)
    
    def _pp_bowtie2_index(self):
        _bin = self.pipe_config.softwareBin['bowtie2-build']
        _in = self.pipe_config.ref['ref']
        _out = self.pipe_config.ref['bowtie2']
        _index_cmd = f"{_bin} {_in} {_out}"
        self._write_cmds(_index_cmd)

    def _pp_gatk_index(self):
        _bin = self.pipe_config.softwareBin['GATK']
        _out = self.pipe_config.ref['GATK']
        _index_cmd = f"{_bin} -Xmx10G CreateSequenceDictionary -R {_out}"
        self._write_cmds(_index_cmd)

    def pp_index(self):
        self._copy_file()
        # Only WGS task need index fasta
        if  self.pipe_config.task == 'WGS':
            self._index_fasta()
        for _software in self.pipe_config.software.values():
            match _software:
                case 'hisat2':
                    self._pp_hisat_index()
                case 'STAR':
                    self._pp_star_index()
                case 'RSEM':
                    self._pp_rsem_index()
                case 'bwa':
                    self._pp_bwa_index()
                case 'bowtie':
                    self._pp_bowtie_index()
                case 'bowtie2':
                    self._pp_bowtie2_index()
                case 'GATK':
                    self._pp_gatk_index()
                case _:
                    # For now, indexing is not required for these software applications:
                    #   htseq-count
                    continue
