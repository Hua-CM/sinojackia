# -*- coding: utf-8 -*-
# @Time    : 2021/8/25 20:52
# @Author  : Zhongyi Hua
# @File    : align.py
# @Usage   : NGS alignments with BWA
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import os
from Utilities.utilities import PP


class PPAlgin(PP):
    def _get_hisat2_cmd(self, _sample):
        """
        :param _sample: {sample:str, RG:str, fq1:str, fq2:str, lane:str}
        :return:
        """
        _bin = os.path.join(self.ainone['config']['align']['bin'], "hisat2-align")
        threads = self.ainone['config']['align']['threads']
        _res_bam = os.path.join(self.ainone['outdir'], f"{_sample['sample']}Aligned.toTranscriptome.out.bam")
        per_sample = (f"{_bin} -p {threads} --dta --rg-id {_sample['RG']} -1 {_sample['fq1']} -2 {_sample['fq2']} | "
                      f"samtools sort -@ {threads} > {_res_bam}")
        return per_sample

    def _get_bwa_mem_cmd(self, _sample):
        """Perform piped alignment of fastq input files, generating sorted output BAM
        :param _sample{sample:str, RG:str, fq1:str, fq2:str}
        """
        ref_index = self.ainone.get('ref')
        bwa_bin = self.ainone['config']['align']['bwa']
        _samblaster = self.ainone['config']['align']['samblaster']
        _sambamba = self.ainone['config']['align']['sambamba']
        threads = self.ainone['config']['align']['threads']
        tmpdir = self.ainone['config']['resources']['temporary']
        _res_bam = os.path.join(self.ainone['outdir'], _sample['sample'] + '.sorted.bam')
        _log1 = os.path.join(self.ainone['logdir'], _sample['sample'] + '.bwa.log')
        _log2 = os.path.join(self.ainone['logdir'], _sample['sample'] + '.marked.log')
        per_sample = (f"{bwa_bin} mem -v 2 -t {threads} -R {_sample['RG']} {ref_index} {_sample['fq1']} {_sample['fq2']} 2> {_log1} |"
                      f"{_samblaster} --acceptDupMarks --addMateTags 2> {_log2} |"
                      f"{_sambamba} view -S -f bam -l 0 /dev/stdin |"
                      f"{_sambamba} sort -t {threads} --tmpdir {tmpdir} -o {_res_bam} /dev/stdin")
        return per_sample

    def _get_star_mem_cmd(self, _sample):
        ref_dir = os.path.join(self.ainone['outdir'], '../00_index')
        star_bin = self.ainone['config']['align']['bin']
        threads = self.ainone['config']['align']['threads']
        _log1 = os.path.join(self.ainone['logdir'], _sample['sample'] + '.star.log')
        _result_prefix = os.path.join(self.ainone['outdir'], _sample['sample'])
        per_sample = (f"{star_bin} --runThreadN {threads} --genomeDir {ref_dir} "
                      f"--readFilesIn {_sample['fq1']} {_sample['fq2']} "
                      f"--outFileNamePrefix {_result_prefix} --genomeLoad LoadAndRemove "
                      f"--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM > {_log1}")
        return per_sample

    def pp_align(self):
        """
        :param ainone: a dict prepared by kwargs
        :param outsh: the path for prepared shell script
        :return:
        """
        _cmds = []
        if self.ainone.get('task') == 'WGS':
            _get_func = self._get_bwa_mem_cmd
        elif self.ainone.get('task') == 'RNA':
            if self.ainone['config']['align']['software'] == 'STAR':
                _get_func = self._get_star_mem_cmd
            elif self.ainone['config']['align']['software'] == 'hisat2':
                _get_func = self._get_hisat2_cmd
        for _sample in self.meta:
            _cmds.append(_get_func(_sample))
        with open(self.ainone.get('outsh'), 'w') as f_out:
            f_out.write('\n'.join(_cmds))
