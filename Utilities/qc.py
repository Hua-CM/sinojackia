# -*- coding: utf-8 -*-
# @Time    : 2021/8/26 15:02
# @Author  : Zhongyi Hua
# @File    : qc.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com
import os
from Utilities.utilities import PP


class PPQC(PP):
    def _get_fastp_cmd(self, _sample):
        """Perform QC of fastq input files, generating fastq
        :param _sample{sample:str, RG:str, fq1:str, fq2:str}
        """
        fastp_bin = self.ainone['config']['QC']['bin']
        threads = self.ainone['config']['QC']['threads']
        _log1 = os.path.join(self.ainone['logdir'], _sample['sample']+'_report.html')
        _log2 = os.path.join(self.ainone['logdir'], _sample['sample']+'_fastp_report.txt')
        per_sample = (f"{fastp_bin} -q 20 -l 50 -g -x -w {threads} "
                      f"-i {_sample['fq1']} -o {self.ainone['outdir']}/{_sample['sample']}_1_clean.fq.gz "
                      f"-I {_sample['fq2']} -O {self.ainone['outdir']}/{_sample['sample']}_2_clean.fq.gz "
                      f"-h {_log1} > {_log2} ")
        return per_sample

    def pp_qc(self):
        """
        :param ainone: a dict prepared by kwargs
        :param outsh: the path for prepared shell script
        :return:
        """
        _cmds = []
        for _sample in self.meta:
            if self.ainone['config']['QC']['software'] == 'fastp':
                _func1 = self._get_fastp_cmd
            _cmds.append(_func1(_sample))
        with open(self.ainone['outsh'], 'w') as f_out:
            f_out.write('\n'.join(_cmds))
            f_out.write('\n')
