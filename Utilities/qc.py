# -*- coding: utf-8 -*-
# @Time    : 2021/8/26 15:02
# @Author  : Zhongyi Hua
# @File    : qc.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com
from Utilities.Pipe import PP


class PPQC(PP):
    def _get_fastp_cmd(self):
        """Perform QC of fastq input files, generating fastq
        :param _sample{sample:str, RG:str, fq1:str, fq2:str}
        """
        _bin = self.pipe_config.softwareBin['fastp']
        params = ""
        if params_dct := self.pipe_config.softwarePara['fastp']:
            _thread = params_dct.get('threads', 1)
            params = f"-w {_thread} "
            if other_params := params_dct.get('OtherParams'):
                params += ' '.join(other_params)
        cmds = []
        for _sample_name in self.outIO.outdct:
            _infq1 = self.inIO.outdct.get(_sample_name).fq1
            _infq2 = self.inIO.outdct.get(_sample_name).fq2
            _outfq1 = self.outIO.outdct.get(_sample_name).cleanfq1
            _outfq2 = self.outIO.outdct.get(_sample_name).cleanfq2
            _log1 = self.outIO.logdir / f"{_sample_name}_report.html"
            _log2 = self.outIO.logdir / f"{_sample_name}_fastp_report.txt"
            per_sample = (f"{_bin} {params} -i {_infq1} -o {_outfq1} -I {_infq2} -O {_outfq2} -h {_log1} > {_log2}")
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')

    def pp_out(self):
        self._get_fastp_cmd()
