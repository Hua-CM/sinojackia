# -*- coding: utf-8 -*-
# @Time    : 2021/8/25 20:52
# @Author  : Zhongyi Hua
# @File    : align.py
# @Usage   : NGS alignments with BWA
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com
from Utilities.Pipe import PP


class PPAlign(PP):
    def __wgs_internal(self, software, cmd_pattern):
        """To streamline the code, make an internal function
        """
        _ref = self.pipe_config.ref[software]
        _bin = self.pipe_config.softwareBin[software]
        _samblaster = self.pipe_config.softwareBin['samblaster']
        _sambamba = self.pipe_config.softwareBin['sambamba']
        _params = ''
        if _threads := self.pipe_config.softwarePara[software].get('threads'):
            _params += f'-t {_threads}'
        if _param := self.pipe_config.softwarePara[software].get('OtherParams'):
            _params += ' '.join(_param)
        _tmpdir = self.pipe_config.tmp
        cmds = []
        for _sample_name in self.outIO.outdct:
            _infq1 = self.inIO.outdct.get(_sample_name).cleanfq1
            _infq2 = self.inIO.outdct.get(_sample_name).cleanfq2
            _rg = self.inIO.outdct.get(_sample_name).RG
            _out = self.outIO.outdct.get(_sample_name).bam
            _log1 = self.outIO.logdir / f"{_sample_name}.bwa.log"
            _log2 = self.outIO.logdir / f"{_sample_name}.rmdup.log"
            sample_cmd = cmd_pattern.format(bin=_bin, ref=_ref, params=_params, rg=_rg,
                                            infq1=_infq1, infq2=_infq2, log1=_log1)
            per_sample = (f"{sample_cmd} | "
                          f"{_samblaster} --acceptDupMarks --addMateTags 2> {_log2} | "
                          f"{_sambamba} view -S -f bam -l 0 /dev/stdin | "
                          f"{_sambamba} sort -t {_threads} --tmpdir {_tmpdir} -o {_out} /dev/stdin")
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')

    def _pp_bwa(self):
        self.__wgs_internal('bwa', "{bin} mem -v 2 {params} -R {rg} {ref} {infq1} {infq2} 2> {log1}")
    
    def _pp_bowtie(self):
        self.__wgs_internal('bowtie', "{bin} -S {params} --sam-RG {rg} -x {ref} -1 {infq1} -2 {infq2} 2> {log1}")
    
    def _pp_bowtie2(self):
        self.__wgs_internal('bowtie2', "{bin} {params} --rg {rg} -x {ref} -1 {infq1} -2 {infq2} 2> {log1}")

    def _pp_hisat2(self):
        """
        :param _sample: {sample:str, RG:str, fq1:str, fq2:str, lane:str}
        :return:
        """
        _ref = self.pipe_config.ref['hisat2']
        _bin = self.pipe_config.softwareBin['hisat2']
        params = ''
        if _threads := self.pipe_config.softwarePara['hisat2'].get('threads'):
            params += f'-t {_threads}'
        if _param := self.pipe_config.softwarePara['hisat2'].get('OtherParams'):
            params += ' '.join(_param)
        cmds = []
        for _sample_name in self.outIO.outdct:
            _infq1 = self.inIO.outdct.get(_sample_name).cleanfq1
            _infq2 = self.inIO.outdct.get(_sample_name).cleanfq2
            _rg = self.inIO.outdct.get(_sample_name).RG
            _out = self.outIO.outdct.get(_sample_name).bam
            _out = self.outIO.outdct.get(_sample_name).bam
            _log = self.outIO.logdir / (_sample_name + '.summary.txt')
            per_sample = (f"{_bin} -p {_threads} --summary-file {_log} "
                          f"--dta --rg-id {_rg} -x {_ref} -1 {_infq1} -2 {_infq2} | "
                          f"samtools sort -@ {_threads} > {_out}")
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')

    def _pp_star(self):
        """_summary_
        STAR not support other params for now.
        """

        params = f"--genomeLoad NoSharedMemory --readFilesCommand zcat  --quantMode TranscriptomeSAM " \
                 f"--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD " \
                 f"--outSAMunmapped Within --outSAMheaderHD @HD VN:1.4 SO:unsorted --outFilterType BySJout " \
                 f"--outFilterMultimapNmax 20 --outFilterMismatchNmax 999 " \
                 f"--outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 " \
                 f"--alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1"
        _ref = self.pipe_config.ref['STAR']
        _bin = self.pipe_config.softwareBin['STAR']
        if params_dct := self.pipe_config.softwarePara.get('STAR'):
            _thread = params_dct.get('threads', 1)
            params += f" --runThreadN {_thread}"
        cmds = []
        for _sample_name in self.outIO.outdct:
            _log = self.outIO.logdir / f"{_sample_name}.log"
            _infq1 = self.inIO.outdct.get(_sample_name).cleanfq1
            _infq2 = self.inIO.outdct.get(_sample_name).cleanfq2
            _out = self.outIO.outdct.get(_sample_name).bam # Take it easy, just use 'bam' attr here, it is a prefix
            per_sample = (f"{_bin} --genomeDir {_ref} --readFilesIn {_infq1} {_infq2} "
                          f"--outFileNamePrefix {_out} {params} > {_log}")
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')


    def pp_out(self):
        for _software in self.pipe_config.software.values():
            match _software:
                case 'hisat2':
                    self._pp_hisat2()
                case 'STAR':
                    self._pp_star()
                case 'bwa':
                    self._pp_bwa()
                case 'bowtie':
                    self._pp_bowtie()
                case 'bowtie2':
                    self._pp_bowtie2()
                case _:
                    continue
