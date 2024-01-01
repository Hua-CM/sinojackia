# -*- coding: utf-8 -*-
# @Time    : 2021/8/26 19:25
# @Author  : Zhongyi Hua
# @File    : variantCalling.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

"""
This is the most complicated module because there are many substeps in 
calling vcf whatever software you used.

For GATK:
1. HaplotypeCaller single sample gvcf
2. Merge samples' gvcf into a big gvcf file
3. Genotyping each chromosone one by one
4. Merge chromosome results

For bcftools
1. Split genomes into mutiple regions (I prefer 50M per file)
2. Genotype each region
3. Merge different regions into one file
"""

from Utilities.Pipe import PP, StepIO


class CallSub:
    def __init__(self, input_io: StepIO, output_io: StepIO, call_bin, call_params, call_ref) -> None:
        self.inIO = input_io
        self.outIO = output_io
        self.bin = call_bin
        self.params = call_params
        self.ref = call_ref

class GVCF(CallSub):
    def pp_out(self):
        cmds = []
        for _sample_name in self.outIO.outdct:
            _in_bam = self.inIO.outdct.get(_sample_name).bam
            _out_gvcf = self.outIO.outdct.get(_sample_name).gvcf
            _log = (self.outIO.logdir / (_sample_name + '.gatkhc.log'))
            if _params := self.params.get('HaplotypeCaller'):
                java_options = "--java-options \"%s\"" % " ".join(_params)
            else:
                java_options = ''
            per_sample = (f"{self.bin} {java_options} HaplotypeCaller -ERC GVCF "
                        f"-R {self.ref} --input {_in_bam} --output {_out_gvcf} > {_log}")
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')

class MergeSample(CallSub):
    def pp_out(self):
        if _params := self.params.get('CombineGVCFs'):
            java_options = "--java-options \"%s\"" % " ".join(_params)
        else:
            java_options = ''
        _out = self.outIO.outdct.get('gvcfDB')
        cmds = f'{self.bin} {java_options} CombineGVCFs --reference {self.ref} --output {_out} \\\n'
        for _sample_name in self.inIO.outdct:
            _in_gvcf = self.inIO.outdct.get(_sample_name).gvcf
            cmds += f'-V {_in_gvcf} \\\n'
        cmds = cmds.strip(' \\\n') + '\n'
        self.outIO.outsh.write_text(cmds)

class GenotypeGATK(CallSub):
    def pp_out(self):
        cmds = []
        if _params := self.params.get('GenotypeGVCFs'):
            java_options = "--java-options \"%s\"" % " ".join(_params)
        else:
            java_options = ''
        for _interval in self.outIO.outdct:
            _out = self.outIO.outdct[_interval].vcf
            _log = self.outIO.logdir / f'{_interval}.log'
            per_chr = (f"{self.bin} {java_options} GenotypeGVCFs " 
                       f"-R {self.ref} -L {_interval} -V {self.inIO.outdct['gvcfDB']} -O {_out} > {_log}")
            cmds.append(per_chr)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')

class MergeSeqGATK(CallSub):
    def pp_out(self):
        if _params := self.params.get('MergeVcfs'):
            java_options = "--java-options \"%s\"" % " ".join(_params)
        else:
            java_options = ''
        _out = self.outIO.outdct['vcf']
        cmd = f'{self.bin} {java_options} MergeVcfs -O {_out} \\\n'
        for _interval_ins in self.inIO.outdct.values():
            vcf_file = _interval_ins.vcf
            cmd += f'-I {vcf_file} \\\n'
        cmd = cmd.strip(' \\\n') + '\n'
        self.outIO.outsh.write_text(cmd)

class GenotypeBcf(CallSub):
    def pp_out(self):
        paras = ""
        if threads := self.params.get('threads'):
            paras = f"--threads {threads}"
        cmds = []
        # generate a bam_lst
        bam_lst = [str(self.inIO.outdct.get(_sample_name).bam) for _sample_name in self.inIO.outdct]
        bam_lst_path = self.outIO.outdir / 'bam.lst'
        bam_lst_path.write_text('\n'.join(bam_lst) + '\n')
        for _part in self.outIO.outdct:
            _out_vcf = self.outIO.outdct.get(_part).vcf
            _region = self.outIO.outdct.get(_part).region
            per_sample = f"{self.bin} mpileup {paras} -r {_region} -Ou -g -f {self.ref} -b {bam_lst_path} | {self.bin} call -mv -Oz -o {_out_vcf}"
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds))

class MergeSeqBcf(CallSub):
    def pp_out(self):
        paras = ""
        if threads := self.params.get('threads'):
            paras = f"--threads {threads}"
        vcf_files = ' '.join([str(_part_ins.vcf) for _part_ins in self.inIO.outdct.values()])
        cmd = f"{self.bin} concat {paras} {vcf_files}"
        self.outIO.outsh.write_text(cmd)

class PPcall(PP):
    def pp_out(self):
        call_software = self.pipe_config.software.get('call')
        call_bin = self.pipe_config.softwareBin.get(call_software)
        call_para = self.pipe_config.softwarePara.get(call_software)
        call_ref = self.pipe_config.ref.get(call_software)
        match call_software:
            case 'GATK':
                class_dct = {'gvcf': GVCF,
                             'merge_sample': MergeSample,
                             'genotype': GenotypeGATK,
                             'merge_seq': MergeSeqGATK}
            case 'bcftools':
                class_dct = {'genotype': GenotypeBcf,
                             'merge_seq':MergeSeqBcf}
        SubClass = class_dct[self.outIO.step_name]
        sub_ins = SubClass(self.inIO, self.outIO, call_bin, call_para, call_ref)
        sub_ins.pp_out()
