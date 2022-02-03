# -*- coding: utf-8 -*-
# @Time    : 2021/8/26 19:25
# @Author  : Zhongyi Hua
# @File    : variantCalling.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import os
from Utilities.utilities import PP


class PPcall(PP):
    def _get_gatk_gvcf_cmd(self, _sample):
        """
        :param _sample: {'sample':sample_name,'bam':bam_path}
        :return:
        """
        ref_index = self.ainone.get('ref')
        gatk_bin = self.ainone['config']['variant']['bin']
        _result_gvcf = os.path.join(self.ainone['outdir'], _sample['sample'] + '.g.vcf.gz')
        _log = os.path.join(self.ainone['logdir'], _sample['sample'] + '.gatkhc.log')
        java_options = "--java-options \"%s\"" % " ".join(self.ainone['config']['variant']['gatk_java']['HaplotypeCaller'])

        per_sample = (f"{gatk_bin} {java_options} HaplotypeCaller -ERC GVCF "
                      f" -R {ref_index} --input {_sample['bam']} --output {_result_gvcf} > {_log}")
        return per_sample

    def _get_samtools_gvcf_cmd(self, _sample):
        pass

    def pp_call(self):
        _cmds = []
        for _sample in self.meta:
            _cmds.append(self._get_gatk_gvcf_cmd(_sample))
        with open(self.ainone.get('outsh'), 'w') as f_out:
            f_out.write('\n'.join(_cmds))


class PPmerge(PP):
    def pp_CombineGVCFs(self):
        gatk_bin = self.ainone['config']['variant']['bin']
        ref_index = self.ainone.get('ref')
        _out_put = self.ainone.get('outfile')
        java_options = "--java-options \"%s\"" % " ".join(self.ainone['config']['variant']['gatk_java']['CombineGVCFs'])

        _cmds = f'{gatk_bin} {java_options} CombineGVCFs --reference {ref_index} --output {_out_put} \\\n'
        _cmds += ''.join('-V %s \\\n' % _sample['gvcf'] for _sample in self.meta)
        _cmds = _cmds.strip(' \\\n')
        with open(self.ainone.get('outsh'), 'w') as f_out:
            f_out.write(_cmds)

    def pp_concat(self):
        # for samtools pipeline
        pass


class PPgenotype(PP):
    def _get_intervals(self):
        with open(self.ainone.get('ref') + '.fai') as f_in:
            all_intervals = [_.split()[0] for _ in f_in.read().splitlines()]
        self.ainone['config']['intervals'] = all_intervals

    def _get_gatk_genotype_cmd(self, _interval):
        ref_index = self.ainone.get('ref')
        gatk_bin = self.ainone['config']['variant']['bin']
        _result_vcf = os.path.join(self.ainone['outdir'], _interval + '.vcf.gz')
        _log = os.path.join(self.ainone['logdir'], _interval + '.gatkgenotype.log')
        java_options = "--java-options \"%s\"" % " ".join(self.ainone['config']['variant']['gatk_java']['GenotypeGVCFs'])

        per_chr = (f"{gatk_bin} {java_options} GenotypeGVCFs"
                   f" -R {ref_index} -V {self.meta} -O {_result_vcf} > {_log}")
        return per_chr

    def _get_merge(self):
        gatk_bin = self.ainone['config']['variant']['bin']
        out_file = os.path.join(self.ainone.get("outdir2"), 'result.vcf.gz')
        _cmd = f'{gatk_bin} MergeVcfs -O {out_file} \\\n'
        for _interval in self.ainone['config']['intervals']:
            vcf_file = os.path.join(self.ainone['outdir'], _interval + '.vcf.gz')
            _cmd += f'-I {vcf_file} \\\n'
        return _cmd

    def pp_geno(self):
        if self.ainone['config'].get('intervals') is None:
            self._get_intervals()
        _cmds = []
        for _interval in self.ainone['config']['intervals']:
            _cmds.append(self._get_gatk_genotype_cmd(_interval))
        with open(self.ainone.get('outsh'), 'w') as f_out:
            f_out.write('\n'.join(_cmds))
        _cmd = self._get_merge()
        with open(self.ainone.get('outsh2'), 'w') as f_out:
            f_out.write(_cmd)
