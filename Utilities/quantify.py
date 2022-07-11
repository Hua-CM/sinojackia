# -*- coding: utf-8 -*-
# @Time : 2021/9/29 0:42
# @Author : Zhongyi Hua
# @FileName: quantify.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com


import os
from Utilities.utilities import PP


class PPquantify(PP):
    def _get_rsem_cmd(self, _sample):
        _ref_index = self.ainone.get('ref')
        _rsem = os.path.join(self.ainone['config']['quantify']['bin'], 'rsem-calculate-expression')
        _out = os.path.join(self.ainone['outdir'], _sample['sample'])
        _log = os.path.join(self.ainone['logdir'], _sample['sample'] + '.log')
        per_sample = f'{_rsem}  --paired-end --alignments {_sample["bam"]} {_ref_index} {_out}'
        return per_sample

    def _get_featurecount_cmd(self, _sample):
        pass

    def pp_quantify(self):
        _cmds = []
        if self.ainone['config']['quantify']['software'] == 'RSEM':
            _get_func = self._get_rsem_cmd
        elif self.ainone['config']['quantify']['software'] == 'featureCount':
            _get_func = self._get_featurecount_cmd
        for _sample in self.meta:
            _cmds.append(_get_func(_sample))
        with open(self.ainone.get('outsh'), 'w') as f_out:
            f_out.write('\n'.join(_cmds))
            f_out.write('\n')
