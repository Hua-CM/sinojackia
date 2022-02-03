# -*- coding: utf-8 -*-
# @Time : 2021/8/27 0:25
# @Author : Zhongyi Hua
# @FileName: WGS.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from Utilities.index import PPIndex
from Utilities.qc import PPQC
from Utilities.align import PPAlgin
from Utilities.variantCalling import PPcall, PPmerge, PPgenotype
from Utilities.utilities import GenConfig, GenMeta


def wgs(args):
    start_dict = {'QC': 0, 'align': 1, 'variant': 2, 'merge': 3, 'genotype': 4}
    start_point = 5
    for unit in args.units:
        start_point = min(start_point, start_dict[unit])
    start_point = list(start_dict.keys())[start_point]
    InsConfig = GenConfig(args.sysconf, args.outdir, args.logdir, args.project, args.units, 'WGS')
    InsMeta = GenMeta(args.meta, args.outdir, args.project, 'WGS', start_point=start_point)
    # index
    InsIndex = PPIndex(InsConfig.pp_index())
    InsIndex.copy_file()
    InsIndex.pp_index()
    if 'QC' in args.units:
        PPQC(InsConfig.pp_qc(), InsMeta.pp_qc()).pp_qc()
    if 'align' in args.units:
        PPAlgin(InsConfig.pp_align(), InsMeta.pp_align()).pp_align()
    if 'variant' in args.units:
        PPcall(InsConfig.pp_call(), InsMeta.pp_quan_var()).pp_call()
    if 'merge' in args.units:
        PPmerge(InsConfig.pp_merge(), InsMeta.pp_merge()).pp_CombineGVCFs()
    if 'genotype' in args.units:
        PPgenotype(InsConfig.pp_genotype(), InsMeta.pp_genotype()).pp_geno()


def rna(args):
    pass
