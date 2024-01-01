# -*- coding: utf-8 -*-
# @Time : 2021/8/27 0:25
# @Author : Zhongyi Hua
# @FileName: interface.py
# @Usage: interface between every step to main
# @Note:
# @E-mail: njbxhzy@hotmail.com
from pathlib import Path
import yaml

from Utilities.index import IndexPipe
from Utilities.qc import PPQC
from Utilities.align import PPAlign
from Utilities.calling import PPcall
from Utilities.quantify import PPquantify
from Utilities.Pipe import ConfigPipe, MetaPipe, pp_out_dir


def interface(args):
    """The total interface

    Args:
        command (str): WGS/RNA
        units (list): ['QC', 'align', 'variant', 'merge', 'genotype']
        sysconf
    """
    sysconf = yaml.safe_load(Path(args.sysconf).read_bytes())
    pp_out_dir(args.outdir)
    pp_out_dir(args.logdir)

    ConfigGlobal = ConfigPipe(sysconf, args.outdir, args.logdir,
                              args.project, args.units, args.command)
    # index
    InsIndex = IndexPipe(ConfigGlobal)
    InsIndex.pp_index()
    MetaGlobal = MetaPipe(args.meta, ConfigGlobal)
    # The first step input is from the file and output the files used as next step input
    initial_step_io = MetaGlobal.pp_meta(ConfigGlobal.units[0])
    last_step_io = initial_step_io

    # Other units
    class_dct = {'QC': PPQC, 'align': PPAlign, 'call': PPcall, 'quantify': PPquantify}
    unit_order = 1
    for _unit in ConfigGlobal.units:
        step_io_lst = MetaGlobal.pp_out(_unit, unit_order)
        for sub_step_io in step_io_lst:
            # Some units have sub steps
            # Each PP* class should handle substeps based on 'step_name' them self
            cur_step_io = sub_step_io
            step_ins = class_dct[_unit](ConfigGlobal, last_step_io, cur_step_io)
            step_ins.pp_out()
            last_step_io = cur_step_io
        unit_order += len(step_io_lst)
