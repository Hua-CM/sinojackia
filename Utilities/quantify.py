# -*- coding: utf-8 -*-
# @Time : 2021/9/29 0:42
# @Author : Zhongyi Hua
# @FileName: quantify.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com


import os


def _parse_args(kwargs):
    """
    metalist: each bam per line
    """
    ainone = dict()
    # load config file
    with open(kwargs.sysconf) as C:
        ainone['config'] = yaml.safe_load(C)
    ainone['outdir'] = kwargs.outdir
    prepare_out_dir(ainone['outdir'], '03_expression')
    ainone['loginfo'] = os.path.join(kwargs.outdir, 'loginfo', '03_expression')
    prepare_out_dir(os.path.join(kwargs.outdir, 'loginfo'), '03_expression')
    if 'align' in kwargs.unit:
        ainone['samplelist'] = align2expression(kwargs.metalist, ainone['config'], ainone['outdir'])
    else:
        with open(kwargs.metalist, 'r') as f_in:
            ainone['samplelist'] = f_in.read().split('\n')
    return ainone


def _get_rsem_expression_cmd(_sample, config, out_dir, log_dir):
    """
     rsem-calculate-expression --paired-end \
                               --alignments \
                               -p 8 \
                               /data/mmliver_paired_end_quals.bam \
                               /ref/mouse_125 \
                               mmliver_paired_end_quals
    """
    pass