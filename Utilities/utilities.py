# -*- coding: utf-8 -*-
# @File    :   utilities.py
# @Time    :   2023/10/27 00:36:42
# @Author  :   Zhongyi Hua 
# @Usage   :   
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com

from collections import defaultdict
from pathlib import Path
from typing import List


def split_interval(fai_path: Path, intervals: List[str] = None, winsize:int = 50_000_000):
    lines = fai_path.read_text().splitlines()
    seq_dct = defaultdict()
    for _line in lines:
        _seq_name, _seq_len = _line.split()[:2]
        seq_dct[_seq_name] = int(_seq_len)
    if intervals:
       seq_dct = {_key: _value for _key, _value in seq_dct.items() if _key in intervals} 
    region_lst = []
    for _seq_name, _seq_len in seq_dct.items():
        beg_pos = 1
        end_pos = winsize
        while end_pos < _seq_len:
            region_lst.append(f"{_seq_name}:{beg_pos}-{end_pos}")
            beg_pos = end_pos + 1
            end_pos += winsize
        region_lst.append(f"{_seq_name}:{beg_pos}-{_seq_len}")    
    return region_lst
