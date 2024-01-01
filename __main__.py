# -*- coding: utf-8 -*-
# @Time    : 2021/8/25 20:05
# @Author  : Zhongyi Hua
# @File    : main.py
# @Usage   :
# @Note    :
# @E-mail  : njbxhzy@hotmail.com

import sys
import argparse
from interface import interface
from pathlib import Path


def getArgs():
    """Parse input commandline arguments, handling multiple cases.
    """
    cmdparser = argparse.ArgumentParser(description="sinajackia: A HTS analysis pipeline generator for plants.")
    commands = cmdparser.add_subparsers(dest="command", title="sinojackia commands")

    # The standard pipeline for WGS.
    wgs_cmd = commands.add_parser("WGS", help="Creating pipeline for WGS(from fastq to genotype VCF)")
    wgs_cmd.add_argument("-c", "--conf", dest="sysconf", required=True, type=Path,
                         help="<File path> YAML configuration file specifying details.")
    wgs_cmd.add_argument("-m", "--meta", dest="meta", required=True, type=Path,
                         help="<File path> File information table (see README/example).")
    wgs_cmd.add_argument("-o", "--out", dest="outdir", required=True,  type=Path,
                         help="<Directory Path> A directory for output results.")
    wgs_cmd.add_argument("-l", "--log", dest="logdir", required=True, type=Path,
                         help="<Directory Path> A directory for output results.")
    wgs_cmd.add_argument("-p", "--project", dest="project", type=str, default="WGS",
                         help="Name of the project. Default value: test")
    wgs_cmd.add_argument("-u", "--units", dest="units", type=str, default='all',
                         help="Specific one or more unit(s) (separated by comma) of WGS processes. "
                              "Possible values: {QC, align, call, all}"
                              "Defualt value: all"
                         )

    # The standard pipeline for RNA-seq
    rna_cmd = commands.add_parser("RNA", help="Creating pipeline for RNA-seq (from fastq to expression matrix)")
    rna_cmd.add_argument("-c", "--conf", dest="sysconf", required=True, type=Path,
                         help="<File path> YAML configuration file specifying details.")
    rna_cmd.add_argument("-m", "--meta", dest="meta", required=True, type=Path,
                         help="<File path> File information table (see README/example).")
    rna_cmd.add_argument("-o", "--out", dest="outdir", required=True, type=Path,
                         help="<Directory Path> A directory for output results.")
    rna_cmd.add_argument("-l", "--log", dest="logdir", required=True, type=Path,
                         help="<Directory Path> A directory for output results.")
    rna_cmd.add_argument("-p", "--project", dest="project", type=str, default="RNA",
                         help="Name of the project. Default value: test")
    rna_cmd.add_argument("-u", "--units", dest="units", type=str, default='all',
                         help="Specific one or more unit(s) (separated by comma) of WGS processes. "
                              "Possible values: {QC, align, quantify, all}. Defualt value:all.",
                         )
    return cmdparser.parse_args()


def main():
    kwargs = getArgs()
    # change to absolute path
    if kwargs.command is None:
        print("Please type: sinojackia -h or sinojackia --help to show the help message.")
        sys.exit(1)
    else:
        kwargs.logdir = kwargs.logdir.resolve()
        kwargs.outdir = kwargs.outdir.resolve()
        kwargs.units = kwargs.units.split(',')
        interface(kwargs)


if __name__ == '__main__':
    main()
