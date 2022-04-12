# -*- coding: utf-8 -*-
# @Time    : 2021/8/26 11:35
# @Author  : Zhongyi Hua
# @File    : utilities.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import os
import sys
import yaml


def pp_out_dir(_output_dir, _sub_dir=''):
    if not os.path.exists(_output_dir):
        os.mkdir(_output_dir)
    else:
        if os.path.isfile(_output_dir):
            print("The output directory seems a file?")
            sys.exit(1)
    # make some output directory
    try:
        os.mkdir(os.path.join(_output_dir, _sub_dir))
    except FileExistsError:
        pass


class GenConfig:
    def __init__(self, config, outdir, logdir, project, units, task):
        self.config = self._pp_config(config, outdir)
        self.outdir = outdir
        pp_out_dir(outdir)
        self.logdir = logdir
        self._ainone_ = {'project': project, 'units': units, 'task': task}

    @staticmethod
    def _pp_config(_file, outdir):
        with open(_file) as C:
            _config = yaml.safe_load(C)
        # set temporary directory
        if _config.get('resources').get('temporary') is None:
            _config['resources']['temporary'] = os.path.join(outdir, 'tmp')
            pp_out_dir(outdir, 'tmp')
        else:
            pp_out_dir(_config['resources']['temporary'])
        return _config

    def pp_index(self):
        self._ainone_.update(
            {'ref': os.path.join(self.outdir,
                                 '00_index',
                                 os.path.split(self.config['resources']['reference'])[1]
                                 )}
        )
        ainone = dict()
        ainone.update(self._ainone_)
        pp_out_dir(self.outdir, '00_index')
        ainone.setdefault('outsh', os.path.join(self.outdir, self._ainone_.get('project') + '.step0.index.sh'))
        ainone.setdefault('outdir', os.path.join(self.outdir, '00_index'))
        ainone.setdefault('config', self.config)
        return ainone

    def pp_qc(self):
        ainone = dict()
        ainone.update(self._ainone_)
        pp_out_dir(self.outdir, '01_cleandata')
        pp_out_dir(self.logdir, '01_QC')
        ainone.setdefault('outsh', os.path.join(self.outdir, self._ainone_.get('project') + '.step1.QC.sh'))
        ainone.setdefault('logdir', os.path.join(self.logdir, '01_QC'))
        ainone.setdefault('outdir', os.path.join(self.outdir, '01_cleandata'))
        ainone.setdefault('config', self.config)
        return ainone

    def pp_align(self):
        ainone = dict()
        ainone.update(self._ainone_)
        ainone.setdefault('outsh', os.path.join(self.outdir, self._ainone_.get('project') + '.step2.align.sh'))
        pp_out_dir(self.outdir, '02_mapping')
        ainone.setdefault('outdir', os.path.join(self.outdir, '02_mapping'))
        pp_out_dir(self.logdir, '02_mapping')
        ainone.setdefault('logdir', os.path.join(self.logdir, '02_mapping'))
        ainone.setdefault('config', self.config)
        return ainone

    def pp_quan(self):
        ainone = dict()
        ainone.update(self._ainone_)
        ainone.setdefault('outsh', os.path.join(self.outdir, self._ainone_.get('project') + '.step3.quan.sh'))
        pp_out_dir(self.outdir, '03_expression')
        ainone.setdefault('outdir', os.path.join(self.outdir, '03_expression'))
        pp_out_dir(self.logdir, '03_expression')
        ainone.setdefault('logdir', os.path.join(self.logdir, '03_expression'))
        ainone.setdefault('config', self.config)
        return ainone

    def pp_call(self):
        ainone = dict()
        ainone.update(self._ainone_)
        ainone.setdefault('outsh', os.path.join(self.outdir, self._ainone_.get('project') + '.step3.calling.sh'))
        pp_out_dir(self.outdir, '03_gvcf')
        ainone.setdefault('outdir', os.path.join(self.outdir, '03_gvcf'))
        pp_out_dir(self.logdir, '03_gvcf')
        ainone.setdefault('logdir', os.path.join(self.logdir, '03_gvcf'))
        ainone.setdefault('config', self.config)
        return ainone

    def pp_merge(self):
        ainone = dict()
        ainone.update(self._ainone_)
        ainone.setdefault('outsh', os.path.join(self.outdir, self._ainone_.get('project') + '.step4.merge.sh'))
        ainone.setdefault('outfile', os.path.join(self.outdir, '04_genomeDB.g.vcf.gz'))
        ainone.setdefault('logfile', os.path.join(self.logdir, '04_genomeDB.log'))
        ainone.setdefault('config', self.config)
        return ainone

    def pp_genotype(self):
        ainone = dict()
        ainone.update(self._ainone_)
        ainone.setdefault('outsh', os.path.join(self.outdir, self._ainone_.get('project') + '.step5.genotype.sh'))
        ainone.setdefault('outsh2', os.path.join(self.outdir, self._ainone_.get('project') + '.step6.merge.sh'))
        pp_out_dir(self.outdir, '05_vcf')
        ainone.setdefault('outdir', os.path.join(self.outdir, '05_vcf'))
        ainone.setdefault('outdir2', self.outdir)
        pp_out_dir(self.logdir, '05_vcf')
        ainone.setdefault('logdir', os.path.join(self.logdir, '05_vcf'))
        ainone.setdefault('config', self.config)
        return ainone


class GenMeta:
    def __init__(self, meta, outdir, project, task, start_point='QC'):
        self.meta = meta
        self.outdir = outdir
        pp_out_dir(outdir)
        self.project = project
        self.sp = start_point
        self.task = task

    @staticmethod
    def _pp_meta1(_file):
        """
        :return [{sample:str, RG:str, fq1:str, fq2:str, lane:str},...]
        """
        _samplelist = []
        with open(_file, 'r') as f_in:
            for _ in f_in.read().splitlines():
                if not _.startswith('#'):
                    _sample = dict()
                    _sample['sample'], _sample['RG'], _sample['fq1'], _sample['fq2'] = _.split()
                    _samplelist.append(_sample)
        return _samplelist

    @staticmethod
    def _pp_meta2(_file, _name='bam'):
        _samplelist = []
        with open(_file, 'r') as f_in:
            for _ in f_in.read().splitlines():
                _samplelist.append(dict(zip(('sample', _name), _.split('\t'))))
        return _samplelist

    def pp_qc(self):
        _samplelst = self._pp_meta1(self.meta)
        return _samplelst

    def pp_align(self):
        if self.sp == 'align':
            _samplelst = self._pp_meta1(self.meta)
        else:
            in_dir = os.path.join(self.outdir, '01_cleandata')
            _samplelst = self.pp_qc()
            for _sample in _samplelst:
                _sample['fq1'] = f"{in_dir}/{_sample['sample']}_1_clean.fq.gz"
                _sample['fq2'] = f"{in_dir}/{_sample['sample']}_2_clean.fq.gz"
        return _samplelst

    def pp_quan_var(self):
        if self.sp == 'quantity':
            _samplelst = self._pp_meta2(self.meta)
        else:
            in_dir = os.path.join(self.outdir, '02_mapping')
            _samplelst_tmp = self.pp_align()
            _samplelst = []
            for _sample in _samplelst_tmp:
                _tmp = dict()
                _tmp.setdefault('sample', _sample['sample'])
                if self.task == 'RNA':
                    _tmp.setdefault('bam', f"{in_dir}/{_sample['sample']}Aligned.toTranscriptome.out.bam")
                else:
                    _tmp.setdefault('bam', f"{in_dir}/{_sample['sample']}.sorted.bam")
                _samplelst.append(_tmp)
        return _samplelst

    def pp_merge(self):
        if self.sp == 'merge':
            _samplelst = self._pp_meta2(self.meta, 'gvcf')
        else:
            in_dir = os.path.join(self.outdir, '03_gvcf')
            _samplelst_tmp = self.pp_quan_var()
            _samplelst = []
            for _sample in _samplelst_tmp:
                _tmp = dict()
                _tmp.setdefault('sample', _sample['sample'])
                _tmp.setdefault('gvcf', f"{in_dir}/{_sample['sample']}.g.vcf.gz")
                _samplelst.append(_tmp)
        return _samplelst

    def pp_genotype(self):
        """
        Interval information was in config
        """
        if self.sp == 'genotype':
            os.symlink(self.meta, os.path.join(self.outdir, os.path.splitext(self.meta)[1]))
        return os.path.join(self.outdir, '04_genomeDB.g.vcf.gz')


class PP:
    def __init__(self, ainone, meta):
        self.ainone = ainone
        self.meta = meta
