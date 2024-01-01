# -*- coding: utf-8 -*-
# @Time : 2021/9/29 0:42
# @Author : Zhongyi Hua
# @FileName: quantify.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com
import inspect
from Utilities.Pipe import PP


class PPquantify(PP):
    def _pp_rsem(self, bam_dct):
        _ref = self.pipe_config.ref['RSEM']
        _bin = self.pipe_config.softwareBin['RSEM'] / 'rsem-calculate-expression'
        cmds = []
        for _sample_name in self.outIO.outdct:
            _in = bam_dct.get(_sample_name)
            _out = self.outIO.outdct.get(_sample_name).tsv
            per_sample = f'{_bin}  --paired-end --alignments {_in} {_ref} {_out}'
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')
    
    def _pp_rsem_merge(self):
        _bin = self.pipe_config.softwareBin['RSEM'] / 'rsem-generate-data-matrix'
        cmd = f"{_bin} {self.inIO.outdir} > {self.outIO.outdct['Alltsv']}"
        self.outIO.outsh.write_text(cmd)

    def _pp_htseqcount(self, bam_dct):
        """
        	1. I recommend use `gffread` to generate gtf file used here.
               A command is `gffread /path/to/gff -E -T -o /path/to/gtf`
	        2. The suffix to expression output file must be standard suffixs, such as  `tsv`, `csv`. 
               The htseq use these suffixes to determine the output file format
        """
        _bin = self.pipe_config.softwareBin['htseq-count']
        _gtf = self.pipe_config.ref['gtf']
        cmds = []
        for _sample_name in self.outIO.outdct:
            _in = bam_dct.get(_sample_name)
            _out = self.outIO.outdct.get(_sample_name).tsv
            per_sample = f'{_bin} --format bam {_in} {_gtf} -c {_out}'
            cmds.append(per_sample)
        self.outIO.outsh.write_text('\n'.join(cmds) + '\n')

    
    def _pp_htseqcount_merge(self):
        python_script = inspect.cleandoc(f"""
            #!/usr/bin/python3
            from pathlib import Path
            from collections import defaultdict
            
            all_dct = defaultdict(list)
            sample_lst = []
            result_path = Path("{self.inIO.outdir}")
            for _sample_path in result_path.glob('*'):
                _sample_path = Path(_sample_path)
                sample_lst.append(_sample_path.stem)
                for _line in _sample_path.read_text().strip().split('\\n'):
                    _gene, _count = _line.split('\\t')
                    all_dct[_gene].append(_count)
            # Generate output
            out_str = '\\t' + '\\t'.join(sample_lst) + '\\n'
            for _gene, _counts in all_dct.items():
                out_str += _gene + '\\t' + '\\t'.join(_counts) + '\\n'
            Path("{self.outIO.outdct['Alltsv']}").write_text(out_str)""")
        self.outIO.outsh.write_text(python_script)


    def pp_out(self):
        if self.inIO.software == 'STAR':
            bam_dct = {}
            for _sample_name, _sample_io in self.inIO.outdct.items():
                bam_dct[_sample_name] = _sample_io.bam.parent / (_sample_io.bam.name + '.toTranscriptome.out.bam')
        else:
            bam_dct = self.inIO.outdct
        if self.outIO.step_name == 'quantify':
            match self.outIO.software:
                case 'RSEM':
                    self._pp_rsem(bam_dct)
                case 'htseq-count':
                    self._pp_htseqcount(bam_dct)
        if self.outIO.step_name == 'MergeQuantify':
            match self.outIO.software:
                case 'RSEM':
                    self._pp_rsem_merge()
                case 'htseq-count':
                    self._pp_htseqcount_merge()
