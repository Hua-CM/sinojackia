# sinojackia
**sinojackia** is a semi-automated pipeline commands generator for RNA/WGS analysis.

## Bigupdae:V2.0

The last update time: 2024/01/01

- RNA-seq now support htseq-count.
- WGS now support bcftools, bowtie, and bowtie2.
- Replace all `os` package code with `pathlib`
- Totally new logical
- Fix some bugs

## Introduction
**sinojackia** was influenced by [**Ilus**](https://github.com/ShujiaHuang/ilus) developed by Shujia Huang. It only generates commands for analysis, instead runing these commands directly. For now, It comprises two modules: 
1. WGS
2. RNA-seq

## Requirements

## Meta file
There are two types of meta file for different starting points. **If you have reads from different lanes, please merge them manually before analysis**

### Type1
**Applicable starting point:** QC, align, 

| sample | RG | fastq1 | fastq2 |
| :---: | :---: | :---: | :---: |
| G01 | "@RG\tID:G01\tSM:G01\tPL:ILLUMINA\tLB:G01" | /path/to/fastq/G01_L4_1_clean.fq | /path/to/fastq/G01_L4_2_clean.fq |

### Type2
**Applicable starting point:** call, quantity

| sample | bam |
| :---: | :---: | 
| G01 | /path/to/bam/G01.bam |

### Type3
If you only run the `genotype`. Just give your `combine.g.vcf.gz` path in command line. 

## Reference
[Ilus](https://github.com/ShujiaHuang/ilus)

## TODO

1. QC section will add FastQC
2. GATK accept will be accept parameters