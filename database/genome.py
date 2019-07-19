#coding:utf-8

"""
Author: Sitao Zhu
Date  : 2019-1-1
Desc  : 基因组建立索引；同时完成重复序列注释；python3.6
"""

import os
import sys

from biosoftware import bowtie2
from biosoftware import bowtie1
from biosoftware import RepeatMasker

def build_index(genome,prefix):
    os.system("{bt1} -q -f {genome} {prefix}".format(bt1=bt1,genome=genome,prefix=prefix))
    os.system("{bt2} -q -f {genome} {prefix}".format(bt2=bt2,genome=genome,prefix=prefix))

def repeat(genome,outdir):
    os.system("{repeatmakser} {genome} -dir {outdir}".format(repeatmakser=RepeatMasker(),genome=genome,outdir=outdir))
    repeatout = outdir
    os.system("mv {outdir}/*.out {outdir}/genome.repeat".format(outdir=outdir))


