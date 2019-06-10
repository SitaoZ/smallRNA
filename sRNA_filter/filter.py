#-*-coding:utf-8-*-
import re,os
from optparse import OptionParser

def filter(fastq,sample_name,clean,adaptor,quality,shell_dir,read_dir):
    path = os.path.dirname(fastq)
    shell = ''
    adt = []
    adaptor_list = ['GTTCAGAGTTCTACAGTCCGACGATC', 'TCGTATGCCGTCTTCTGCTTGT', 'ATCTCGTATGCCGTCTTCTGCTTG', 'TGGAATTCTCGGGTGCCAAGG', 'AGATCGGAAGAGCACACGTCT']
    if os.path.exists(path,'1.fqcheck'):
        shell += 'tile=`perl $Bin/findNtile.pl -fq1 $opts{f} -fc1 $path/1.fqcheck` \n'
    else:
        shell += 'tile=9999'
    if re.match(r'(\w+)\,(\w+)',adaptor):
        adt = re.match(r'(\w+)\,(\w+)',adaptor).group(0)
    elif adaptor == 5:
        adt = ['GAACGACATGGCTACGATCCGACTT', 'AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA']
    elif adaptor_list[adaptor]:
        adt = adaptor_list[0],adaptor_list[adaptor]
        
