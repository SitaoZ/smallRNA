#coding:utf-8

import gzip
import os,sys
import glob
from biosoftware import samtools,bowtie2,cmsearch
from optparse import OptionParser
from aasra_fq import aasra_fq

"""
Author: Sitao Zhu
Date : 2019-1-17
Desc : 注释各类小RNA

"""


sncrna_mapping_para = "-q -L 16 -p 6 --phred64 --sensitive"
genome_mapping_para = "-q -L 16 -p 6 --phred64 --rdg 1,10 --rfg 1,10"
cmsearch_para = "--cpu 6 --noali"
mismatch_for_bowtie2 = 1
unknown_tag_count_threshold = 10




def main():
    """
    %prog [options]
    """
    parser = OptionParser()
    parser.add_option('-d', '--database', help="database for the project")
    parser.add_option('-k','--sample',help="sample name")
    parser.add_option('-f','--fastq',help="fastq file")
    opts, args = parser.parse_args()
    if (opts.database == None or opts.sample == None or opts.fastq == None):
        print('\033[0;31;40m%s\033[0m' % "Warning: database, sample and fastq must be given\n")
        sys.exit(parser.print_help())
    database = opts.database
    sample = opts.sample
    fastq = opts.fastq
    shell_outdir = '.'
    annotation_outdir = '.'
    align(database,sample,fastq,shell_outdir,annotation_outdir)


def align(database,sample,fastq,shell_dir,annotation_dir):
    os.makedirs(shell_outdir)
    os.makedirs(annotation_outdir)
    for data in glob.glob(database+"/d_*"):
        name = os.path.basename(data).replace('d_','')
        SH = open(shell_dir+"/"+sample+"_"+name+".sh",'w')
        if name == 'sncrna':
            shell = 'python {dir}/aasra_fq.py {fq} {annot_dir}/{sample}.fq.gz'.format(dir=os.getcwd(),aasra_fq=aasra_fq,fq=fastq,annot_dir=annotation_dir,sample=sample)
            shell += "{bowtie2} {para} -x {data}/{name} -U {annot_dir}/{sample}.fq.gz | {samtools} view -bS -o {annot_dir}/{sample}_{name}.bam".format(bowtie2=bowtie2(),para=sncrna_mapping_para,data=data,name=name,annot_dir=annotation_dir,sample=sample,samtools=samtools())
            shell += 'rm -rf {annot_dir}/{sample}.fq.gz'.format(annot_dir=annotation_dir,sample=sample)
        elif name == 'Rfam':
            shell = 'python {dir}/fq2fa.py -i {annot_dir}/{sample}.fq.gz -o {annot_dir}/{sample}.fa'.format(dir=os.getcwd(),annot_dir=annotation_dir,sample=sample)
            shell += '{cmsearch} {para} {annot_dir}/{sample}.xls {data}/{name}.cm {annot_dir}/{sample}.fa'.format(cmsearch=cmsearch(),para=cmsearch_para,data=data,name=name)
        elif name == 'genome':
            # align filter read to genome
            shell = '{bowtie2} {para} -x {data}/{name} -U {annot_dir}/{sample}.fq.gz'
            shell += 'python {dir}/genome_stat.py {annot_dir}/{sample} {mismatch_for_bowtie2} {database}'.format(dir=os.getcwd(),annot_dir=annotation_dir,sample=sample,mismatch_for_bowtie2=mismatch_for_bowtie2,database=database)
        SH.close()
    with open(shell_dir+"/"+sample+"_catalog.sh",'w') as F:
        shell = 'python {dir}/catalog.py {database} {annot_dir}/{sample} {mismatch_for_bowtie2} {unknown_tag_count_threshold}'
        shell += "python {dir}/catalog_stat.py {annot_dir}/{sample}"


if __name__=="__main__":
    main()


