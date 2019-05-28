#coding:utf-8

import gzip
import os,sys
import glob
from biosoftware import samtools,bowtie2
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
        shell = ''
        if name == 'sncrna':
            shell = 'python {dir}/{aasra_fq} {fq} {annot_dir}/{sample}.fq.gz'.format(dir=os.getcwd(),aasra_fq=aasra_fq,fq=fastq,annot_dir=annotation_dir,sample=sample)
            shell += "{bowtie2} {para} -x {data}/{name} -U {annot_dir}/{sample}.fq.gz | {samtools} view -bS -o {annot_dir}/{sample}_{name}.bam".format(bowtie2=bowtie2(),para=sncrna_mapping_para,data=data,name=name,annot_dir=annotation_dir,sample=sample,samtools=samtools())
            shell += 'rm -rf {annot_dir}/{sample}.fq.gz'.format(annot_dir=annotation_dir,sample=sample)
        elif name == 'Rfam':
            shell = 'python {dir}/fq2fa.py -i {annot_dir}/{sample}.fq.gz -o {annot_dir}/{sample}.fa'.format(dir=os.getcwd(),annot_dir=annotation_dir,sample=sample)

def Genome_mapping(self,out_bam):
    """ bowtie2 mapping """
    genome_prefix = self.database_path+"/d_genome/genome"
    bowtie2 = bowtie2()
    samtools = samtools()
    self._log('Start smallRNA annotation: genome_mapping task', showTime=True)
    cmd = [bowtie2,self.genome_mapping_para,genome_prefix,self.fastq,samtools,out_bam]
    self.runSysCommandSubProcess(cmd, showInLog=True, justTest=True)
    self._log('Finish smallRNA annotation: genome_mapping task', showTime=True)

def Genome_unknown(self,unknown_fa):
    """获取 基因组 不知道的fasta文件"""
    pass







def Rfam_cmsearch(self,in_fa,outxls):
    """ 注释 Rfam"""
    cmsearch = cmsearch()
    Rfam = database.Rfam()
    self._log('Start smallRNA annotation: Rfam task', showTime=True)
    cmd = [cmsearch,'--cpu 6 --noali -o',outxls,Rfam,in_fa]
    self.runSysCommandSubProcess(cmd, showInLog=True, justTest=True)
    self._log('Finish smallRNA annotation: Rfam task', showTime=True)


def Rfam_annotation(self):
    """将cmsearch的结果与Rfam.des数据库中的比较，得到注释信息"""
    pass



def sncrna_mapping(self,out_bam):
    """ 将过滤后的reads比对回去sncrna.fa ，之后用于注释"""
    sncrna_prefix = self.database_path+"/d_sncrna/sncrna"
    bowtie2 = bowtie2()
    samtools = samtools()
    self._log('Start smallRNA annotation: sncrna task', showTime=True)
    cmd = [bowtie2, self.sncrna_mapping_para, '-x', sncrna_prefix, '-U', self.fastq, '|', samtools, "view -bS -o", out_bam]
    self.runSysCommandSubProcess(cmd, showInLog=True, justTest=True)
    self._log('Finish smallRNA annotation: sncrna task', showTime=True)
def sncrna_anotation(self):
    pass

if __name__=="__main__":
    main()


