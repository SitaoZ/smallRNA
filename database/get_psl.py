#coding:utf-8

"""
Author： Sitao Zhu
Date  :  2019-1-1
Desc  :  利用blat 将转录本序列比对回基因组得倒基因组区域信息,python版本是3.6
         psl file formate https://asia.ensembl.org/info/website/upload/psl.html
"""

import os,re


def psl(blat,genome,transcript,pslfile):
    """ genome : genome fasta file; transcript fasta file"""
    cmd = [blat,'-noHead',genome,transcript,pslfile]
    return cmd

def psl_filter(pslfile,outdir):
    OT = open(outdir+"/tmp.psl",'w')
    with open(pslfile,'r') as F:
        gene = 'na'
        psl = {}
        for line in F.readlines():
            if len(line) == 0:
                continue
            line = line.strip()
            tmp = line.split()
            if gene != tmp[9] and gene != 'na':
                array = sorted(psl.keys())
                OT.writelines(array[0]+"\n")
                psl = {}
            gene = tmp[9]
            psl[line] = int(tmp[0]) - int(tmp[1])
    array = sorted(psl.keys())
    OT.writelines(array[0]+"\n")
    OT.close()
    return outdir+"/tmp.psl"

def genome_region(tmp_psl,outdir):

    OT = open(outdir+"/genome.region",'w')
    with open(tmp_psl,'r') as F:
        for line in F.readlines():
            line = line.strip()
            tmp = line.split()
            tmp[18] = re.sub(',$','',tmp[18])
            tmp[20] = re.sub(',$','',tmp[20])
            length = list(map(int,tmp[18].split(',')))
            start = list(map(int,tmp[20].split(',')))
            for i in range(len(start)):
                lineout = "\t".join([tmp[13],str(start[i]+1),str(start[i]+length[i]),tmp[8],"exon"])
                OT.writelines(lineout+"\n")
                if i < len(start)-1:
                    lineout = "\t".join([tmp[13],str(start[i]+length[i]+1),str(start[i+1]),tmp[8],"intron"])
                    OT.writelines(lineout+"\n")
    return outdir+"/genome.region"



