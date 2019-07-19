#conding:utf-8

"""
Author： Sitao Zhu
Date  :  2019-1-1
Desc  :  提取数据库micBase中的小RNA ,python版本是3.6
"""

import os
import sys
import re

def known(spe,micbase,outdir):
    """ 从官方的mirBase数据库提取micBase已知的 mature.fa miRNA.fa """
    speciesPattern = re.compile(r'(\w+)-')
    mature_fa = outdir+"/mature.fa"
    miRNA_fa = outdir + "/miRNA.fa"
    miRNA_des = outdir + "/miRNA.des"

    mature = open(mature_fa,'w')
    miRNA = open(miRNA_fa,'w')
    miRNAdes = open(miRNA_des,'w')

    header = "\t".join(["#miRNA id","family","mature","hairpin","end pos(mature on harpin)"])
    miRNAdes.writelines(header+"\n")
    dele = {}
    with open(micbase,'r') as F:
        for line in F.readlines():
            line = line.strip()
            tmp = line.split()
            if speciesPattern.match(tmp[0]):
                species = speciesPattern.match(tmp[0]).group(1)
            end = tmp[4].rindex(tmp[3])
            if tmp[0] in dele or species != spe:
                continue
            #miRNA 由初级转录本经过一系列核酸酶剪切而成，长度几百到几千个碱基不等，有一个至数个发卡结构
            line1 = "\n".join([">"+tmp[0],tmp[4]])
            miRNA.writelines(line1+"\n")
            #miRNA describe
            length = str(end + len(tmp[3]))
            line2 = "\t".join([tmp[0],tmp[2],tmp[3],tmp[4],length])
            miRNAdes.writelines(line2+"\n")
            #mature 成熟的发夹结构
            line3 = "\n".join([">"+tmp[0],tmp[3]])
            mature.writelines(line3+"\n")
    mature.close()
    miRNA.close()
    miRNAdes.close()
    return mature_fa,miRNA_fa,miRNA_des


def species():
    pass

            
