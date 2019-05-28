#coding:utf-8

"""
Author： Sitao Zhu
Date  :  2019-1-1
Desc  :  提取sncRNA 中的非编码RNA,python版本是3.6
            采用牟定比对的方法 https://github.com/biogramming/AASRA
"""
import os

from biosoftware import bowtie2

def Aasra(species,dbs,outdir):
    """ extract sncRNA 提取所有的小的非编码RNA，包括snoRNA（小核仁RNA），mature， hairpin """

    OF=open(outdir+"/sncrna.fa",'w')
    idlist = []
    seqlist = []
    db_class = ["snoRNA", "mature", "hairpin"]
    for i in range(len(dbs)):
        if not os.path.exists(dbs[i]):
            continue
        with open(dbs[i], 'r') as F:
            for line in F.readlines():
                line = line.strip()
                if line.startswith(">"):
                    rid = line
                    name = rid + "_" + db_class[i]
                    idlist.append(name)
                else:
                    sq = line
                    # https://github.com/biogramming/AASRA/blob/master/README.md
                    seq = "C" * 10 + sq + "G" * 10
                    seqlist.append(seq)

    for i in range(len(idlist)):
        lineout = "\n".join([idlist[i],seqlist[i]])
        OF.writelines(lineout+"\n")

    OF.close()




