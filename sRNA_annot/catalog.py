#-*-coding:utf-8-*-
import os,re
import sys
from optparse import OptionParser
from biosoftware import samtools


def main():
    """
    %prog [options]
    """
    parser = OptionParser()
    parser.add_option('-d', '--database', help="database for the project")
    parser.add_option('-p','--prefix',help="the prefix")
    parser.add_option('-m','--mismatch',help="mismatch")
    parser.add_option('-n', '--novel_threshold', help="novel count threshold")
    opts, args = parser.parse_args()
    if (opts.database == None or opts.prefix == None or opts.mismatch == None):
        print('\033[0;31;40m%s\033[0m' % "Warning: database, prefix, mismatch and novel_threshold must be given\n")
        sys.exit(parser.print_help())




def reverse_complement(seq:str):
    COMPLEMENT_TRANS = str.maketrans('TAGCtagc', 'ATCGATCG')
    return seq[::-1].translate(COMPLEMENT_TRANS)

def annot_ancrna(file_handle,database,prefix):
    des = {} # miRNA description file
    aln = {} # aligner bam dict
    dat = {}
    # 处理数据库中的miRNA.des文件，存入字典
    if os.path.exists(database+"/d_sncrna/miRNA.des"):
        with open(database+"/d_sncrna/miRNA.des",'r') as F:
            for line in F.readlines():
                line = line.strip()
                if line.startswith('#'):continue
                tmp = line.split('\t')
                des[tmp[0]] = tmp[1]
    # 处理sncrna比对bam文件
    bam = os.popen("{samtools} view {prefix}_sncrna.bam".format(samtools=samtools(),prefix=prefix))
    for line in bam.readlines():
        line = line.strip()
        if line.startswith("@"):continue
        a = line.split()
        MD= 0
        if re.match('MD:Z:(\S+)',line):
            MD = re.match('MD:Z:(\S+)',line).group(1)
        if tmp[1] == '16':
            tmp[9] = reverse_complement(tmp[9]) # 反向互补负链
        tmp[9].replace('^CCCCC','').replace('GGGGG$','') # 除去加上的碱基
        if tmp[9] in aln.keys():
            dat[aln[tmp[9]]]


        dat[tmp[2]] += 1

        aln[tmp[9]] = tmp[2]





def catalog(database,prefix,mimatch,novel):
    OA = open(prefix+".catalog.xls",'w')
    OA.writelines("sRNA id\tCount\tType\tDescription\n")
    annot_ancrna(OA,database,prefix)





def annot_rfam():
    pass



