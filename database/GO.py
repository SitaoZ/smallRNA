#coding:utf-8
"""
Author:Sitao Zhu
Date : 20190125
Desc: the GO.py was used to blast the gene fasta with kegg database
"""
import os,re,sys
import glob
import datetime
from fasta import Fasta
from subprocess import Popen,PIPE

def split_go(godatabase,outdir,tag='go',split_number=30):
    """ 切分go数据库"""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fasta_obj = Fasta(godatabase)
    icount = fasta_obj.id_num()
    piece = icount // split_number + 1
    num, flag = 0, 1
    tmp_outdir = outdir + "/tmp_seq_" + tag
    if not os.path.exists(tmp_outdir):
        os.makedirs(tmp_outdir)
    tmp_fa = os.path.join(tmp_outdir,"%s%d.fa"%(tag,flag))
    OT = open(tmp_fa,'w')
    for key in fasta_obj._fasta.keys():
        if num == piece:
            flag += 1
            num = 0
            tmp_fa = os.path.join(tmp_outdir, "%s%d.fa" % (tag, flag))
            OT = open(tmp_fa,'w')
        num +=1
        OT.writelines(fasta_obj.printf(key,fasta_obj._fasta[key],number=80))
    OT.close()
    go_fasta_files = glob.glob(r'%s/*.fa' % tmp_outdir)
    return go_fasta_files

def run_log(content, showTime=True):
    ##### Display text to defined stream, add timestamp
    logStream = sys.stdout
    if isinstance(content, list):
        modifText = " ".join(content)
    else:
        modifText = str(content)

    if showTime == True:
        modifText += "    " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    else:
        return modifText
def runSysCommandSubProcess(command, showInLog=True, justTest=False):
    if showInLog == True:
        run_log(command, showTime=False)
        if justTest == True:
            return
    p = Popen(" ".join(command), shell=True, stderr=PIPE)
    return p.stderr.read().decode()

def blast(blast_path,go,query,prefix,outdir,evalue='1e-5'):
    """进行blast比对"""
    # step formatdb
    formatdb = os.path.join(blast_path, 'formatdb')
    cmd1 = [formatdb, '-p T -i', go]
    runSysCommandSubProcess(cmd1, showInLog=True, justTest=False)
    # step 2 blastall
    blastall = os.path.join(blast_path, 'blastall')
    xml_result = os.path.join(outdir,prefix + ".xml")
    cmd2 = [blastall, '-d', go, '-i', query, '-p', 'blastx -F F -m 7 -a 12 -v 2 -b 2 -e', evalue, '-o', xml_result]
    runSysCommandSubProcess(cmd2, showInLog=True, justTest=False)
    return xml_result


def read_blocks(file):
    """ Python Reads the file with the specified delimiter """
    block = ''
    for line in file:
        if line.startswith('<Iteration>') and len(block) > 0:
            yield block
            block = ''
        block += line
    yield block

def blast2go(prefix):
    """解析blast 的比对结果,解析 xml文件"""
    Hit_id = re.compile(r'<Hit_id>.*<\/Hit_id>')
    Hit_def = re.compile(r'<Hit_def>(\S+)\s.*<\/Hit_def>')

    with open(xml_result,'r') as f:
        flag = 0
        txml = open(prefix+".txml",'w')
        for line in f.readlines():
            line = line.strip()
            if line.startswith('<Iteration>'):
                flag += 1
                if flag > 2000:
                    # 将大的blast比对文件拆分成小的，执行blast2go java 程序
                    # 执行程序
                    txml.writelines("</BlastOutput_iterations>\n</BlastOutput>\n")
                    txml = open(prefix + ".txml", 'w')
                    flag = 0
                else:
                    txml.writelines(line)
            elif Hit_id.match(line):
                continue
            elif Hit_def.match(line):
                key = Hit_def.match(line).group(1)
                txml.writelines("<Hit_id>{0}<\/Hit_id>".format(key))
                txml.writelines(line)
            else:
                txml.writelines(line)




def nr_desc():
    pass


if __name__=="__main__":
    #split_go('/ifs4/BC_PUB/biosoft/database/PUB/Pub/nr/RNA/20160219/Plants.fa',"/nascngb/gccnt2/ST_BI/zhusitao/05.test/small/00.data/go")
    blast("/share/app/blast-2.2.26/bin","/nascngb/gccnt2/ST_BI/zhusitao/05.test/small/00.data/go/tmp_seq_go/go1.fa",'Piper_nigrum.annotation.cds','go1',".")