#coding:utf-8
"""
Author:Sitao Zhu
Date : 20190125
Desc: the KO.py was used to blast the gene fasta with kegg database
"""
import os,re,sys
import glob
import datetime
from subprocess import Popen,PIPE
from multiprocessing import Process
sys.path.append('')
from fasta import Fasta



def split_ko(kodatabase,outdir,tag='ko',split_number=30):
    """ 切分ko数据库"""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fasta_obj = Fasta(kodatabase)
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
        OT.writelines(fasta_obj.printf(key,fasta_obj._fasta[key]))
    OT.close()
    ko_fasta_files = glob.glob(r'%s/*.fa' % tmp_outdir)
    return ko_fasta_files


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

def parse_ko(kegg):
    """解析kegg数据库"""
    kos = {}
    pattern = re.compile(r'^>(?P<name>\S+)\s.*(?P<ko>K\d+)\s(?P<function>[^\;]+).*')
    with open(kegg, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if ">" not in line:
                continue
            m = pattern.search(line)
            if m:
                #print (m.group('function'))
                kos[m.group('name')] = m.group('ko')+"|"+m.group('function')
    return kos

def bsort(ele):
    """根据特定identity值进行排序"""
    tmp = ele.split('|')
    return float(ele[4])

def parse_m8(m8,kos:dict):
    """解析blast m8的结果"""
    from collections import defaultdict
    dat = {}
    with open(m8,'r') as f:
        for line in f.readlines():
            tmp = line.strip().split()
            if tmp[1] not in kos.keys():
                continue
            tmp[-1] = re.sub(r'\s','',tmp[-1])
            tmp[0] = re.sub(r'\s+','',tmp[0])
            if tmp[0] not in dat.keys():
                dat[tmp[0]] = '|'.join([kos[tmp[1]],tmp[1],tmp[-2],tmp[-1]+"_kos_"])
            else:
                #存储多个比对结果
                dat[tmp[0]] += '|'.join([kos[tmp[1]],tmp[1],tmp[-2],tmp[-1]+"_kos_"])
    ## sort 根据比对的blast identity的值进行排序
    for key in dat.keys():
        dat[key] = re.sub(r'_kos_$','',dat[key])
        flag = [1,'']
        tmp = dat[key].split('_kos_')
        ##对多个比对结果进行排序，选取前五
        for i in sorted(tmp,key=lambda ele_tmp: float(ele_tmp.split('|')[4]),reverse=True):
            inf = i.split('|')
            flag[1] += inf[0]+"|"+str(flag[0])+"|"+inf[3]+"|"+inf[4]+"|"+inf[2]+"|"+inf[1]+'!'
            flag[0] += 1
            if flag[0] > 5:
                break
        dat[key] = re.sub(r'\!$','',flag[1])
    return dat

def blast(blast_path,kegg,query,prefix,outdir,evalue='1e-5'):
    """blast 主程序，解析blast kegg 的结果，进行提取，得到基因对应kegg编号的对应关系"""
    pattern = re.compile(r'>(?P<id>\S+).*')
    #step1 formatdb
    formatdb = os.path.join(blast_path,'formatdb')
    cmd1 = [formatdb,'-p T -i', kegg]
    runSysCommandSubProcess(cmd1,showInLog=True,justTest=False)
    #step 2 blastall
    blastall = os.path.join(blast_path,'blastall')
    m8result = os.path.join(outdir,prefix + ".m8")
    cmd2 = [blastall, '-d', kegg, '-i', query, '-p', 'blastx -m 8 -a 12 -e', evalue, '-o', m8result]
    runSysCommandSubProcess(cmd2,showInLog=True,justTest=False)
    #step 3 deal with result
    kos = parse_ko(kegg)
    dat = parse_m8(m8result,kos)
    #print(dat)
    ko_result = os.path.join(outdir,prefix+".ko")
    OT=open(ko_result,'w')
    OT.writelines("#query\tko:rank:evalue:score:gene:definition\n")
    with open(query,'r') as f:
        for line in f.readlines():
            line = line.strip()
            if pattern.match(line):
                ID = pattern.match(line).group('id')
                if ID not in dat.keys():
                    dat[ID] = ''
                outline = ID+"\t"+dat[ID]
                OT.writelines(outline+"\n")
    OT.close()

def annotation(blast_path,kegg,query,outdir,tag='ko'):
    """解析全部的ko.fa的比对结果，进行汇总,上面定义的函数都是处理单个ko.fa的"""
    ko_fasta_files = split_ko(kegg,outdir,tag='ko',split_number=30)
    print (ko_fasta_files)
    #p = Process(target=blast(blast_path,query,),)
    for fa in ko_fasta_files:
        prefix = fa.replace('.fa','')
        blast(blast_path,fa,query,prefix,outdir,evalue='1e-5')





if __name__=="__main__":
    #split_ko("/ifs4/BC_PUB/biosoft/database/PUB/Pub/kegg/RNA/79/plant.fa","/nascngb/gccnt2/ST_BI/zhusitao/05.test/small/00.data/ko3")
    #blast("/share/app/blast-2.2.26/bin","ko8.fa",'Piper_nigrum.annotation.cds','ko8',".")
    annotation('/share/app/blast-2.2.26/bin','/ifs4/BC_PUB/biosoft/database/PUB/Pub/kegg/RNA/79/plant.fa','/hwfssz5/ST_BIGDATA/USER/zhusitao/Project/micro/00.data/v2/Piper_nigrum.annotation.cds','/nascngb/gccnt2/ST_BI/zhusitao/05.test/small/00.data/ko',tag='ko')
