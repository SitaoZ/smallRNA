#-*-coding:utf-8-*-
import os,re
import sys
from optparse import OptionParser
from biosoftware import samtools


"""
http://chagall.med.cornell.edu/Mason_NGScourse/SAM.pdf
The MD field aims to achieve SNP/indel calling without looking at the reference. SOAP and Eland SNP callers prefer
such information. For example, a string “10A5^AC6” means from the leftmost reference base in the alignment, there
are 10 matches followed by an A on the reference which is different from the aligned read base; the next 5 reference
bases are matches followed by a 2bp deletion from the reference; the deleted sequence is AC; the last 6 bases are
matches. The MD field should match the CIGAR string, although an SAM parser may not check this optional field.

"""

def reverse_complement(seq:str):
    COMPLEMENT_TRANS = str.maketrans('TAGCtagc', 'ATCGATCG')
    return seq[::-1].translate(COMPLEMENT_TRANS)

def annot_sncrna(file_handle,database,prefix,misMatch,**ann):
    des = {} # miRNA description file
    dat = {} # given chrom count
    aln = {} # alignment information
    dele = {} # unmapped infor
    # process miRNA.des in database and store in des
    if os.path.exists(database+"/d_sncrna/miRNA.des"):
        with open(database+"/d_sncrna/miRNA.des",'r') as F:
            for line in F.readlines():
                line = line.strip()
                if line.startswith('#'):continue
                tmp = line.split('\t')
                des[tmp[0]] = tmp[1]
    # 处理sncrna比对bam文件 dat aln dele
    bam = os.popen("{samtools} view {prefix}_sncrna.bam".format(samtools=samtools(),prefix=prefix))
    for line in bam.readlines():
        line = line.strip()
        if line.startswith("@"):continue
        a = line.split()
        MD = ''
        if re.match('MD:Z:(\S+)',line):
            MD = re.match('MD:Z:(\S+)',line).group(1)  # MD:Z:0N0N52 - String for mismatching positions in the format of [0-9]+(([ACGTN]|\^[ACGTN]+)[0-9]+)*
        mismatch = MD.count('A') + MD.count('T') + MD.count('G') + MD.count('C') + MD.count('N')
        if tmp[1] == '16':
            tmp[9] = reverse_complement(tmp[9]) # revers and complemet
        tmp[9].replace('^CCCCC','').replace('GGGGG$','') # remove 3' CCCCC  5' GGGGG

        if aln.get(tmp[9]):
            dat[aln[tmp[9]]] += 1
            continue
        # put infor into  del,include unmapped, non mapped chromsome and mismatch count greater than misMatch
        if tmp[1] == 4 or tmp[2] == '*' or ann.get(tmp[9]) or mismatch > misMatch:
            dele[tmp[9]] += 1
            continue
        dat[tmp[2]] += 1  # save infor
        aln[tmp[9]] = tmp[2] # save infor
    # select infor in dele
    for de in dele.keys():
        if not aln.get(de): continue
        dat[aln[de]] += dele[de]
    # output
    pt = re.compile('(.*)\_([^\_]+)')
    for key in dat.keys():
        if pt.match(key):
            id ,db = pt.match(key).group(1),pt.match(key).group(2)
            if not des[id]:
                des[id] = id
            file_handle.writelines("\t".join([id,dat[key],db,des[id]])+"\n")
    for key in aln.keys():
        ann[key] = ''
    dat.clear()
    aln.clear()
    dele.clear()
    des.clear()

def blockread(fh, sep):
    """read a file by given separator
    :param fh:
    :param sep:
    :return:
    """
    buf = ""
    while True:
        while sep in buf:
            pos = buf.index(sep)
            yield buf[:pos]
            buf = buf[pos + len(sep):]
        chunk = fh.read(4096)
        if not chunk:
            yield buf
            break
        buf += chunk

def annot_rfam(file_handle,database,prefix,**ann):
    seq = {}
    des = {}
    tag = {}
    dat = {}
    with open(prefix+".fa",'r') as F:
        for line in F.readlines():
            line = line.strip()
            if line.startswith('>'):
                hd = line.replace('>','')
            else:
                seq[hd] = line
    with open(database+"/d_Rfam/Rfam.des",'r') as F:
        for line in F.readlines():
            if len(line) == 0 :continue
            line = line.strip()
            tmp = line.split()
            typ = 'Rfam other sncRNA'
            if 'miRNA' in tmp[3]:
                typ = 'precursor'
            elif 'snoRNA' in tmp[3]:
                typ = 'snoRNA'
            elif 'snRNA' in tmp[3]:
                typ = 'snRNA'
            elif 'rRNA' in tmp[3]:
                typ = 'rRNA'
            des[tmp[1]] = '\t'.join([typ,tmp[4]])
    Rfam_xls = open(prefix+"_Rfam.xls",'r')
    for line in blockread(Rfam_xls,'//'):
        tmp = line.split('\n')
        query,acc = '',''
        if re.match('Query:\s+(\S+).*',line):
            query = re.match('Query:\s+(\S+).*',line).group(1)
        if re.match('Accession:\s+(\S+).*',line):
            acc = re.match('Accession:\s+(\S+).*',line).group(1)
        if re.match('^\s+\(\d+\).*',line):
            tmp = line.split('\s')
            if not ann.get(tmp[6]):
                if not tag.get(tmp[6]):
                    tag[tmp[6]] = "\t".join([acc,tmp[4]])
                else:
                    aacc = tag[tmp[6]].split('\t')
                    if float(tmp[4]) >= float(aacc[1]):
                        tag[tmp[6]] = "\t".join([acc,tmp[4]])
    for key in tag.keys():
        atag = key.split('_')
        aacc = tag[key].split('\t')
        dat[aacc[0]] += atag[1]
        ann[seq[key]] = ''
    for key in dat.keys():
        file_handle.writelines("\t".join([key,dat[key],des[key]])+"\n")
    dat.clear()
    tag.clear()
    seq.clear()
    des.clear()


def region_load(repeat_region,reg,flag,RT):
    with open(repeat_region,'r') as F:
        for line in F.readlines():
            line = line.strip()
            if flag == 'repeat':
                # repeat
                if re.match(r'^\s+?\d+',line):
                    tmp = line.split()
                    strand = "+"
                    if tmp[9] == 'C':
                        strand = '-'
                    RT.writelines("\t".join([tmp[5],tmp[6],tmp[7],tmp[9],'repeat']))
            else:
                # region
                RT.writelines(line+"\n")

def catalog(database,prefix,mismatch=1,novel=10):
    OA = open(prefix+".catalog.xls",'w')
    OA.writelines("sRNA id\tCount\tType\tDescription\n")
    ann ={}
    annot_sncrna(OA,database,prefix,mismatch,**ann)  # annotation sncrna, sncrna classfication
    annot_rfam(OA,database,prefix,**ann)    # annotation rfam database, rfam classfication
    reg,dat,bit,aln = {},{},{},{}
    aln = {} # save all filter read alignment genome information
    # repeat and region
    RT = open(prefix+".region",'w')
    repeat = opts.database+"/genome_annot/genome.repeat"
    region = opts.database+"/genome_annot/genome.region"
    region_load(repeat,reg,'repeat',RT)
    region_load(region,reg,'na',RT)
    RT.close()
    #global variant
    ann = {}
    unmap = 0
    OU = open(prefix+".unknown.fa",'w') # tag1_61301_Pn7 19068351        -
    # the unknown defination: if the read mapping to genome  more than 10 positions,will be unknown fa
    bam = os.popen("{samtools} view {prefix}_genome.bam".format(samtools=samtools(), prefix=prefix))
    for line in bam:
        line = line.strip()
        if line.startswith('^@'):continue
        tmp = line.split('\t')   #bam column
        strand = '+'
        if tmp[1] == '16':
            reverse_complement(tmp[9])
            strand = '-'
        if ann.get(tmp[9]):continue
        MD = ''
        if re.match('MD:Z:(\S+)', line):
            MD = re.match('MD:Z:(\S+)', line).group(1)
        mismatch_c = MD.count('A') + MD.count('T') + MD.count('G') + MD.count('C') + MD.count('N')
        if tmp[1] == 4 or tmp[2] == '*' or mismatch_c > mismatch:
            unmap += 1
            continue
        end = int(tmp[3]) + len(tmp[9]) - 1
        if not aln.get(tmp[9]):
            aln[tmp[9]] = ["\t".join([tmp[2],tmp[3],end,strand]),1] # chromosome,position,end,strand,count
        aln[tmp[9]][1] += 1  #count accumulation
    TM = open(prefix+'.map.xls','w')
    flag = 1
    for seq in sorted(aln.keys(),key=lambda ele:ele[1]):
        inf = aln[seq][0].split('\t')
        if aln[seq][1] >= novel:
            #>tagid_count_chromosome position strand    example: tag1_61301_Pn7 19068351        -
            # seq
            OU.writelines(">tag{a}_{b}_{c}\t{d}\t{e}\n{f}\n".format(a=flag,b=aln[seq][1],c=inf[0],d=inf[1],e=inf[3],f=seq))
        TM.writelines("{f}_{b}\t{c}\n".format(f=seq,b=aln[seq][1],c=aln[seq][0]))
        flag += 1
    TM.close()
    os.system("TGR {prefix}.map.xls {prefix}.region {prefix}.aln".format(prefix=prefix)) # according map ,repeat and region to generate aln
    saln = {}
    IA=open(prefix+".aln",'r')
    for line in IA.readlines():
        tmp = line.strip().split()
        tag = tmp[0].split('_')
        if aln.get(tag[0]):
            if saln.get(tmp[2]):
               saln[tmp[2]] += tag[1] # tmp[2]:mirRNA id   tag[1]:count
            else:
                saln[tmp[2]] = tag[1]
        del aln[tag[0]] # delete read not in region
    for key in saln.keys():
        OA.writelines("{a}\t{b}\t{c}\t-\n".format(a=key,b=saln[key],c=key)) ## exon intron repeat classfication
    inter = 0
    for min in aln.keys():
        inter += float(aln[min][1]) # print intergenic
    OA.writelines("intergenic\t{inter}\tintergenic\t-\n".format(inter=inter)) # intergenic classfication
    OA.writelines("unmap\t{unmap}\tunmap\t-\n".format(unmap=unmap))           # unmap classfication
    OU.close()
    OA.close()
    os.system('rm -rf {prefix}.map.xls {prefix}.region {prefix}.aln'.format(prefix=prefix))


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
    if (opts.database == None or opts.prefix == None or opts.mismatch == None,opts.novel_threshold == None):
        print('\033[0;31;40m%s\033[0m' % "Warning: database, prefix, mismatch and novel_threshold must be given\n")
        sys.exit(parser.print_help())
    misMatch = opts.mismatch
    dataBase = opts.database
    prefix = opts.prefix
    novel_threshold = opts.novel_threshold
    catalog(dataBase,prefix,misMatch,novel_threshold)


if __name__=="__main__":
    main()









