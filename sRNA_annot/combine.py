#-*-coding:utf-8-*-
#!/usr/bin/python

import os,re,sys
import subprocess
from optparse import OptionParser
from util.fasta import Fasta
from collections import defaultdict




def main():
	"""
	%prog [options]
	:return:
	"""
	parser = OptionParser()
	parser.add_option('-k','--purpose',type='int',help="purpose, 0:unknown 1:express 2:express_aasra")
	parser.add_option('-f','--sample',type='string',help="prefix of each sample annotation,separated bu comma")
	parser.add_option('-p','--predict',type='string',help="prediction outdir")
	parser.add_option('-d','--dir',type='string',help="expression outdir for degs")

	opts,args = parser.parse_args()
	if opts.purpose == None or opts.sample == None or opts.predict ==None :
		print('\033[0;31;40m%s\033[0m' % "Warning: purpose, samples and predict_dir must be given\n")
		sys.exit(parser.print_help())
	purpose = opts.purpose
	samples = opts.sample
	predict_dir = opts.predict
	express_dir = opts.dir
	if purpose == 0:
		combine_unkonwn(predict_dir,samples)
	if purpose == 1:
		print(type(purpose))
		combine_express(predict_dir,samples,express_dir)



def combine_unkonwn(predict_dir,samples):
	"""
	combine all sample unknown read which map multiply(10 locites) to genome according to seequnce
	unknown.fa source from catalog.py(mapping to genome.fa with 10bp and not annotated by ncrna database)
	>tagid_count_chromosome position strand    example: tag1_61301_Pn7 19068351        -
	:param predict_dir: predict outdir
	:param samples: samples prefix partition by comma
	:return: none
	"""
	if not os.path.exists(predict_dir):
		os.system('mkdir %s'%(predict_dir))
	dat = defaultdict(dict)
	key = []
	## in
	for prefix in samples.split(','):
		sample = os.path.basename(prefix)
		key.append(sample)
		fasta_obj = Fasta(prefix+'.unknown.fa')
		for record in fasta_obj:
			name = record.id
			seq = record.seq
			hd = name.split('_') # name formt #>tagid_count_chromosome position strand . count is read mapping position count
			chr_start_strand = hd[2].split('\t')
			end = int(chr_start_strand[1]) + len(seq) - 1
			chr_start_strand.insert(2,str(end)) # add end
			dat[seq]['map'] = "\t".join(chr_start_strand)
			dat[seq][sample] = int(hd[1]) # count

	## out
	OF = open(predict_dir+"/unknown.fa",'w')  # Take the intersection of multiple samples
	OD = open(predict_dir+"/unknown.xls",'w') # print alignment information for each unknown tag
	OD.writelines("Tag\tChr\tStart\tEnd\tStrand\t"+"\t".join(key)+"\tSequence\n") # key 顺序固定
	flag = 1
	for seq in dat.keys():
		count = 0
		OD.writelines("tag{flag}\t{map}".format(flag=flag,map=dat[seq]['map'])) #tag chr start end strand
		for sample in key:
			if not dat[seq].get(sample):
				dat[seq][sample] = 0
			count += dat[seq][sample] # 每一个序列特定下，把所有的个体的该序列的条数累加
			OD.writelines("\t"+str(dat[seq][sample])) # write count of every sample
		OD.writelines("\t"+seq+"\n") # sequence
		OF.writelines(">tag{flag}_{count}\n{seq}\n".format(flag=flag,count=count,seq=seq)) # fasta
		flag += 1
	OF.close()
	OD.close()

def predicted_sRNA(file_handle,exp,tag,key,typ):
	"""
	Add predict sRNA information to sRNA which had been annotated by database
	:param file_handle: predict miRNA or piRNA or siRNA file
	:param exp: a dict data from combine_express
	:param tag: a dixt data from combine_express
	:param key: a list data from combine_express
	:param typ: a string from combine_express
	:return: modified exp tag key typ
	"""
	with open(file_handle,'r') as F:
		for line in F.readlines():
			line = line.strip()
			if re.search('Tag info',line):continue
			a = line.split('\t')
			if typ == 'miRNA':
				des = a[1]+":"+":"+a[2]+":"+a[5]+"_"+a[6]
			if typ == 'piRNA':
				des = a[1]+"_"+a[2]+":"+a[3]+":"+a[4]+"_"+a[5]
			ttag = []
			for t in a[-1].split(','):
				t = re.sub('\_.*','',t)
				if t not in tag:continue
				ttag.append(t)
			if len(ttag) == 0 :continue
			i = 0
			while(i < len(key)):
				c = 0
				for b in ttag:
					c += tag[b][i]
				if c == 0 : continue
				exp[key[i][1]].append("{0}\t{1}\t{2}\t{3}".format(a[0],str(c),typ,des))
				exp[key[i][2]] += c
			for i in ttag:
				del tag[i]




def combine_express(predict_dir,samples,express_dir):
	"""
	combine known and unknown expression file
	:param predict_dir: predict unknown out dir
	:param samples: samples prefix partition comma,string
	:param express_dir: expression output dir,string
	:return:none
	"""
	tag = {}
	keys = []
	exp = defaultdict(dict)
	if not os.path.exists(express_dir):os.mkdir(express_dir)
	if os.path.exists(predict_dir+"/unknown.xls"):
		with open(predict_dir+"/unknown.xls",'r') as F:
			for line in F.readlines():
				a = line.strip().split()
				if a[0] == 'Tag':
					keys = a[5:]
				else:
					tag[a[0]] = a[5:]

	if os.path.exists(predict_dir+"/miRNA/predicted_miRNA.xls"):
		predicted_sRNA(predict_dir+"/miRNA/predicted_miRNA.xls",exp,tag,keys,'miRNA')
	if os.path.exists(predict_dir + "/miRNA/predicted_piRNA.xls"):
		predicted_sRNA(predict_dir + "/miRNA/predicted_piRNA.xls", exp, tag, keys, 'piRNA')
	if os.path.exists(predict_dir + "/miRNA/predicted_siRNA.xls"):
		predicted_sRNA(predict_dir + "/miRNA/predicted_siRNA.xls", exp, tag, keys, 'siRNA')
	flag= defaultdict(int)
	for prefix in samples.split(','):
		base = os.path.basename(prefix)
		exp[base][1] = []
		exp[base][2] = 0
		with open(prefix+".catalog.xls",'r') as F:
			for line in F.readlines():
				line = line.strip()
				a = line.split('\t')
				if a[0] == 'sRNA id' or a[0] == 'unmap':continue
				exp[base][1].append(line) #save all samples expression
				exp[base][2] += float(a[1])
		fh = {}
		for i in ['miRNA','piRNA','siRNA']:
			os.makedirs('{0}/expr_{1}'.format(express_dir,i),mode=0o777,exist_ok=True)
			fh[i] = open("{0}/expr_{1}/{2}.expr2.xls".format(express_dir,i,base),'w')
			fh[i].writelines("sRNA id\tCount({count})\tType\tDescription\tTPM\n".format(count=exp[base][2]))
		OT =  open("{prefix}.expr2.xls".format(prefix=prefix),'w')
		OT.writelines("sRNA id\tCount({count})\tType\tDescription\tTPM\n".format(count=exp[base][2]))
		for inf in exp[base][1]:
			info = inf.split()
			if info[2] == 'mature':
				info[2] = 'miRNA'
			tpm = str(float(info[1]) / exp[base][2])
			if fh.get(info[2]):
				flag[info[2]] += 1
				fh[info[2]].writelines("%s\t%s\n"%(inf,tpm)) # smallRNA TPM = count/total
			OT.writelines("%s\t%s\n"%(inf,tpm))
		OT.close()
		for key in fh.keys():
			fh[key].close()
		for key in flag.keys():
			if flag[key] == 0:
				subprocess.check_call("rm -rf {express_dir}/expr_{key}".format(express_dir=express_dir,key=key))


















if __name__=="__main__":
	main()