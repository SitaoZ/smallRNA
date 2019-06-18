#-*-coding:utf-8-*-
#!/usr/bin/python

import os,sys
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
	parser.add_option('-k','--purpose',help="purpose, 0:unknown 1:express 2:express_aasra")
	parser.add_option('-f','--sample',help="prefix of each sample annotation,separated bu comma")
	parser.add_option('-p','--predict',help="prediction outdir")
	parser.add_option('-d','--dir',help="expression outdir for degs")

	opts,args = parser.parse_args()
	if opts.purpose == None or opts.sample == None or opts.predict ==None :
		print('\033[0;31;40m%s\033[0m' % "Warning: purpose, samples and predict_dir must be given\n")
		sys.exit(parser.print_help())
	purpose = opts.purpose
	samples = opts.sample
	predict_dir = opts.predict
	combine_unkonwn(predict_dir,samples)



def combine_unkonwn(predict_dir,samples):
	"""
	combine all sample unknown read which map multiply(10 locites) to genome
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
		base = os.path.basename(prefix)
		key.append(base)
		fasta_obj = Fasta(prefix+'.unknown.fa')
		for record in fasta_obj:
			name = record.id
			seq = record.seq
			hd = name.split('_') # name formt #>tagid_count_chromosome position strand
			chr_start_strand = hd[2].split('\t')
			end = len(seq) + int(chr_start_strand[1]) - 1
			chr_start_strand.insert(2,str(end)) # add end
			dat[seq]['map'] = "\t".join(chr_start_strand)
			dat[seq][base] = int(hd[1]) # count

	## out
	OF = open(predict_dir+"/unknown.fa",'w')  # Take the intersection of multiple samples
	OD = open(predict_dir+"/unknown.xls",'w') # print alignment information for each unknown tag
	OD.writelines("Tag\tChr\tStart\tEnd\tStrand\t"+"\t".join(key)+"\tSequence\n") # key 顺序固定
	flag = 1
	for sq in dat.keys():
		count = 0
		OD.writelines("tag{flag}\t{map}".format(flag=flag,map=dat[sq]['map'])) #tag chr start end strand
		for k in key:
			if not dat[sq].get(k):
				dat[sq][k] = 0
			count += dat[sq][k] # 每一个序列特定下，把所有的个体的该序列的条数累加
			OD.writelines("\t"+str(dat[sq][k])) # write count of every sample
		OD.writelines("\t"+sq+"\n") # sequence
		OF.writelines(">tag{flag}_{count}\n{seq}\n".format(flag=flag,count=count,seq=sq)) # fasta
		flag += 1
	OF.close()
	OD.close()

def predict_sRNA(file_handle,exp,tag,key,typ):
	pass

def combine_express(predict_dir,samples,express_dir):
	"""
	combine known and unknown expression file
	:param predict_dir: predict unknown out dir
	:param samples: samples prefix partition comma
	:param express_dir: expression output dir
	:return:none
	"""
	tag = {}
	keys = []
	exp = defaultdict(dict)
	if os.path.exists(predict_dir+"/unknown.xls"):
		with open(predict_dir+"/unknown.xls",'r') as F:
			for line in F.readlines():
				a = line.strip().split()
				if a[0] == 'Tag':
					key = a[5:]
				else:
					tag[a[0]] = a[5:]

	if os.path.exists(predict_dir+"/miRNA/predicted_miRNA.xls"):
		predicted_sRNA(predict_dir+"/miRNA/predicted_miRNA.xls",exp,tag,key,'miRNA')
	if os.path.exists(predict_dir + "/miRNA/predicted_piRNA.xls"):
		predicted_sRNA(predict_dir + "/miRNA/predicted_piRNA.xls", exp, tag, key, 'piRNA')
	if os.path.exists(predict_dir + "/miRNA/predicted_siRNA.xls"):
		predicted_sRNA(predict_dir + "/miRNA/predicted_siRNA.xls", exp, tag, key, 'siRNA')
	flag= defaultdict(dict)
	for prefix in samples.split(','):
		base = os.path.basename(prefix)
		with open(prefix+".catalog.xls",'r') as F:
			for line in F.readlines():
				a = line.strip().split()
				if a[0] == 'sRNA id' or a[0] == 'unmap':
					exp[key][1] += line #save all samples expression
					exp[key][2] += float(a[1])
		fh = {}
		for i in ['miRNA','piRNA','siRNA']:
			subprocess.check_call('mkdir -p {a}/expr_{b}'.format(a=express_dir,b=i))
			fh[i] = open("{1}/expr_{2}/{3}.expr.xls")
			fh[i].writelines("sRNA id\tCount({count})\tType\tDescription\tTPM\n".format(count=exp[key][2]))
		OT =  open("{prefix}.expr.xls".format(prefix=prefix))
		OT.writelines("sRNA id\tCount({count})\tType\tDescription\tTPM\n".format(count=exp[key][2]))
		for inf in exp[key][1].split("\n"):
			info = inf.split('\t')
			if info[2] == 'mature':
				info[2] = 'miRNA'
			if fh.get(info[2]):
				flag[info[2]] += 1
			fh[i].writelines("%s\t%.2f\n"%(inf,float(info[1])/exp[key][2]))
		OT.writelines("%s\t%.2f\n"%(inf,float(info[1])/exp[key][2]))
		OT.close()
		for key in fh.keys():
			fh[key].close()
		for key in flag.keys():
			if flag[key] == 0:
				subprocess.check_call("rm -rf {express_dir}/expr_{key}".format(express_dir=express_dir,key=key))


















if __name__=="__main__":
	main()