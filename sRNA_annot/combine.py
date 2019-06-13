#-*-coding:utf-8-*-
#!/usr/bin/python

import os,sys
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
	purpose = opts.purpose
	samples = opts.sample
	predict_dir = opts.predict


# unknown.fa source from catalog.py(mapping to genome.fa with 10bp and not annotated by ncrna database)
# >tagid_count_chromosome position strand    example: tag1_61301_Pn7 19068351        -
def combine_unkonwn(predict_dir,samples):
	os.system(['mkdir',predict_dir])
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
			hd = name.split('_')
			chr_start_strand = hd[2].split('\t')
			end = len(seq) + int(chr_start_strand[1]) - 1
			chr_start_strand.insert(2,str(end))
			dat[seq]['map'] = "\t".join(chr_start_strand)
			dat[seq][base] = int(hd[1]) # count

	## out
	OF = open(predict_dir+"/unknown.fa")  # Take the intersection of multiple samples
	OD = open(predict_dir+"/unknown.xls") # print alignment information for each unknown tag
	file_head = "Tag\tChr\tStart\tEnd\tStrand\t"+"\t".join(key)+"\tSequence\n"
	OD.writelines(file_head)
	flag = 1
	for sq in dat.keys():
		count = 0
		OD.writelines("tag{flag}\t{map}".format(flag=flag,map=dat[sq]['map']))
		for k in key:
			if not dat[sq][k]:
				dat[sq][k] = 0
			count += dat[sq][k]
			OD.writelines("\t"+dat[sq][k])
		OD.writelines("\t"+sq+"\n")
		OF.writelines(">tag{flag}_{count}\n{seq}\n".format(flag=flag,count=count,seq=sq))
		flag += 1
	OF.close()
	OD.close()

def combine_express():
	pass













if __name__=="__main__":
	main()