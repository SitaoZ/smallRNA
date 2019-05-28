#-*-coding:utf-8-*-
import os
import sys
from multiprocessing import Pool

def main():
	"""
	%prog -i input_fq -o output_fa
	:return:
	"""
	parser = OptionParser()
	parser.add_option('-i','--input_fq',help='input fastq file')
	parser.add_option('-o','--output_fq',help="output fastq with CCCCC GGGGG")
	opts,args = parser.parse_args()
	if opts.input_fq == None or opts.output_fq == None:
		sys.exit(parser.print_help())
	aasra_fq(opts.input_fq,opts.output_fq)

def aasra_fq(in_fq,out_fq):
	IN = open(in_fq, 'r') if not in_fq.startswith(".gz") else gzip.open(in_fq, 'rb')
	OUT = open(out_fq,'w')
	flag = 0
	for line in IN.readlines():
		line = line.strip()
		if flag == 1:
			out = "C"*5+line+"G"*5
			OUT.writelines(out+"\n")
		elif flag == 3:
			out = 'i'*5+line+'i'*5
			OUT.writelines(out + "\n")
			flag = 0
		else:
			out = line+"\n"
			OUT.writelines(out + "\n")
		flag += 1
	IN.close()
	OUT.close()




