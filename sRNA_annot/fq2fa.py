#-*-coding:utf-8-*-
import os
import sys
from optparse import OptionParser
def main():
	"""
	%prog -i input_fq -o output_fa
	:return:
	"""
	parser = OptionParser()
	parser.add_option('-i','--input_fq',help='input fastq file')
	parser.add_option('-o','--output_fa',help="output fasta")
	opts,args = parser.parse_args()
	if opts.input_fq == None or opts.output_fa == None:
		sys.exit(parser.print_help())
	fq2fa(opts.input_fq,opts.output_fa)

def fq2fa(in_fq,out_fa):
	"""将fq变fa，便于cmsearch进行Rfam注释"""
	flag=0
	IN=open(in_fq,'r') if not in_fq.startswith(".gz") else gzip.open(in_fq,'rb')
	Dict = {}
	for line in IN.readlines():
		line = line.strip()
		if flag == 1:
			if line not in Dict.keys():
				Dict[line] = 1
			else:
				Dict[line] += 1
		flag += 1
		if flag == 4:
			flag = 0
	IN.close()
	line_id = 1
	with open(out_fa,'w') as OUT:
		for key in sorted(Dict.keys()):
			OUT.writelines('>tag%d_%d\n%s\n'%(line_id,Dict[key],key))





if __name__=="__main__":
	main()