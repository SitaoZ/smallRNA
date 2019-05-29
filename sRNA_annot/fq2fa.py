#-*-coding:utf-8-*-
import os
import sys
import subprocess
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
	tmp_file = out_fa.replace('.fa','.tmp')
	TMP=open(tmp_file,'w')
	for line in IN.readlines():
		if flag == 1:
			TMP.writelines(line)
		flag += 1
		if flag == 4:
			flag = 0
	IN.close()
	TMP.close()

	tmp = os.popen("sort --buffer-size=500M --temporary-directory='.' %s" % (tmp_file))
	OUT=open(out_fa,'w')
	seq,num,flag = 'na',0,1
	for line in tmp.readlines():
		if line != seq:
			if seq != 'na':
				OUT.writelines(">tag{flag}_{num}\n{seq}".format(flag=flag,num=num,seq=seq))
				flag += 1
			seq,num = line,0
		num += 1
	OUT.writelines(">tag{flag}_{num}\n{seq}".format(flag=flag, num=num, seq=seq))
	tmp.close()
	OUT.close()






if __name__=="__main__":
	main()