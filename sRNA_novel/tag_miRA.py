#-*-coding:utf-8-*-
import sys
from util.fasta import Fasta
if sys.argv !=3:
	print ('usage:python3 get_miRA.py <unknown> <tag> <count>')
	sys.exit(1)

unknon_fasta = Fasta(sys.argv[1]).fastaDict()
OUT = open(sys.argv[2],'w')

for header,seq in unknon_fasta.items():
	hd = header.split("_")
	if had[1] > sys.argv[3]:
		OUT.writelines(">{header}_{count}\n{seq}\n".format(header=hd[0],count=hd[1]),seq=seq)

OUT.close()


