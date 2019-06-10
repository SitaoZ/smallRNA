#-*-coding:utf-8-*-
import os,re,sys
import glob
from collections import defaultdict
from optparse import OptionParser


def main():
	"""
	%prog [options]
	:return:
	"""
	parser =OptionParser()
	parser.add_option('-i','--indir',help="the expression input dir")
	parser.add_option('-d','--diffList',help="the diff compare list")
	parser.add_option('-l','--log2',help="log2 for expression")
	parser.add_option('-p','--padj',help="padj")
	parser.add_option('-o','--outdir',help="output dir")

	opts,args = parser.parse_args()


	if opts.log2 == None:
		opts.log2 = 1
	if opts.padj == None:
		opts.padj = 0.1
	if opts.outdir == None:
		opts.outdir = os.getcwd()
	print(opts.log2, opts.padj, opts.outdir)
	if opts.indir == None or opts.diffList == None:
		print('\033[0;31;40m%s\033[0m' % "Warning: indir must be given\n")
		sys.exit(parser.print_help())
	indir = opts.indir
	diffList = opts.diffList
	log2 = opts.log2
	padj = opts.padj
	outdir = opts.outdir
	DESeq2(indir,diffList,log2,padj,outdir)

def DESeq2(indir,diffList,log2,padj,outdir):
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	#format reads count DESeq.xls
	#GeneID  Uniq_reads_num(42330999)        Length  Coverage        RPKM
	header = "GeneID"
	samples2 = []
	results = defaultdict()
	descs = {}
	len = {}
	for file in glob.glob(indir+"/*.expr.xls"):
		keyname = os.path.basename(file).replace('.expr.xls','')
		print("keyname",keyname)
		if re.match('^\d+',keyname):
			keyname = "X"+keyname
		samples2.append(keyname)
		with open(file,'r') as F:
			temp = F.readline()
			tabs = temp.split('\t',6)
			header += "\t%s"%(keyname)
			for line in F.readlines():
				line = line.strip()
				tabs = line.split('\t',6)
				results[tabs[0]][keyname] = tabs[1]
				len[tabs[0]] = tabs[2]
	with open(outdir+"/DESeq.xls",'w') as F:
		F.writelines(header+"\n")
		for gene in results.keys():
			F.writelines(gene+"\t")
			out = []
			for sample in samples2:
				if results[gene].get(sample):
					out.append(results[gene][sample])
				else:
					out.append('0')
			F.writelines("\t".join(out)+"\n")
	id = {}
	xls = {}
	with open(outdir+"/DESeq.xls",'r') as F:
		line2 = F.readline().strip()
		sample_name = line2.split()[1:]
		gene_num = 0
		for line in F.readlines():
			line = line.strip()
			tmp = line.split()
			id[gene_num] = tmp[0]
			gene_num += 1
			print("sample_name:",sample_name)
			for j in range(1):
				if not xls.get(sample_name[j]):
					xls[sample_name[j]] = {tmp[0]:tmp[1]}
	groups = {}
	samples = {}
	# groupA:A1,A2&groupB:B1,B2;C:c1,c2&D:d1
	all_diff = diffList.split(';')
	for diff in all_diff:
		one_diff = diff.split('&')
		for one in one_diff:
			one_group = one.split(':')
			if not groups.get(one_group[0]):
				groups[one_group[0]] = one_group[0]
				sampleaa = one_group[1].split(',')
				for sample in sampleaa:
					if re.match('^\d+',sample):
						sample = "X"+sample
					groups[one_group[0]] = {sample:sample}
					samples[sample] =sample
	# print all sample name
	cds = ''
	for key in samples:
		key = key.replace('-','.')
		cds += key+","
	cds = cs[:-1]
	print(cds)
	for diff in all_diff:
		if ',' in diff:
			sys.stderr.write()
			continue
		ga = diff.split('&')[0]
		gan = ga.split(':')[0]
		gas = ga.split(':')[1].split(',')

		gb = diff.split('&')[1]
		gbn = gb.split(':')[0]
		gbs = gb.split(':')[1].split(',')
		for sample in sample_name:
			sample_name[sample] = "C"
		for i in gas:
			if re.match('^\d+',i):
				i = "X"+i
			sample_name[i] = "A"
		for j in gbs:
			if re.match('^\d+',j):
				j = "X"+j
			sample_name[j] = "B"
		conds, tsp = '',''
		for sample in sample_name:
			conds += "\"%s\","%(sample_name[sample])
			tsp += "\'%s\',"%(sample)
		conds = conds[:-1]
		tsp = re.sub(',$','',tsp)

		OUT = open(outdir+"/"+gan+"-vs-"+gbn+".DESeq.R")
		R_script = """
		library(\"DESeq2\")
		countdata <- read.table(\"{outdir}/DESeq.xls\",skip=1)
		len <- length(countdata)
		rownames(countdata) <- countdata[,1]
		countdata <- countdata[,2:len]
		type <- c({conds})
		coldata <- data.frame(type)
		rownames(coldata) <- c({tsp})
		colnames(countdata) <- c({tsp})

		dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ type)
		dds <- DESeq(dds,quiet=TRUE)
		sizefactor <- sizeFactors(dds)
		result <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE, pAdjustMethod=\"BH\")
		write.table(result, file=\"{outdir}/tmp.output_a\", quote=FALSE, sep=\"\\t\")
		write.table(sizefactor, file=\"{outdir}/tmp.output_b\", quote=FALSE, sep=\"\\t\")";
		""".format(outdir=outdir,conds=conds,tsp=tsp)
		OUT.writelines(R_script)

		os.system("/ifs4/BC_PUB/biosoft/pipeline/Package/R-3.3.1/bin/Rscript {outdir}/{gan}-vs-{gbn}.DESeq.R".format(outdir=outdir,gan=gan,gbn=gbn))
		fh_diffexp = open("{outdir}/{gan}-vs-{gbn}_DESeq2.diffexp.xls".format(outdir=outdir,gan=gan,gbn=gbn),'w')
		fh_diffexpfilter = open("{outdir}/{gan}-vs-{gbn}_DESeq2.diffexpfilter.xls".format(outdir=outdir,gan=gan,gbn=gbn),'w')
		diff_header = "miRNA id\tExpression({gan})\tExpression({gbn})\tlog2FoldChange({gbn}/{gan})\tPvalue\tPadj\tUp/Down-Regulation\n".format(gan=gan,gbn=gbn)
		fh_diffexp.writelines(diff_header)
		fh_diffexpfilter.writelines(diff_header)
		sizefactors = {}
		with open("{outdir}/tmp.output_b".format(outdir=outdir)) as F:
			for line in F.readlines():
				a = line.strip().split()
				sizefactors[a[0]] = a[1]
		with open("{outdir}/tmp.output_a".format(outdir=outdir),'w') as F:
			a = line.strip().split()
			control_mean, treat_mean = 0,0
			ratio = 2 ** float(a[2])
			sum = 2 * float(a[1])
			control_mean = sum/(ratio+1)
			treat_mean = (sum*ratio)/(ratio+1)
			up_down = 'Up' if float(a[2]) > 0 else 'Down'
			print_out = '\t'.join([a[0],control_mean,treat_mean,a[2],a[5],a[6],up_down])
			fh_diffexp.write(print_out+"\n")
			if abs(float(a[2])) >= log2 and float(a[6]) <= padj:
				fh_diffexpfilter.writelines(print_out+'\n')
		fh_diffexp.close()
		fh_diffexpfilter.close()
		os.system('rm -rf {outdir}/tmp.output_b {outdir}/tmp.output_a {outdir}/{gan}-vs-{gbn}.DESeq.R'.format(outdir=outdir,gan=gan,gbn=gbn))
	os.system('rm -rf {outdir}/DESeq.xls'.format(outdir=outdir))


if __name__=="__main__":
	main()

