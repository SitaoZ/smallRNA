#-*-coding:utf-8-*-
import sys,os


def bowtie(bowtie_software_path,genome,tag_fa,outdir):
	"""
	mapping tag.fa to genome
	:return: sam file
	"""
	para = "-f -l 16 -p 6 "
	sam=os.path.join(outdir,'mipa.sam')
	shell = "{bowtie} {para} {reference} {tag_fa} {sam}".format(bowtie=bowtie_software_path,\
																para = para,\
																reference = genome,\
																sam=sam)
	shell += "{samtools} view -hS -F 4 {sam} | {samtools} sort -O SAM -o {mipa_sorted.sam}".format()
	shell += "{miRA} full {miRA_config} {mipa_sorted.sam} {reference} {output} ".format()


