#-*-coding:utf-8-*-
import re,os,sys
from optparse import OptionParser



def main():
    """
    %prog [options]
    :return:
    """
    parser = OptionParser()
    parser.add_option('-f','--fastq',help="the fastq file")
    parser.add_option('-k','--sample_id',help="sample name")
    parser.add_option('-c','--clean',help="clean read number")
    parser.add_option('-a','--adaptor',help="a adaptor type")
    parser.add_option('-q','--quality_system',help="quality system :illumina/sanger")
    parser.add_option('-s','--shell',help="shell outdir")
    parser.add_option('-r','--read_outdir',help="read outdir")
    opts,args = parser.parse_args()
    if opts.fastq == None or opts.sample_id == None or opts.clean == None or opts.adaptor == None or opts.quality_system == None or opts.shell == None or opts.read_outdir == None:
        print('\033[0;31;40m%s\033[0m' % "Warning: all options must be given\n")
        sys.exit(parser.print_help())
    fastq = opts.fastq
    sample_id = opts.sample_id
    clean = opts.clean
    adaptor = opts.adaptor
    quality = opts.quality_system
    shell_dir = opts.shell
    read_dir = opts.read_outdir
    filter(fastq,sample_id,clean,adaptor,quality,shell_dir,read_dir)




def filter(fastq,sample_name,clean,adaptor,quality,shell_dir,read_dir):
    path = os.path.dirname(fastq)
    shell = ''
    adt = []
    adaptor_list = ['GTTCAGAGTTCTACAGTCCGACGATC', 'TCGTATGCCGTCTTCTGCTTGT', 'ATCTCGTATGCCGTCTTCTGCTTG', 'TGGAATTCTCGGGTGCCAAGG', 'AGATCGGAAGAGCACACGTCT']
    if os.path.exists(path,'1.fqcheck'):
        shell += 'tile=`perl $Bin/findNtile.pl -fq1 $opts{f} -fc1 $path/1.fqcheck` \n'
    else:
        shell += 'tile=9999\n'
    if re.match(r'(\w+)\,(\w+)',adaptor):
        adt = re.match(r'(\w+)\,(\w+)',adaptor).group(0)
    elif adaptor == 5:
        adt = ['GAACGACATGGCTACGATCCGACTT', 'AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA']
    elif adaptor_list[adaptor]:
        adt = adaptor_list[0],adaptor_list[adaptor]
    clean = clean*1024
    shell += "soapnuke filtersRNA -f {fastq} -F Raw.{sample} -3 {adt1} -5 {adt2} --tile $tile -o {read_dir} -x {sample} -Q {quality} -q -c {clean}\n".format(fastq=fastq,sample=sample_name,adt1=adt[1],adt2=adt[0],read_dir=read_dir,quality=quality,clean=clean)
    os.popen('mkdir %s %s '%(shell_dir,read_dir))
    with open(shell_dir+"/"+sample_name+"_filter.sh",'w') as F:
        F.writelines(shell)
if __name__=="__main__":
    main()