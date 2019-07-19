#coding:utf-8
#!/usr/bin/python
import sys,os
import datetime
from biosoftware import blat,bowtie2,RepeatMasker
from subprocess import Popen, PIPE


class database(object):
    """ the first step for smallRNA """
    def __init__(self,workdir,genome,transcript,verbose = False):
        if not os.path.exists(workdir):
            os.makedirs(workdir,0o755,True)
            self.workdir = workdir

        if os.path.exists(genome):
            self.genome = genome
        else:
            print("check the genome path")
            sys.exit(1)
        if os.path.exists(transcript):
            self.transcript = transcript
        else:
            print("check the transcript")
            sys.exit(1)
        self.verbose = verbose
        self.DataBase = {"snoRNA":"/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_smallRNA/RNA_smallRNA_2016a/database/snoRNA/ath/ath.fa",
                         "mirBase":"/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_smallRNA/RNA_smallRNA_version2.2_beta/mirbase/21/mirbase.txt",
                         "organisms":"/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_smallRNA/RNA_smallRNA_version2.2_beta/mirbase/21/organisms.txt"}


    def _log(self, content, showTime=True):
        ##### Display text to defined stream, add timestamp
        logStream = sys.stdout
        if isinstance(content, list):
            modifText = " ".join(content)
        else:
            modifText = str(content)

        if showTime == True:
            modifText += "    " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if self.verbose == True:
            logStream.write(modifText + "\n")
        else:
            return modifText
    def runSysCommandSubProcess(self, command, showInLog=True, justTest=False):
        if showInLog == True:
            self._log(command, showTime=False)
            if justTest == True:
                return
        p = Popen(" ".join(command), shell=True, stderr=PIPE)
        return p.stderr.read().decode()

    def species(self,spe):
        organisms = self.DataBase['organisms']
        organisms_dict = {}
        with open(organisms,'r') as F:
            for line in F.readlines():
                tmp = line.strip().split()
                organism,division,name,tree,NCBI_taxid = tmp[0:5]
                organisms_dict[organism] = NCBI_taxid
        if spe in organisms_dict.keys():
            return True
        else:
            return False


    def genome_index(self):
        d_genome = self.workdir + "/d_genome"
        if not os.path.exists(d_genome):
            os.makedirs(d_genome, 0o755, True)
        new_refer_path = d_genome + "/genome.fa"
        if not os.path.islink(new_refer_path):
            os.symlink(self.genome,new_refer_path,target_is_directory=False)

        prefix_refer = os.path.join(d_genome,"genome")

        self._log('Start genome index task', showTime=True)
        cmd = [bowtie2(), '-q -f', new_refer_path, prefix_refer]
        self.runSysCommandSubProcess(cmd, showInLog=True, justTest=False)
        self._log('Finish genome index task', showTime=True)


    def genome_repeatMasker(self):
        d_genome = self.workdir + "/d_genome"
        new_refer_path = d_genome + "/genome.fa"
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir+"/genome_annot")
        repeatmasker_path = RepeatMasker()
        repeatmasker_outpath = self.workdir+"/genome_annot"
        self._log('Start genome repeat task', showTime=True)
        cmd1 = [repeatmasker_path, new_refer_path, '-dir', repeatmasker_outpath]
        self.runSysCommandSubProcess(cmd1, showInLog=True, justTest=False)
        cmd2 = ["mv", repeatmasker_outpath+"/*out", repeatmasker_outpath+"/genome.repeat"]
        self.runSysCommandSubProcess(cmd2, showInLog=True, justTest=False)
        self._log('Finish genome repeat task',showTime=True)


    def genome_region(self):
        """ 获取基因组区域 """
        from database.get_psl import psl, psl_filter, genome_region
        cwd = self.workdir + "/genome_annot"
        if not os.path.exists(cwd):
            os.makedirs(cwd)
        pslfile = cwd + "/genome.psl"
        self._log('Start genome region task', showTime=True)
        # step1 blat 比对
        cmd = psl(blat(),self.genome,self.transcript,pslfile)
        self.runSysCommandSubProcess(cmd,showInLog=True,justTest=False)
        # step 2 psl 过滤
        psl_tmpfile = psl_filter(pslfile,cwd)
        # step 3 genome region result
        genome_region_file = genome_region(psl_tmpfile,cwd)
        self._log('Finish genome region task', showTime=True)




    def ncgb(self):
        """non conding rna from genebank, ftp://ftp.ncbi.nih.gov/genbank/"""
        pass

    def mRNA(self):
        pass

    def mirBase(self,spe):
        """
        从官方的mirBase数据库提取micBase已知的 mature.fa miRNA.fa
        ftp://mirbase.org/pub/mirbase/
        """
        self._log('Start Extract mature miRNA des in mirBase ', showTime=True)
        outdir = self.workdir + "/d_sncrna"
        if not os.path.exists(outdir):
            os.makedirs(outdir,0o755,True)
        if self.species(spe):
            self._log('the species is in organise ', showTime=True)
        else:
            self._log('the species is not in organise ', showTime=True)
            sys.exit('1')

        from database.mirbase import known
        mirbase = self.DataBase["mirBase"]

        mature_fa, miRNA_fa, miRNA_des = known(spe,mirbase,outdir)
        self._log('Finish Extract mature miRNA des in mirBase ', showTime=True)
        return mature_fa, miRNA_fa, miRNA_des

    def run_aasra(self,spe):
        """ extract sncRNA 提取所有的小的非编码RNA，包括snoRNA（小核仁RNA），mature， hairpin """
        self._log('Start merge snoRNA,mature,hairpin ', showTime=True)
        outdir = self.workdir+"/d_sncrna"
        snoRNA = self.DataBase["snoRNA"]
        mature = outdir+"/mature.fa"
        hairpin = outdir+"/miRNA.fa"
        dbs=[snoRNA,mature,hairpin]
        from sRNA_annot.aasra import Aasra
        Aasra(spe,dbs,outdir)
        self._log('Start build index for sncrna ', showTime=True)
        cmd = [bowtie2(), '-q -f', outdir+"/sncrna.fa", outdir+"/sncrna"]
        self.runSysCommandSubProcess(cmd,showInLog=True,justTest=False)
        self._log('Finish build index for sncrna ', showTime=True)
        self._log('Finish merge snoRNA,mature,hairpin ', showTime=True)
        return outdir + "/sncrna.fa"




    def Rfam(self):
        """
        ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.0/
        """
        return "cm_path"
    def piRNA(self):
        """
        piRBase  http://www.regulatoryrna.org/database/piRNA/download.html
        piRNABank  http://pirnabank.ibab.ac.in/request.html
        """
        pass
    def snoRNA(self):
        """
        RNAdb  http://research.imb.uq.edu.au/rnadb/rnadb2_archive.htm
        snoRNA Database  http://lowelab.ucsc.edu/snoRNAdb/
        Plant snoRNA Database  http://bioinf.scri.sari.ac.uk/cgi-bin/plant_snorna/home
        """
        pass



    def target(self):
        """
        http://mirdb.org/miRDB/download.html
        http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index
        """
        pass

    def taxonomy(self):
        """
        NCBI /ifs1/pub/database/ftp.ncbi.nih.gov/pub/taxonomy/20160229/taxdump.tar.gz
        """
        pass


if __name__ == "__main__":
    build_database =database("/nascngb/gccnt2/ST_BI/zhusitao/05.test/small/00.data/test","/nascngb/gccnt2/ST_BI/zhusitao/05.test/small/00.data/small.fa","/hwfssz5/ST_BIGDATA/USER/zhusitao/Project/micro/00.data/v2/Piper_nigrum.annotation.cds",verbose = True)
    build_database.species("ath")
    build_database.genome_index()
    build_database.genome_repeatMasker()
    build_database.genome_region()
    build_database.mirBase('ath')
    build_database.run_aasra("ath")

