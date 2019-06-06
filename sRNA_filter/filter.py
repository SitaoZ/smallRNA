import os,sys
import datetime
from subprocess import Popen, PIPE, check_output
import platform
"""
Author: Sitao Zhu
Date  : 20190101
Desc  : a class for RNA fqstq filter,python version 3.7
"""
class filter(object):
    def __init__(self,fastq,outdir,soapnuke_path,qualSys,verbose = False):
        self.fastq = fastq
        self.outdir = outdir
        self.soapnuke = soapnuke_path
        self.qualSys = qualSys
        self.platform = sys.platform
        self.adapter = '-3 AGATCGGAAGAGCACACGTCT -5 GTTCAGAGTTCTACAGTCCGACGATC'
        self.verbose = verbose


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

    def runSysCommand(self, command, showInLog=True, justTest=False):
        if showInLog == True:
            self._log(command, showTime=False)
            if justTest == True:
                return
        try:
            # returnVal = os.popen(" ".join(command))
            p = Popen(" ".join(command), shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
        except Exception as e:
            pass
        return 0

    def runSysCommandSubProcess(self, command, showInLog=True, justTest=False):
        if showInLog == True:
            self._log(command, showTime=False)
            if justTest == True:
                return
        p = Popen(" ".join(command), shell=True, stderr=PIPE)
        return p.stderr.read().decode()

    def _runSOAPnuke(self,raw_pre,clean_pre):
        self.raw_pre = raw_pre
        self.clean_pre = clean_pre
        self._log('Start soapnuke filtersRNA task',showTime=True)
        cmd = [self.soapnuke,'filtersRNA','-f',self.fastq,'-F',self.raw_pre,self.adapter,'-o',self.outdir,'-x',self.clean_pre,'-Q',self.qualSys,'-q -c 11821.2725806991']
        self.runSysCommandSubProcess(cmd,showInLog=True,justTest=True)
        self._log('Finish soapnuke filtersRNA task ',showTime=True)

if __name__=="__main__":
    a = filter('/hwfssz5/ST_BIGDATA/USER/zhusitao/Project/micro/00.data/rawdata/D1_raw.fq.gz','/nascngb/gccnt2/ST_BI/zhusitao/05.test/small/01.flow/process/Shell/D1/test/test','/ifs4/BC_PUB/biosoft/pipeline/Package/SOAPnuke/SOAPnuke1.5.0','2',verbose=True)
    a._runSOAPnuke('raw','clean')
