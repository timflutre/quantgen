#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: perform the bioinformatics aspects of genotyping-by-sequencing
# Copyright (C) 2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

# TODO:
# - check external tools are in PATH
# - add option to ignore R2 files
# - add option to give VCF of known indels
# - turn Job & co into https://docs.python.org/2/tutorial/modules.html#packages

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import sys
import os
import getopt
import time
import datetime
from subprocess import Popen, PIPE
import math
import gzip
import shlex
import warnings
import shutil
import glob
import random
import string
import stat
import xml.dom.minidom
# import sqlite3
# import numpy as np
# import scipy as sp

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Data.IUPACData import ambiguous_dna_values

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "0.1.9" # http://semver.org/


class TimUtils(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def user_input(msg):
        if sys.version_info[0] == 2:
            return raw_input(msg)
        elif sys.version_info[0] == 3:
            return input(msg)
        else:
            msg = "Python's major version should be 2 or 3"
            raise ValueError(msg)
        
    @staticmethod
    def uniq_alphanum(length):
        return "".join(random.choice(string.letters+string.digits) \
                       for i in xrange(length))
    
    
class Job(object):
    
    def __init__(self, groupId, name, cmd=None, bashFile=None):
        self.groupId = groupId
        self.name = name
        self.cmd = cmd
        self.bashFile = bashFile
        self.queue = None # set by JobGroup upon insertion
        self.duration = None # set by JobGroup upon insertion
        self.memory = None # set by JobGroup upon insertion
        self.id = None # set right after submission
        self.node = None # not used yet
        
    def submit(self):
        qsubCmd = "qsub -cwd -j y -V -q %s -N %s" % (self.queue, self.name)
        if self.duration:
            pass
        if self.memory:
            pass
        
        cmd = ""
        if self.bashFile:
            cmd += "%s %s" % (qsubCmd, self.bashFile)
        elif self.cmd:
            cmd += "echo -e '%s'" % self.cmd.encode('unicode-escape')
            cmd += " | %s" % qsubCmd
        else:
            msg = "try to submit job '%s' with neither cmd nor bash file" \
                  % self.name
            raise ValueError(msg)
        
        p = Popen(cmd, shell=True, stdout=PIPE).communicate()
        p = p[0].split()[2]
        self.id = int(p)
        
        return self.id
    
    
class JobGroup(object):
    
    def __init__(self, groupId, scheduler, queue):
        self.id = groupId
        self.scheduler = scheduler
        self.checkScheduler()
        self.queue = queue
        self.checkQueue()
        self.lJobs = [] # filled via self.insert()
        self.lJobNames = [] # filled via self.insert()
        self.lJobIds = [] # filled via self.submit()
        
    def checkScheduler(self):
        if self.scheduler not in ["SGE"]:
            msg = "unknown scheduler '%s'" % self.scheduler
            raise ValueError(msg)
        
    def checkQueue(self):
        if self.scheduler == "SGE":
            p = Popen(["qconf", "-sql"], shell=False, stdout=PIPE).communicate()
            p = p[0].split("\n")
            if self.queue not in p:
                msg = "unknown queue '%s'" % self.queue
                raise ValueError(msg)
            
    def insert(self, iJob):
        self.lJobs.append(iJob)
        self.lJobs[-1].scheduler = self.scheduler
        self.lJobs[-1].queue = self.queue
        self.lJobNames.append(iJob.name)
        
    def submit(self):
        for i in range(len(self.lJobs)):
            jobId = self.lJobs[i].submit()
            self.lJobIds.append(jobId)
            
    def getUnfinishedJobIds(self, method="oneliner"):
        if method not in ["oneliner", "xml"]:
            msg = "unknown method '%s'" % method
            raise ValueError(msg)
        
        lUnfinishedJobIds = []
        cmd = "qstat -u '%s'" % os.getlogin()
        cmd += " -q %s" % self.queue
        
        if method == "oneliner":
            cmd += " | sed 1,2d"
            cmd += " | awk '{print $1}'"
            # print(cmd) # debug
            p = Popen(cmd, shell=True, stdout=PIPE).communicate()
            p = p[0].split("\n")[:-1]
            # print(p) # debug
            lUnfinishedJobIds = [int(jobId) for jobId in p
                                 if int(jobId) in self.lJobIds]
        elif method == "xml": # http://stackoverflow.com/a/26104540/597069
            cmd += " -r -xml"
            f = os.popen(cmd)
            dom = xml.dom.minidom.parse(f)
            jobs = dom.getElementsByTagName('job_info')
            print(jobs) # debug
            for job in jobs:
                jobname = job.getElementsByTagName('JB_name')[0].childNodes[0].data
                jobown = job.getElementsByTagName('JB_owner')[0].childNodes[0].data
                jobstate = job.getElementsByTagName('state')[0].childNodes[0].data
                jobnum = job.getElementsByTagName('JB_job_number')[0].childNodes[0].data
                if jobname in self.lJobNames:
                    lUnfinishedJobIds.append(jobnum)
            print(lUnfinishedJobIds) # debug
            
        return lUnfinishedJobIds
    
    def wait(self, verbose=1):
        if verbose > 0:
            msg = "nb of jobs: %i (first=%i last=%i)" % (len(self.lJobIds),
                                                         self.lJobIds[0],
                                                         self.lJobIds[-1])
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        time.sleep(2)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        time.sleep(5)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        time.sleep(10)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        time.sleep(30)
        if len(self.getUnfinishedJobIds()) == 0:
            return
        while True:
            time.sleep(60)
            if len(self.getUnfinishedJobIds()) == 0:
                return
            
            
class JobManager(object):
    
    def __init__(self, projectId):
        self.projectId = projectId
        self.dbPath = ""
        self.groupId2group = {}
        
    def __getitem__(self, jobGroupId):
        return self.groupId2group[jobGroupId]
    
    def insert(self, iJobGroup):
        self.groupId2group[iJobGroup.id] = iJobGroup
        
        
class GbsSample(object):
    """
    Corresponds to some DNA present on a given lane. This DNA was extracted 
    from an individual which can also have some DNA on other lanes, which will 
    then be considered as other "samples".
    """
    
    def __init__(self, id):
        self.id = id # individual + "_" + flowcell + "_" + lane
        self.individual = ""
        self.generation = ""
        self.species = ""
        self.barcode = ""
        self.seqCenter = ""
        self.seqPlatform = ""
        self.seqPlatformModel = ""
        self.flowcell = ""
        self.lane = ""
        self.date = ""
        self.initFastqFile1 = ""
        self.initFastqFile2 = None
        self.dDemultiplexedFastqFiles = {} # for step 3
        self.dCleanedFastqFiles = {} # for step 4
        self.initialBamFile = "" # for step 5
        
    def __str__(self):
        txt = "id=%s" % self.id
        txt += ";individual=%s" % self.individual
        txt += ";generation=%s" % self.generation
        txt += ";species=%s" % self.species
        txt += ";barcode=%s" % self.barcode
        txt += ";seqCenter=%s" % self.seqCenter
        txt += ";seqPlatform=%s" % self.seqPlatform
        txt += ";seqPlatformModel=%s" % self.seqPlatformModel
        txt += ";flowcell=%s" % self.flowcell
        txt += ";lane=%s" % self.lane
        txt += ";date=%s" % self.date
        txt += ";initFastqFile1=%s" %self.initFastqFile1
        txt += ";initFastqFile2=%s" %self.initFastqFile2
        return txt
    
    def setDemultiplexedFastqFiles(self, pathToDir):
        lFilesR1 = glob.glob("%s/*_%s*_R1.fastq.gz" % (pathToDir,
                                                       self.individual))
        if len(lFilesR1) != 1:
            msg = "can't find R1 file for sample '%s' in '%s'" % (self.id,
                                                                  pathToDir)
            raise ValueError(msg)
        self.dDemultiplexedFastqFiles["R1"] = lFilesR1[0]
        lFilesR2 = glob.glob("%s/*_%s*_R2.fastq.gz" % (pathToDir,
                                                       self.individual))
        if len(lFilesR2) == 1:
            self.dDemultiplexedFastqFiles["R2"] = lFilesR2[0]
            
    def clean(self, adpR1, adpR2, outDir, iJobGroup, useBashScript=True):
        """
        https://cutadapt.readthedocs.org/en/stable/
        """
        cmd = "time cutadapt"
        cmd += " -a %s" % reverse_complement(str(adpR2)) # to be removed from R1 reads
        cmd += " -A %s" % reverse_complement(str(adpR1) + str(self.barcode)) # idem from R2 reads
        cmd += " -o %s/%s_clean_R1.fastq.gz" % (outDir, self.individual)
        cmd += " -p %s/%s_clean_R2.fastq.gz" % (outDir, self.individual)
        cmd += " -e 0.2" # error tolerance
        cmd += " -O 3" # min overlap length btw reads and seq passed to -a/-A
        cmd += " -m 35" # min read length
        # cmd += " -U 3" # fixed nb of bases removed from starts of R2 reads
        cmd += " -q 20,20" # quality trimming
        cmd += " --max-n 0.2"
        # cmd += " --maximum-length 150"
        cmd += " %s" % self.dDemultiplexedFastqFiles["R1"]
        if "R2" in self.dDemultiplexedFastqFiles:
            cmd += " %s" % self.dDemultiplexedFastqFiles["R2"]
        # print(cmd) # debug
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = None
        if useBashScript:
            bashFile = "job_%s_%s.bash" % (iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(iJobGroup.id, jobName, bashFile=bashFile)
        else:
            iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
    def cleanQc(self, outDir, iJobGroup):
        for Ri in ["R1", "R2"]:
            cmd = "time fastqc -o %s %s/%s_clean_%s.fastq.gz" \
                  % (outDir, outDir, self.individual, Ri)
            jobName = "stdout_%s_%s_%s" % (iJobGroup.id, self.id, Ri)
            iJob = Job(iJobGroup.id, jobName, cmd)
            iJobGroup.insert(iJob)
            
    def setCleanedFastqFiles(self, pathToDir):
        lFilesR1 = glob.glob("%s/%s_clean_R1.fastq.gz" % (pathToDir,
                                                          self.individual))
        if len(lFilesR1) != 1:
            msg = "can't find R1 file for sample '%s' in '%s'" % (self.id,
                                                                  pathToDir)
            raise ValueError(msg)
        self.dCleanedFastqFiles["R1"] = lFilesR1[0]
        lFilesR2 = glob.glob("%s/%s_clean_R2.fastq.gz" % (pathToDir,
                                                          self.individual))
        if len(lFilesR2) == 1:
            self.dCleanedFastqFiles["R2"] = lFilesR2[0]
            
    def align(self, pathToPrefixRefGenome, tmpDir, dictFile, outDir,
              iJobGroup):
        """
        http://www.htslib.org/workflow/#mapping_to_variant
        https://www.broadinstitute.org/gatk/guide/best-practices
        Need to give bashFile to Job to keep <tab> -R given to BWA
        """
        bashFile = "job_%s_%s.bash" % (iJobGroup.id, self.id)
        bashHandle = open(bashFile, "w")
        txt = "#!/usr/bin/env bash"
        txt += "\ndate"
        txt += "\noutDir=\"%s\"" % outDir
        txt += "\n\necho \"align, fixmate and sort ...\""
        bashHandle.write("%s\n" % txt)
        
        txt = "bwa mem"
        txt += " -R \'@RG"
        txt += "\tID:%s" % self.id
        txt += "\tCN:%s" % self.seqCenter
        txt += "\tDT:%s" % self.date
        txt += "\tLB:%s" % self.individual
        txt += "\tPL:%s" % self.seqPlatform
        txt += "\tPM:%s" % self.seqPlatformModel
        txt += "\tPU:%s_%s" % (self.flowcell, self.lane)
        txt += "\tSM:%s" % self.individual
        txt += "\'"
        txt += " -M %s" % pathToPrefixRefGenome
        txt += " %s" % self.dCleanedFastqFiles["R1"]
        if "R2" in self.dCleanedFastqFiles:
            txt += " %s" % self.dCleanedFastqFiles["R2"]
        bashHandle.write("%s" % txt.encode('unicode-escape'))
        
        tmpBamFile = "tmp_%s.bam" % self.id
        txt = " | samtools fixmate"
        txt += " -O bam"
        txt += " -" # stdin
        txt += " -" # stdout
        txt += " | samtools sort"
        txt += " -o ${outDir}/%s" % tmpBamFile
        txt += " -O bam"
        txt += " -T %s/tmp%s_%s" % (tmpDir, TimUtils.uniq_alphanum(5), self.id)
        txt += " -" # stdin
        bashHandle.write("%s\n" % txt)
        
        # update the header with @SQ from dictFile
        tmpHeaderFile = "tmp_%s_header.sam" % self.id
        txt = "\necho \"update header ...\""
        txt += "\ncat"
        txt += " <(samtools view -H ${outDir}/tmp_%s.bam" % self.id
        txt += " | grep -v '@SQ')"
        txt += " <(grep '@SQ' %s)" % dictFile
        txt += " > ${outDir}/%s" % tmpHeaderFile
        txt += "\ntime samtools reheader"
        txt += " ${outDir}/%s" % tmpHeaderFile
        txt += " ${outDir}/%s" % tmpBamFile
        txt += " > ${outDir}/%s.bam" % self.id
        txt += "\nrm ${outDir}/%s ${outDir}/%s" % (tmpBamFile, tmpHeaderFile)
        bashHandle.write("%s\n" % txt)
        
        # index
        txt = "\necho \"index\""
        txt += "\nsamtools index"
        txt += " ${outDir}/%s.bam" % self.id
        bashHandle.write("%s\n" % txt)
        
        # basic stats
        txt = "\necho \"flagstat\""
        txt += "\ntime samtools flagstat"
        txt += " ${outDir}/%s.bam" % self.id
        txt += " >& ${outDir}/flagstat_%s.txt" % self.id
        txt += "\ndate"
        bashHandle.write("%s\n" % txt)
        
        bashHandle.close()
        os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
        
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(iJobGroup.id, jobName, bashFile=bashFile)
        iJobGroup.insert(iJob)
        
    def setInitialBamFile(self, pathToDir):
        self.initialBamFile = "%s/%s.bam" % (pathToDir, self.id)
        
    def localRealign(self, memJvm, pathToPrefixRefGenome, outDir, iJobGroup,
                     useBashScript=True):
        cmd1 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd1 += " -T RealignerTargetCreator"
        cmd1 += " -R %s.fa" % pathToPrefixRefGenome
        cmd1 += " -I %s" % self.initialBamFile
        # cmd1 += " --known indels.vcf"
        cmd1 += " -o %s/%s.intervals" % (outDir, self.id)
        cmd2 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd2 += " -T IndelRealigner"
        cmd2 += " -R %s.fa" % pathToPrefixRefGenome
        cmd2 += " -I %s" % self.initialBamFile
        # cmd2 += " --known indels.vcf"
        cmd2 += " -targetIntervals %s/%s.intervals" % (outDir, self.id)
        cmd2 += " -o %s/%s_realigned.bam" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = None
        if useBashScript:
            bashFile = "job_%s_%s.bash" % (iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd1
            txt += "\n%s" % cmd2
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(iJobGroup.id, jobName, bashFile=bashFile)
        else:
            cmd = "%s; %s" % (cmd1, cmd2)
            iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
    def baseQualityRecalibrate(self, memJvm, pathToPrefixRefGenome, knownFile,
                               outDir, iJobGroup, useBashScript=True):
        cmd1 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd1 += " -T BaseRecalibrator"
        cmd1 += " -R %s.fa" % pathToPrefixRefGenome
        cmd1 += " -I %s/%s_realigned.bam" % (outDir, self.id)
        if knownFile != "":
            cmd1 += " --known %s" % knownFile
        else:
            cmd1 += " --run_without_dbsnp_potentially_ruining_quality"
        cmd1 += " -o %s/%s_recal.table" % (outDir, self.id)
        cmd2 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd2 += " -T PrintReads"
        cmd2 += " -R %s.fa" % pathToPrefixRefGenome
        cmd2 += " -I %s/%s_realigned.bam" % (outDir, self.id)
        cmd2 += " --BQSR %s/%s_recal.table" % (outDir, self.id)
        cmd2 += " -o %s/%s_recal.bam" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = None
        if useBashScript:
            bashFile = "job_%s_%s.bash" % (iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd1
            txt += "\n%s" % cmd2
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(iJobGroup.id, jobName, bashFile=bashFile)
        else:
            cmd = "%s; %s" % (cmd1, cmd2)
            iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
        
class GbsLane(object):
    
    def __init__(self, laneId, flowcell, number):
        self.id = laneId # flowcell identifier + "_" + lane number
        self.flowcell = flowcell
        self.number = number
        self.dSamples = {}
        self.dInitFastqFileSymlinks = {} # key(s): R1 (and R2, optional)
        
    def insert(self, iSample):
        if iSample.id in self.dSamples:
            msg = "sample '%s' already in lane '%s'" % (iSample.id, self.id)
            raise ValueError(msg)
        self.dSamples[iSample.id] = iSample
        
    def getInitFastqFiles(self):
        lFiles = []
        for sampleId,iSample in self.dSamples.items():
            if iSample.initFastqFile1 not in lFiles:
                lFiles.append(iSample.initFastqFile1)
            if iSample.initFastqFile2 and iSample.initFastqFile2 not in lFiles:
                lFiles.append(iSample.initFastqFile2)
        if len(lFiles) == 0:
            msg = "lane '%s' has no initial fastq file" % self.id
            raise ValueError(msg)
        if len(lFiles) > 2:
            msg = "lane '%s' has more than 2 initial fastq files" % self.id
            raise ValueError(msg)
        return lFiles
    
    def setInitFastqFileSymlinks(self, lInitFastqFiles, pathToDir):
        self.dInitFastqFileSymlinks["R1"] = [lInitFastqFiles[0],
                                             "%s/%s_R1.fastq.gz" \
                                             % (pathToDir, self.id)]
        if len(lInitFastqFiles) > 1:
            self.dInitFastqFileSymlinks["R2"] = [lInitFastqFiles[1],
                                                 "%s/%s_R2.fastq.gz" \
                                                 % (pathToDir, self.id)]
            
    def initQc(self, outDir, iJobGroup):
        for Ri,lFiles in self.dInitFastqFileSymlinks.items():
            cmd = "time fastqc -o %s %s" % (outDir, lFiles[1])
            jobName = "stdout_%s_%s_%s" % (iJobGroup.id, self.id, Ri)
            iJob = Job(iJobGroup.id, jobName, cmd)
            iJobGroup.insert(iJob)
            
    def saveBarcodeFile(self, format="fasta"):
        fileName = "barcodes_%s.fa" % self.id
        fileHandle = open(fileName, "w")
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sample in lSamples:
            iSample = self.dSamples[sample]
            fileHandle.write(">%s\n%s\n" % (iSample.individual,
                                            iSample.barcode))
        fileHandle.close()
        return fileName
        
    def demultiplex(self, outDir, enzyme, iJobGroup):
        if len(self.dInitFastqFileSymlinks) < 2:
            msg = "can't demultiplex (yet) lane '%s' if single-end" % self.id
            raise ValueError(msg)
        cmd = "demultiplex.py"
        cmd += " --idir %s" % os.path.dirname(self.dInitFastqFileSymlinks["R1"][1])
        cmd += " --ifq1 %s" % os.path.basename(self.dInitFastqFileSymlinks["R1"][1])
        cmd += " --ifq2 %s" % os.path.basename(self.dInitFastqFileSymlinks["R2"][1])
        cmd += " --it %s" % self.saveBarcodeFile()
        cmd += " --ofqp %s/%s" % (outDir, self.id)
        cmd += " --met %s" % "4c"
        cmd += " --re %s" % enzyme
        cmd += " --chim 1"
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
    def setDemultiplexedFastqFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setDemultiplexedFastqFiles("%s/%s" % (pathToDir, self.id))
            
    def clean(self, adpR1, adpR2, outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.clean(adpR1, adpR2, outDir, iJobGroup)
            
    def cleanQc(self, outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.cleanQc(outDir, iJobGroup)
            
    def setCleanedFastqFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setCleanedFastqFiles("%s/%s" % (pathToDir, self.id))
            
    def align(self, pathToPrefixRefGenome, tmpDir, dictFile, outDir,
              iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.align(pathToPrefixRefGenome, tmpDir, dictFile, outDir,
                          iJobGroup)
            
    def gather(self, memJvm, outDir, iJobGroup):
        cmd = "time samtools merge -f"
        cmd += " %s/%s.bam" % (outDir, self.id)
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            cmd += " %s/%s/%s.bam" % (outDir, self.id,
                                      self.dSamples[sampleId].id)
        cmd += "; samtools index"
        cmd += " %s/%s.bam" % (outDir, self.id)
        cmd += "; java -Xmx%ig -jar `which picard.jar`" % memJvm
        cmd += " CollectInsertSizeMetrics"
        cmd += " HISTOGRAM_FILE=%s/hist_insert-sizes_picard_%s.pdf" \
               % (outDir, self.id)
        cmd += " INPUT=%s/%s.bam" % (outDir, self.id)
        cmd += " OUTPUT=%s/insert-sizes_picard_%s.txt" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
    def setInitialBamFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setInitialBamFile("%s/%s" % (pathToDir, self.id))
            
    def localRealign(self, memJvm, pathToPrefixRefGenome, outDir, iJobGroup,
                     useBashScript=True):
        for sampleId,iSample in self.dSamples.items():
            iSample.localRealign(memJvm, pathToPrefixRefGenome, outDir,
                                 iJobGroup)
            
    def baseQualityRecalibrate(self, memJvm, pathToPrefixRefGenome, knownFile,
                               outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.baseQualityRecalibrate(memJvm, pathToPrefixRefGenome,
                                           knownFile, outDir, iJobGroup)
            
            
class GbsInd(object):
    
    def __init__(self, indId, flowcell, lane):
        self.id = indId
        self.flowcell = flowcell
        self.lane = lane
        self.dSamples = {}
        self.lPreprocessedBamFiles = []
        
    def insert(self, iSample):
        if iSample.id in self.dSamples:
            msg = "sample '%s' already in individual '%s'" % (iSample.id, self.id)
            raise ValueError(msg)
        self.dSamples[iSample.id] = iSample
        
    def setPreprocessedBamFiles(self, pathToDir):
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            iSample = self.dSamples[sampleId]
            self.lPreprocessedBamFiles.append("%s/%s_%s/%s_recal.bam" \
                                              % (pathToDir,
                                                 iSample.flowcell,
                                                 iSample.lane,
                                                 iSample.id))
            
    def variantCalling(self, memJvm, pathToPrefixRefGenome, knownFile, outDir,
                       iJobGroup, useBashScript=True):
        cmd = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd += " -T HaplotypeCaller"
        cmd += " -R %s.fa" % pathToPrefixRefGenome
        if len(self.lPreprocessedBamFiles) == 1:
            cmd += " -I %s" % self.lPreprocessedBamFiles[0]
        else:
            listBamsFile = "%s/bams_%s.list" % (outDir, self.id)
            listBamsHandle = open(listBamsFile, "w")
            for f in self.lPreprocessedBamFiles:
                listBamsHandle.write("%s\n" % f)
            listBamsHandle.close()
            cmd += " -I %s" % listBamsFile
        cmd += " --emitRefConfidence GVCF"
        cmd += " --variant_index_type LINEAR"
        cmd += " --variant_index_parameter 128000"
        if knownFile != "":
            cmd += " --known %s" % knownFile
        cmd += " -o %s/%s.g.vcf.gz" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = None
        if useBashScript:
            bashFile = "job_%s_%s.bash" % (iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(iJobGroup.id, jobName, bashFile=bashFile)
        else:
            cmd = "%s; %s" % (cmd1, cmd2)
            iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
        
class Gbs(object):
    
    def __init__(self):
        self.verbose = 1
        self.lSteps = []
        self.infoFile = ""
        self.scheduler = "SGE"
        self.queue = "normal.q"
        self.projectId = ""
        self.enzyme = "ApeKI"
        self.jobManager = None
        self.infoCol2idx = {"individual": None,
                            "generation": None,
                            "species": None,
                            "barcode": None,
                            "seq_center": None,
                            "seq_platform": None,
                            "seq_platform_model": None,
                            "flowcell": None,
                            "lane": None,
                            "date": None,
                            "fastq_file_R1": None,
                            "fastq_file_R2": None}
        self.dSamples = {}
        self.dLanes = {}
        self.dInds = {}
        self.lDirSteps = ["lane_qualities",
                          "demultiplexing_lanes",
                          "cleaning_samples",
                          "aligning_samples",
                          "realign_recalib",
                          "calling_variants"]
        self.adpFile = ""
        self.adapters = {}
        self.pathToPrefixRefGenome = ""
        self.dictFile = ""
        self.memJvm = 4
        self.knownFile = ""
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' performs the bioinformatics aspects of genotyping-by-sequencing.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "      --proj\tname of the project\n"
        msg += "      --schdlr\tname of the cluster scheduler (default=SGE)\n"
        msg += "      --queue\tname of the cluster queue (default=normal.q)\n"
        msg += "      --step\tstep(s) to perform (1/2/3/4, can be 1-2)\n"
        msg += "\t\t1: raw read quality (with FastQC v >= 0.11.2)\n"
        msg += "\t\t2: demultiplexing\n"
        msg += "\t\t3: cleaning (with CutAdapt v >= 1.8)\n"
        msg += "\t\t4: aligning (with BWA MEM v >= 0.7.12 and Samtools v >= 1.1)\n"
        msg += "\t\t5: local realignment and BQSR (with GATK v >= 3.3)\n"
        msg += "\t\t6: SNP and genotype calling (with GATK HC and/or BcfTools)\n"
        msg += "      --info\tpath to the 'information' file\n"
        msg += "\t\tthe file should be encoded in ASCII\n"
        msg += "\t\tthe first row should be a header with column names\n"
        msg += "\t\teach sample should have one and only one row\n"
        msg += "\t\tany two columns should be separated with one tabulation\n"
        msg += "\t\tcolumns can be in any order\n"
        msg += "\t\t12 columns are compulsory:\n"
        msg += "\t\t individual (e.g. 'Col-0', but no underscore '_')\n"
        msg += "\t\t generation (empty or 0/1, with 0 for parents)\n"
        msg += "\t\t species (e.g. 'Arabidopsis thaliana')\n"
        msg += "\t\t barcode (e.g. 'ATGG')\n"
        msg += "\t\t seq_center (e.g. 'Broad Institute', 'GenoToul', etc)\n"
        msg += "\t\t seq_platform (e.g. 'ILLUMINA', see SAM format specification)\n"
        msg += "\t\t seq_platform_model (e.g. 'HiSeq 2000')\n"
        msg += "\t\t flowcell (e.g. 'C5YMDACXX')\n"
        msg += "\t\t lane (e.g. '3')\n"
        msg += "\t\t date (e.g. '2015-01-15', see SAM format specification)\n"
        msg += "\t\t fastq_file_R1 (one per lane, gzip-compressed)\n"
        msg += "\t\t fastq_file_R2 (one per lane, gzip-compressed)\n"
        msg += "      --adp\tpath to the file containing the adapters\n"
        msg += "\t\tsame format as FastQC: name<tab>sequence\n"
        msg += "\t\tname: at least 'adpR1' (also 'adpR2' if paired-end)\n"
        msg += "\t\tsequence: from 5' (left) to 3' (right)\n"
        msg += "      --enz\tname of the restriction enzyme (default=ApeKI)\n"
        msg += "      --ref\tpath to the prefix of files for the reference genome\n"
        msg += "\t\t'/data/Atha_v2' for '/data/Atha_v2.fa', '/data/Atha_v2.bwt', etc\n"
        msg += "\t\tthese files are produced via 'bwa index ...'\n"
        msg += "      --dict\tpath to the 'dict' file (SAM header with @SQ tags)\n"
        msg += "      --jvm\tmemory given to the Java Virtual Machine (default=4, in Gb)\n"
        msg += "\t\tused in steps 4, 5 and 6 for Picard and GATK\n"
        msg += "      --known\tpath to a VCF file with known sites (e.g. from dbSNP)\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --step 1 --info info_gbs.txt\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        
        The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012).
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2015 Institut National de la Recherche Agronomique.\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:i:",
                                        ["help", "version", "verbose=",
                                         "proj=", "step=", "info=", "dict=",
                                         "schdlr=", "queue=", "enz=", "adp=",
                                         "ref=", "jvm=", "known="])
        except getopt.GetoptError as err:
            sys.stderr.write("%s\n\n" % str(err))
            self.help()
            sys.exit(2)
        for o, a in opts:
            if o == "-h" or o == "--help":
                self.help()
                sys.exit(0)
            elif o == "-V" or o == "--version":
                self.version()
                sys.exit(0)
            elif o == "-v" or o == "--verbose":
                self.verbose = int(a)
            elif o == "--proj":
                self.projectId = a
            elif o == "--schdlr":
                self.scheduler = a
            elif o == "--queue":
                self.queue = a
            elif o == "--step":
                self.lSteps = a.split("-")
            elif o == "--info":
                 self.infoFile = a
            elif o == "--enz":
                self.enzyme = a

            elif o == "--adp":
                self.adpFile = a
            elif o == "--ref":
                self.pathToPrefixRefGenome = a
            elif o == "--dict":
                 self.dictFile = a
            elif o == "--jvm":
                self.memJvm = int(a)
            elif o == "--known":
                self.knownFile = a
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if self.projectId == "":
            msg = "ERROR: missing compulsory option --proj"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "_" in self.projectId:
            msg = "ERROR: forbidden underscore '_' in project identifier '%s'" \
                  % self.projectId
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.infoFile == "":
            msg = "ERROR: missing compulsory option --info"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.infoFile):
            msg = "ERROR: can't find file %s" % self.infoFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.scheduler == "":
            msg = "ERROR: missing compulsory option --schdlr"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.scheduler == "OGE":
            self.scheduler = "SGE"
        if self.queue == "":
            msg = "ERROR: missing compulsory option --queue"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        self.jobManager = JobManager(self.projectId)
        if self.lSteps == []:
            msg = "ERROR: missing compulsory option --step"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        for i in range(len(self.lDirSteps)):
            self.lDirSteps[i] = "%s/%s_%s" % (os.getcwd(), self.projectId,
                                              self.lDirSteps[i])
        if "3" in self.lSteps:
            if self.adpFile == "":
                msg = "ERROR: missing compulsory option --adp"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not os.path.exists(self.adpFile):
                msg = "ERROR: can't find file %s" % self.adpFile
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        if "4" in self.lSteps:
            if self.dictFile == "":
                msg = "ERROR: missing compulsory option --dict"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not os.path.exists(self.dictFile):
                msg = "ERROR: can't find file %s" % self.dictFile
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if os.path.dirname(self.dictFile) == '':
                self.dictFile = "%s/%s" % (os.getcwd(), self.dictFile)
        if "4" in self.lSteps or "5" in self.lSteps or "6" in self.lSteps:
            if self.pathToPrefixRefGenome == "":
                msg = "ERROR: missing compulsory option --ref"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not os.path.exists("%s.bwt" % self.pathToPrefixRefGenome):
                msg = "ERROR: can't find file %s.bwt" % self.pathToPrefixRefGenome
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not os.path.exists("%s.fa.fai" % self.pathToPrefixRefGenome):
                msg = "ERROR: can't find file %s.fa.fai" % self.pathToPrefixRefGenome
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if os.path.dirname(self.pathToPrefixRefGenome) == '':
                self.pathToPrefixRefGenome = "%s/%s" % (os.getcwd(),
                                                        self.pathToPrefixRefGenome)
        if "5" in self.lSteps:
            if self.knownFile != "" and not os.path.exists(self.knownFile):
                msg = "ERROR: can't find file %s" % self.knownFile
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
                
                
    def loadHeaderInfoFile(self, line):
        """
        Set values to self.infoCol2idx.
        """
        try:
            line.decode('ascii')
        except UnicodeDecodeError as err:
            raise
        tokens = line.rstrip("\n").split("\t")
        if len(tokens) < 12:
            msg = "header should have at least 12 tab-separated columns"
            raise ValueError(msg)
        for idx,tok in enumerate(tokens):
            if tok in self.infoCol2idx:
                self.infoCol2idx[tok] = idx
        for infoCol,idx in self.infoCol2idx.items():
            if idx is None:
                msg = "column '%s' not found in information file" % infoCol
                raise ValueError(msg)
            
            
    def loadContentInfoFile(self, lines):
        for line in lines:
            tokens = line.rstrip("\n").split("\t")
            
            # create and fill a "GbsSample" object
            if "_" in tokens[self.infoCol2idx["individual"]]:
                msg = "underscore in individual '%s'" \
                      % tokens[self.infoCol2idx["individual"]]
                raise ValueError(msg)
            individual = tokens[self.infoCol2idx["individual"]]
            if "_" in tokens[self.infoCol2idx["flowcell"]]:
                msg = "underscore in flowcell '%s'" \
                      % tokens[self.infoCol2idx["flowcell"]]
                raise ValueError(msg)
            flowcell = tokens[self.infoCol2idx["flowcell"]]
            laneNum = int(tokens[self.infoCol2idx["lane"]])
            sampleId = "%s_%s_%i" % (individual, flowcell, laneNum)
            iSample = GbsSample(sampleId)
            iSample.individual = individual
            iSample.flowcell = flowcell
            iSample.lane = laneNum
            iSample.generation = tokens[self.infoCol2idx["generation"]]
            iSample.species = tokens[self.infoCol2idx["species"]]
            iSample.barcode = tokens[self.infoCol2idx["barcode"]]
            iSample.seqCenter = tokens[self.infoCol2idx["seq_center"]]
            iSample.seqPlatform = tokens[self.infoCol2idx["seq_platform"]]
            iSample.seqPlatformModel = tokens[self.infoCol2idx["seq_platform_model"]]
            iSample.date = tokens[self.infoCol2idx["date"]]
            if not os.path.exists(tokens[self.infoCol2idx["fastq_file_R1"]]):
                msg = "file '%s' doesn't exist" \
                      % tokens[self.infoCol2idx["fastq_file_R1"]]
                raise OSError(msg)
            if not tokens[self.infoCol2idx["fastq_file_R1"]].endswith(".gz"):
                msg = "file '%s' should be gzipped" \
                      % tokens[self.infoCol2idx["fastq_file_R1"]]
                raise ValueError(msg)
            iSample.initFastqFile1 = tokens[self.infoCol2idx["fastq_file_R1"]]
            if tokens[self.infoCol2idx["fastq_file_R2"]] != "":
                if not os.path.exists(tokens[self.infoCol2idx["fastq_file_R2"]]):
                    msg = "file '%s' doesn't exist" \
                          % tokens[self.infoCol2idx["fastq_file_R2"]]
                    raise OSError(msg)
                if not tokens[self.infoCol2idx["fastq_file_R2"]].endswith(".gz"):
                    msg = "file '%s' should be gzipped" \
                          % tokens[self.infoCol2idx["fastq_file_R2"]]
                    raise ValueError(msg)
                iSample.initFastqFile2 = tokens[self.infoCol2idx["fastq_file_R2"]]
            self.dSamples[iSample.id] = iSample
            
            laneId = "%s_%i" % (flowcell, laneNum)
            if laneId not in self.dLanes:
                self.dLanes[laneId] = GbsLane(laneId, flowcell, laneNum)
            self.dLanes[laneId].insert(iSample)
            
            if individual not in self.dInds:
                self.dInds[individual] = GbsInd(individual, flowcell, laneNum)
            self.dInds[individual].insert(iSample)
            
            
    def loadInfoFile(self):
        if self.verbose > 0:
            msg = "load information file '%s' ..." % self.infoFile
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        infoHandle = open(self.infoFile)
        
        lines = infoHandle.readlines()
        self.loadHeaderInfoFile(lines[0])
        self.loadContentInfoFile(lines[1:])
        
        infoHandle.close()
        if self.verbose > 0:
            msg = "nb of samples: %i" % len(self.dSamples)
            msg += "\nnb of lanes: %i" % len(self.dLanes)
            msg += "\nnb of individuals: %i" % len(self.dInds)
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def enterStepDir(self, dirStep):
        if os.path.isdir(dirStep):
            msg = "directory '%s' already exists" % dirStep
            warnings.warn(msg, UserWarning)
            wantRmvDir = TimUtils.user_input("Do you want to remove the directory '%s'? [y/n] " % dirStep)
            if wantRmvDir == "y":
                shutil.rmtree(dirStep)
            else:
                raise OSError("can't continue")
        os.mkdir(dirStep)
        cwd = os.getcwd()
        os.chdir(dirStep)
        return cwd


    def beginStep(self, stepNum):
        if self.verbose > 0:
            msg = "perform step %i ..." % stepNum
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        cwd = self.enterStepDir(self.lDirSteps[stepNum - 1])
        return cwd
    
    
    def setPathsToInputFiles(self, stepNum):
        if stepNum in [1, 2]:
            for laneId,iLane in self.dLanes.items():
                lFiles = iLane.getInitFastqFiles()
                iLane.setInitFastqFileSymlinks(lFiles, self.lDirSteps[0])
        if stepNum == 3:
            for laneId,iLane in self.dLanes.items():
                iLane.setDemultiplexedFastqFiles(self.lDirSteps[1])
        if stepNum == 4:
            for laneId,iLane in self.dLanes.items():
                iLane.setCleanedFastqFiles(self.lDirSteps[2])
        if stepNum == 5:
            for laneId,iLane in self.dLanes.items():
                iLane.setInitialBamFiles(self.lDirSteps[3])
        if stepNum == 6:
            for indId,iInd in self.dInds.items():
                iInd.setPreprocessedBamFiles(self.lDirSteps[4])
                
                
    def endStep(self, stepNum, cwd):
        os.chdir(cwd)
        if self.verbose > 0:
            msg = "step %i is done (see '%s/')" \
                  % (stepNum,
                     os.path.basename(self.lDirSteps[stepNum - 1]))
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def makeSymlinksToInputFastqFiles(self):
        """
        Make one symlink per lane for R1 (and R2 if necessary).
        """
        for lane,iLane in self.dLanes.items():
            for r,lFiles in iLane.dInitFastqFileSymlinks.items():
                os.symlink(lFiles[0], lFiles[1])
                
                
    def makeGroupJobId(self, prefix):
        return "%s_%s" % (prefix, TimUtils.uniq_alphanum(5))
    
    
    def launchFastqcOnInputFastqFiles(self):
        if self.verbose > 0:
            msg = "assess quality per lane ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step1-fastqc" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            iLane.initQc(self.lDirSteps[0], iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all quality jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def combineFastqcResults(self):
        # if self.verbose > 0:
        #     msg = "combine FastQC results ..."
        #     sys.stdout.write("%s\n" % msg)
        #     sys.stdout.flush()
        pass
    
    
    def step1(self):
        self.setPathsToInputFiles(1)
        cwd = self.beginStep(1)
        self.makeSymlinksToInputFastqFiles()
        self.launchFastqcOnInputFastqFiles()
        self.combineFastqcResults()
        self.endStep(1, cwd)
        
        
    def loadAdpFile(self):
        if self.verbose > 0:
            msg = "load adapter file ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        adpHandle = open(self.adpFile)
        while True:
            line = adpHandle.readline()
            if line == "":
                break
            tokens = line.rstrip("\n").split("\t")
            if len(tokens) != 2:
                msg = "adapter file '%s' should have two columns separated by a tabulation" \
                      % self.adpFile
                raise ValueError(msg)
            self.adapters[tokens[0]] = tokens[1]
        adpHandle.close()
        if "adpR1" not in self.adapters or "adpR2" not in self.adapters:
            msg = "'adpR1' or 'adpR2' missing from adapter file '%s'" % \
                  self.adpFile
            raise ValueError(msg)
        
        if self.verbose > 0:
            msg = "nb of adapters: %i" % len(self.adapters)
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def demultiplexInputFastqFiles(self):
        if self.verbose > 0:
            msg = "demultiplex samples per lane ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step2-demultiplex" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[1], laneId)
            os.mkdir(outDir)
            iLane.demultiplex(outDir, self.enzyme, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all demultiplexing jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step2(self):
        self.setPathsToInputFiles(2)
        cwd = self.beginStep(2)
        self.demultiplexInputFastqFiles()
        self.endStep(2, cwd)
        
        
    def cleanDemultiplexedFiles(self):
        if self.verbose > 0:
            msg = "clean reads per sample ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step3-clean" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[2], laneId)
            os.mkdir(outDir)
            iLane.clean(self.adapters["adpR1"], self.adapters["adpR2"],
                        outDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all cleaning jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def launchFastqcOnCleanFastqFiles(self):
        if self.verbose > 0:
            msg = "assess quality per sample ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step3-fastqc" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[2], laneId)
            iLane.cleanQc(outDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all quality jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step3(self):
        self.loadAdpFile()
        self.setPathsToInputFiles(3)
        cwd = self.beginStep(3)
        self.cleanDemultiplexedFiles()
        self.launchFastqcOnCleanFastqFiles()
        self.endStep(3, cwd)
        
        
    def alignCleanedReads(self):
        if self.verbose > 0:
            msg = "align reads per sample ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step4-align" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[3], laneId)
            os.mkdir(outDir)
            iLane.align(self.pathToPrefixRefGenome, ".", self.dictFile, outDir,
                        iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all mapping jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def gatherSamplesPerLane(self):
        if self.verbose > 0:
            msg = "gather samples per lane ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step4-gather" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            iLane.gather(self.memJvm, self.lDirSteps[3], iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all gathering jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step4(self):
        self.setPathsToInputFiles(4)
        cwd = self.beginStep(4)
        self.alignCleanedReads()
        self.gatherSamplesPerLane()
        self.endStep(4, cwd)
        
        
    def localRealignment(self):
        if self.verbose > 0:
            msg = "locally realign reads per sample ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step5-realign" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[4], laneId)
            os.mkdir(outDir)
            iLane.localRealign(self.memJvm, self.pathToPrefixRefGenome,
                               outDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all locally realignment jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def baseQualityRecalibration(self):
        if self.verbose > 0:
            msg = "recalibrate base qualities per sample ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step5-recalib" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[4], laneId)
            iLane.baseQualityRecalibrate(self.memJvm, self.pathToPrefixRefGenome,
                                         self.knownFile, outDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all recalibration jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step5(self):
        self.setPathsToInputFiles(5)
        cwd = self.beginStep(5)
        self.localRealignment()
        self.baseQualityRecalibration()
        self.endStep(5, cwd)
        
        
    def variantCallingPerIndividual(self):
        if self.verbose > 0:
            msg = "call variants per individual ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step6-indcall" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        lIndIds = self.dInds.keys()
        lIndIds.sort()
        for i,indId in enumerate(lIndIds):
            iInd = self.dInds[indId]
            iInd.variantCalling(self.memJvm, self.pathToPrefixRefGenome,
                                self.knownFile, self.lDirSteps[5], iJobGroup)
            if i == 2:
                break
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "all variant calling jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def jointGenotyping(self):
        if self.verbose > 0:
            msg = "combine variant calls across individuals ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step6-joincall" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        cmd = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % self.memJvm
        cmd += " -T GenotypeGVCFs"
        cmd += " -R %s.fa" % self.pathToPrefixRefGenome
        lIndIds = self.dInds.keys()
        lIndIds.sort()
        for i,indId in enumerate(lIndIds):
            cmd += " --variant %s/%s.g.vcf.gz" % (self.lDirSteps[5], indId)
            if i == 2:
                break
        cmd += " -o %s/%s_raw.vcf.gz" % (self.lDirSteps[5], self.projectId)
        if self.verbose > 1:
            print(cmd)
        jobName = "stdout_%s" % (iJobGroup.id)
        iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "joint calling job finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def variantRecalibration(self):
        if self.verbose > 0:
            msg = "recalibrate variant calls across individuals ..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step6-recalibcall" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        cmd1 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % self.memJvm
        cmd1 += " -T VariantRecalibrator"
        cmd1 += " -R %s.fa" % self.pathToPrefixRefGenome
        cmd1 += " -input %s/%s_raw.vcf.gz" % (self.lDirSteps[5], self.projectId)
        # cmd1 += " -resource:%s" % TODO
        cmd1 += " -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff"
        cmd1 += " -mode BOTH"
        cmd1 += " -recalFile %s/%s.recal" % (self.lDirSteps[5], self.projectId)
        cmd1 += " -tranchesFile %s/%s.tranches" % (self.lDirSteps[5], self.projectId)
        cmd1 += " -rscriptFile %s/plots_%s.R" % (self.lDirSteps[5], self.projectId)
        if self.verbose > 1:
            print(cmd1)
        cmd2 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % self.memJvm
        cmd2 += " -T ApplyRecalibration"
        cmd2 += " -R %s.fa" % self.pathToPrefixRefGenome
        cmd2 += " -input %s/%s_raw.vcf.gz" % (self.lDirSteps[5], self.projectId)
        cmd2 += " -tranchesFile %s/%s.tranches" % (self.lDirSteps[5], self.projectId)
        cmd2 += " -recalFile %s/%s.recal" % (self.lDirSteps[5], self.projectId)
        cmd2 += " -mode BOTH"
        cmd2 += " -o %s/%s.vcf.gz" % (self.lDirSteps[5], self.projectId)
        if self.verbose > 1:
            print(cmd2)
        cmd = "%s; %s" % (cmd1, cmd2)
        jobName = "stdout_%s" % (iJobGroup.id)
        iJob = Job(iJobGroup.id, jobName, cmd)
        iJobGroup.insert(iJob)
        
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait()
        
        if self.verbose > 0:
            msg = "joint calling job finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step6(self):
        self.setPathsToInputFiles(6)
        cwd = self.beginStep(6)
        self.variantCallingPerIndividual()
        self.jointGenotyping()
        self.variantRecalibration()
        self.endStep(6, cwd)
        
        
    def run(self):
        self.loadInfoFile()
        
        if "1" in self.lSteps:
            self.step1()
            
        if "2" in self.lSteps:
            self.step2()
            
        if "3" in self.lSteps:
            self.step3()
            
        if "4" in self.lSteps:
            self.step4()
            
        if "5" in self.lSteps:
            self.step5()
            
        if "6" in self.lSteps:
            self.step6()
            
            
if __name__ == "__main__":
    i = Gbs()
    
    i.setAttributesFromCmdLine()
    
    i.checkAttributes()
    
    if i.verbose > 0:
        startTime = time.time()
        msg = "START %s %s %s" % (os.path.basename(sys.argv[0]),
                                  progVersion,
                                  time.strftime("%Y-%m-%d %H:%M:%S"))
        msg += "\ncmd-line: %s" % ' '.join(sys.argv)
        msg += "\ncwd: %s" % os.getcwd()
        print(msg); sys.stdout.flush()
        
    i.run()
    
    if i.verbose > 0:
        msg = "END %s %s %s" % (os.path.basename(sys.argv[0]),
                                progVersion,
                                time.strftime("%Y-%m-%d %H:%M:%S"))
        endTime = time.time()
        runLength = datetime.timedelta(seconds=
                                       math.floor(endTime - startTime))
        msg += " (%s" % str(runLength)
        if "linux" in sys.platform:
            p = Popen(["grep", "VmHWM", "/proc/%s/status" % os.getpid()],
                      shell=False, stdout=PIPE).communicate()
            maxMem = p[0].split()[1]
            msg += "; %s kB)" % maxMem
        else:
            msg += ")"
        print(msg); sys.stdout.flush()
