#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: perform the bioinformatics aspects of genotyping-by-sequencing
# Copyright (C) 2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

# TODO:
# - see where to introduce BQSR
# - catch exception if step dir already exists, and skip step for the given lane(s)
# - add option to load pedigree (PED format)
# - check version of external programs
# - check that dates in samples file agree with SAM specification
# - add option to ignore R2 files
# - add option to give VCF of known indels (for local realign)
# - turn Job & co into https://docs.python.org/2/tutorial/modules.html#packages
# - try sgeparse https://pypi.python.org/pypi/sgeparse https://github.com/mindriot101/sgeparse

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import sys
import os
import getopt
import time
import datetime
import subprocess
import math
import gzip
import shlex
import warnings
import shutil
import glob
import string
import stat
# import xml.dom.minidom
# import sqlite3

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Data.IUPACData import ambiguous_dna_values

from pyutilstimflutre import Utils, Job, JobGroup, JobManager, Fastqc, \
    SamtoolsFlagstat

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "0.2.0" # http://semver.org/


class GbsSample(object):
    """
    A GbsSample corresponds to a unique triplet (individual,flowcell,lane),
    given that the 'individual' is the focus of the analysis.
    """
    
    def __init__(self, ind, flowcell, laneNum):
        self.individual = ind
        self.flowcell = flowcell
        self.lane = laneNum
        self.id = "%s_%s_%s" % (self.individual, self.flowcell, self.lane)
        self.species = ""
        self.library = ""
        self.barcode = ""
        self.seqCenter = ""
        self.seqPlatform = ""
        self.seqPlatformModel = ""
        self.date = ""
        self.initFastqFile1 = None
        self.initFastqFile2 = None
        self.dDemultiplexedFastqFiles = {}
        self.dCleanedFastqFiles = {}
        self.initialBamFile = ""
        
    def __str__(self):
        txt = "id=%s" % self.id
        txt += ";individual=%s" % self.individual
        txt += ";flowcell=%s" % self.flowcell
        txt += ";lane=%s" % self.lane
        txt += ";species=%s" % self.species
        txt += ";library=%s" % self.library
        txt += ";barcode=%s" % self.barcode
        txt += ";seqCenter=%s" % self.seqCenter
        txt += ";seqPlatform=%s" % self.seqPlatform
        txt += ";seqPlatformModel=%s" % self.seqPlatformModel
        txt += ";date=%s" % self.date
        txt += ";initFastqFile1=%s" %self.initFastqFile1
        txt += ";initFastqFile2=%s" %self.initFastqFile2
        return txt
    
    def setDemultiplexedFastqFiles(self, pathToDir):
        # set R1 file:
        lFilesR1 = glob.glob("%s/*_%s*_R1.fastq.gz" % (pathToDir,
                                                       self.individual))
        if len(lFilesR1) == 0:
            msg = "no demultiplexed R1 file found for sample '%s' in '%s'" % \
                  (self.id, pathToDir)
            raise ValueError(msg)
        if len(lFilesR1) > 1:
            msg = "%i demultiplexed R1 files found for sample '%s' in '%s'" % \
                  (len(lFilesR1), self.id, pathToDir)
            raise ValueError(msg)
        self.dDemultiplexedFastqFiles["R1"] = lFilesR1[0]
        # set R2 file, if necessary:
        lFilesR2 = glob.glob("%s/*_%s*_R2.fastq.gz" % (pathToDir,
                                                       self.individual))
        if len(lFilesR2) > 1:
            msg = "%i demultiplexed R2 files found for sample '%s' in '%s'" % \
                  (len(lFilesR2), self.id, pathToDir)
            raise ValueError(msg)
        if len(lFilesR2) == 1:
            self.dDemultiplexedFastqFiles["R2"] = lFilesR2[0]
            
    def clean(self, adpR1, adpR2, outDir, iJobGroup, useBashScript=True):
        """
        https://cutadapt.readthedocs.org/en/stable/
        """
        cmd = "time cutadapt"
        cmd += " -a %s" % reverse_complement(str(adpR2)) # to be removed from R1 reads
        cmd += " -A %s" % reverse_complement(str(adpR1) + str(self.barcode)) # idem from R2 reads
        cmd += " -o %s/%s_clean_R1.fastq.gz" % (outDir, self.id)
        cmd += " -p %s/%s_clean_R2.fastq.gz" % (outDir, self.id)
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
            bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(groupId=iJobGroup.id, name=jobName, bashFile=bashFile,
                       dir=outDir)
        else:
            iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
        iJobGroup.insert(iJob)
        
    def setCleanedFastqFiles(self, pathToDir):
        # set R1 file
        lFilesR1 = glob.glob("%s/%s_clean_R1.fastq.gz" % (pathToDir, self.id))
        if len(lFilesR1) == 0:
            msg = "no cleaned R1 file found for sample '%s' in '%s'" % \
                  (self.id, pathToDir)
            raise ValueError(msg)
        if len(lFilesR1) > 1:
            msg = "%i cleaned R1 files found for sample '%s' in '%s'" % \
                  (len(lFilesR1), self.id, pathToDir)
            raise ValueError(msg)
        self.dCleanedFastqFiles["R1"] = lFilesR1[0]
        # set R2 file, if necessary:
        lFilesR2 = glob.glob("%s/%s_clean_R2.fastq.gz" % (pathToDir, self.id))
        if len(lFilesR2) > 1:
            msg = "%i cleaned R2 files found for sample '%s' in '%s'" % \
                  (len(lFilesR2), self.id, pathToDir)
            raise ValueError(msg)
        if len(lFilesR2) == 1:
            self.dCleanedFastqFiles["R2"] = lFilesR2[0]
            
    def cleanQc(self, outDir, iJobGroup):
        """
        http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
        """
        for Ri in ["R1", "R2"]:
            cmd = "time fastqc -o %s %s" % (outDir,
                                            self.dCleanedFastqFiles[Ri])
            jobName = "stdout_%s_%s_%s" % (iJobGroup.id, self.id, Ri)
            iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
            iJobGroup.insert(iJob)
            
    def saveNbReadsFromFastqc(self, inDir, outHandle):
        inFile = "%s/%s_clean_R1_fastqc.zip" % (inDir, self.id)
        iFqc = Fastqc(inFile)
        txt = "%s" % self.individual
        txt += "\t%s" % self.flowcell
        txt += "\t%s" % self.lane
        txt += "\tR1"
        txt += "\t%s" % iFqc.lStats[1]["content"][3]["value"]
        outHandle.write("%s\n" % txt)
        if"R2" in self.dCleanedFastqFiles:
            inFile = "%s/%s_clean_R2_fastqc.zip" % (inDir, self.id)
            iFqc = Fastqc(inFile)
            txt = "%s" % self.individual
            txt += "\t%s" % self.flowcell
            txt += "\t%s" % self.lane
            txt += "\tR2"
            txt += "\t%s" % iFqc.lStats[1]["content"][3]["value"]
            outHandle.write("%s\n" % txt)
            
    def align(self, pathToPrefixRefGenome, tmpDir, dictFile, outDir,
              iJobGroup):
        """
        http://www.htslib.org/workflow/#mapping_to_variant
        https://www.broadinstitute.org/gatk/guide/best-practices
        Need to give bashFile to Job to keep <tab> in -R given to BWA
        """
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        bashHandle = open(bashFile, "w")
        txt = "#!/usr/bin/env bash"
        txt += "\ndate"
        txt += "\noutDir=\"%s\"" % outDir
        txt += "\n\necho \"align, fixmate and sort...\""
        bashHandle.write("%s\n" % txt)
        
        txt = "bwa mem"
        txt += " -R \'@RG"
        txt += "\tID:%s" % self.id
        txt += "\tCN:%s" % self.seqCenter
        txt += "\tDT:%s" % self.date
        txt += "\tLB:%s" % self.library
        txt += "\tPL:%s" % self.seqPlatform
        txt += "\tPM:%s" % self.seqPlatformModel
        txt += "\tPU:%s-%s.%s" % (self.flowcell, self.barcode, self.lane)
        txt += "\tSM:%s" % self.individual # see GATK's FAQ for what a sample is
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
        txt += " -T %s/tmp%s_%s" % (tmpDir, Utils.uniq_alphanum(5), self.id)
        txt += " -" # stdin
        bashHandle.write("%s\n" % txt)
        
        # update the header with @SQ from dictFile
        tmpHeaderFile = "tmp_%s_header.sam" % self.id
        txt = "\necho \"update header...\""
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
        txt += "\necho \"reads in proper pairs and primary alignments\""
        txt += "\necho -e \"%s" % self.individual
        txt += "\t%s" % self.flowcell
        txt += "\t%s" % self.lane
        txt += "\t\"$(samtools view -f 0x0002 -F 0x0100 -q 5"
        txt += " ${outDir}/%s.bam" % self.id
        txt += " | cut -f1 | sort | uniq | wc -l)"
        txt += " > ${outDir}/reads-proppaired-primaln-q5_%s.txt" % self.id
        txt += "\ndate"
        bashHandle.write("%s\n" % txt)
        
        bashHandle.close()
        os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
        
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, bashFile=bashFile,
                   dir=outDir)
        iJobGroup.insert(iJob)
        
    def saveSamtoolsFlagstat(self, inDir, outHandle):
        inFile = "%s/flagstat_%s.txt" % (inDir, self.id)
        iSf = SamtoolsFlagstat(inFile)
        txt = "%s" % self.individual
        txt += "\t%s" % self.flowcell
        txt += "\t%s" % self.lane
        txt += "\t%s" % iSf.getTxtToWrite()
        outHandle.write("%s\n" % txt)
        
    def setInitialBamFile(self, pathToDir):
        self.initialBamFile = "%s/%s.bam" % (pathToDir, self.id)
        
    def localRealign(self, memJvm, pathToPrefixRefGenome, knownIndelsFile,
                     outDir, iJobGroup, useBashScript=True):
        cmd1 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd1 += " -T RealignerTargetCreator"
        cmd1 += " -R %s.fa" % pathToPrefixRefGenome
        cmd1 += " -I %s" % self.initialBamFile
        if knownIndelsFile:
            cmd1 += " --known %s" % knownIndelsFile
        cmd1 += " -o %s/%s.intervals" % (outDir, self.id)
        cmd2 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd2 += " -T IndelRealigner"
        cmd2 += " -R %s.fa" % pathToPrefixRefGenome
        cmd2 += " -I %s" % self.initialBamFile
        if knownIndelsFile:
            cmd2 += " --known %s" % knownIndelsFile
        cmd2 += " -targetIntervals %s/%s.intervals" % (outDir, self.id)
        cmd2 += " -o %s/%s_realn.bam" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = None
        if useBashScript:
            bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd1
            txt += "\n%s" % cmd2
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(groupId=iJobGroup.id, name=jobName, bashFile=bashFile,
                       dir=outDir)
        else:
            cmd = "%s; %s" % (cmd1, cmd2)
            iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
        iJobGroup.insert(iJob)
        
    def baseQualityRecalibrate(self, memJvm, pathToPrefixRefGenome, knownFile,
                               outDir, iJobGroup, useBashScript=True):
        cmd1 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd1 += " -T BaseRecalibrator"
        cmd1 += " -R %s.fa" % pathToPrefixRefGenome
        cmd1 += " -I %s/%s_realn.bam" % (outDir, self.id)
        if knownFile != "":
            cmd1 += " --known %s" % knownFile
        else:
            cmd1 += " --run_without_dbsnp_potentially_ruining_quality"
        cmd1 += " -o %s/%s_recal.table" % (outDir, self.id)
        cmd2 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd2 += " -T PrintReads"
        cmd2 += " -R %s.fa" % pathToPrefixRefGenome
        cmd2 += " -I %s/%s_realn.bam" % (outDir, self.id)
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
        self.dir = None
        self.dSamples = {}
        self.dInitFastqFiles = {} # key(s): R1 (and R2, optional)
                                  # values: [init, symlink]
        
    def insert(self, iSample):
        if iSample.id in self.dSamples:
            msg = "sample '%s' already in lane '%s'" % (iSample.id, self.id)
            raise ValueError(msg)
        self.dSamples[iSample.id] = iSample
        
    def setInitFastqFiles(self):
        for sampleId,iSample in self.dSamples.items():
            if "R1" not in self.dInitFastqFiles:
                self.dInitFastqFiles["R1"] = [iSample.initFastqFile1]
            if iSample.initFastqFile2 and "R2" not in self.dInitFastqFiles:
                self.dInitFastqFiles["R2"] = [iSample.initFastqFile2]
                
    def makeInitFastqFileSymlinks(self):
        if len(self.dInitFastqFiles["R1"]) == 1:
            if not os.path.exists(self.dInitFastqFiles["R1"][0]):
                msg = "file '%s' doesn't exist" \
                      % self.dInitFastqFiles["R1"][0]
                raise OSError(msg)
            if not self.dInitFastqFiles["R1"][0].endswith(".gz"):
                msg = "fastq file '%s' should be gzipped" \
                      % self.dInitFastqFiles["R1"][0]
                raise ValueError(msg)
            self.dInitFastqFiles["R1"].append("%s/%s_R1.fastq.gz" %
                                              (self.dir, self.id))
        if not os.path.exists(self.dInitFastqFiles["R1"][1]):
            os.symlink(self.dInitFastqFiles["R1"][0],
                       self.dInitFastqFiles["R1"][1])
            
        if "R2" in self.dInitFastqFiles:
            if len(self.dInitFastqFiles["R2"]) == 1:
                if not os.path.exists(self.dInitFastqFiles["R2"][0]):
                    msg = "file '%s' doesn't exist" \
                          % self.dInitFastqFiles["R2"][0]
                    raise OSError(msg)
                if not self.dInitFastqFiles["R2"][0].endswith(".gz"):
                    msg = "fastq file '%s' should be gzipped" \
                          % self.dInitFastqFiles["R2"][0]
                    raise ValueError(msg)
                self.dInitFastqFiles["R2"].append("%s/%s_R2.fastq.gz" %
                                                  (self.dir, self.id))
            if not os.path.exists(self.dInitFastqFiles["R2"][1]):
                os.symlink(self.dInitFastqFiles["R2"][0],
                           self.dInitFastqFiles["R2"][1])
                
    def initQc(self, outDir, iJobGroup):
        """
        http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
        """
        for Ri,lFiles in self.dInitFastqFiles.items():
            cmd = "time fastqc -o %s %s" % (outDir, lFiles[1])
            jobName = "stdout_%s_%s_%s" % (iJobGroup.id, self.id, Ri)
            iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
            iJobGroup.insert(iJob)
            
    def saveBarcodeFile(self, outDir, format="fasta"):
        fileName = "%s/barcodes_%s.fa" % (outDir, self.id)
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
        if len(self.dInitFastqFiles) < 2:
            msg = "can't demultiplex (yet) lane '%s' if single-end" % self.id
            raise ValueError(msg)
        cmd = "demultiplex.py"
        cmd += " --idir %s" % os.path.dirname(self.dInitFastqFiles["R1"][1])
        cmd += " --ifq1 %s" % os.path.basename(self.dInitFastqFiles["R1"][1])
        cmd += " --ifq2 %s" % os.path.basename(self.dInitFastqFiles["R2"][1])
        cmd += " --it %s" % self.saveBarcodeFile(outDir, "fasta")
        cmd += " --ofqp %s/%s" % (outDir, self.id)
        cmd += " --met %s" % "4c"
        cmd += " --re %s" % enzyme
        cmd += " --chim 1"
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
        iJobGroup.insert(iJob)
        
    def setDemultiplexedFastqFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setDemultiplexedFastqFiles(pathToDir)
            
    def clean(self, adpR1, adpR2, outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.clean(adpR1, adpR2, outDir, iJobGroup)
            
    def setCleanedFastqFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setCleanedFastqFiles(pathToDir)
            
    def cleanQc(self, outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.cleanQc(outDir, iJobGroup)
            
    def saveNbReadsFromFastqc(self, inDir, outHandle):
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            iSample = self.dSamples[sampleId]
            iSample.saveNbReadsFromFastqc(inDir, outHandle)
            
    def align(self, pathToPrefixRefGenome, tmpDir, dictFile, outDir,
              iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.align(pathToPrefixRefGenome, tmpDir, dictFile, outDir,
                          iJobGroup)
            
    def setInitialBamFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setInitialBamFile(pathToDir)
            
    def gather(self, memJvm, outDir, iJobGroup):
        cmd = "echo \"merge and index...\""
        cmd += "; time samtools merge -f"
        cmd += " %s/%s.bam" % (outDir, self.id)
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            cmd += " %s" % self.dSamples[sampleId].initialBamFile
        cmd += "; samtools index"
        cmd += " %s/%s.bam" % (outDir, self.id)
        cmd += "; java -Xmx%ig -jar `which picard.jar`" % memJvm
        cmd += " CollectInsertSizeMetrics"
        cmd += " HISTOGRAM_FILE=%s/hist_insert-sizes_picard_%s.pdf" \
               % (outDir, self.id)
        cmd += " INPUT=%s/%s.bam" % (outDir, self.id)
        cmd += " OUTPUT=%s/insert-sizes_picard_%s.txt" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
        iJobGroup.insert(iJob)
        
    def saveSamtoolsFlagstat(self, inDir, outHandle):
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            iSample = self.dSamples[sampleId]
            iSample.saveSamtoolsFlagstat(inDir, outHandle)
            
    def localRealignSamples(self, memJvm, pathToPrefixRefGenome,
                            knownIndelsFile, outDir,
                            iJobGroup, useBashScript=True):
        for sampleId,iSample in self.dSamples.items():
            iSample.localRealign(memJvm, pathToPrefixRefGenome,
                                 knownIndelsFile, outDir, iJobGroup)
            
    def baseQualityRecalibrate(self, memJvm, pathToPrefixRefGenome, knownFile,
                               outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.baseQualityRecalibrate(memJvm, pathToPrefixRefGenome,
                                           knownFile, outDir, iJobGroup)
            
            
class GbsInd(object):
    
    def __init__(self, ind, flowcell, lane):
        self.id = ind
        self.flowcell = flowcell
        self.lane = lane
        self.dir = None
        self.dSamples = {}
        self.lRealignedSampleBamFiles = []
        self.realignedIndBamFile = ""
        
    def insert(self, iSample):
        if iSample.id in self.dSamples:
            msg = "sample '%s' already in individual '%s'" % (iSample.id, self.id)
            raise ValueError(msg)
        self.dSamples[iSample.id] = iSample
        
    def setRealignedSampleBamFiles(self, allLanesDir, stepDir):
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            iSample = self.dSamples[sampleId]
            self.lRealignedSampleBamFiles.append("%s/lane_%s_%s/%s/%s_realn.bam" \
                                                 % (allLanesDir,
                                                    iSample.flowcell,
                                                    iSample.lane,
                                                    stepDir,
                                                    iSample.id))
            
    def localRealign(self, memJvm, pathToPrefixRefGenome, knownIndelsFile,
                     outDir, iJobGroup, useBashScript=True):
        cmd1 = ""
        cmd2 = ""
        if len(self.lRealignedSampleBamFiles) == 1:
            pathWoExt = os.path.splitext(self.lRealignedSampleBamFiles[0])[0]
            cmd1 = "ln -s %s.bam" % pathWoExt
            cmd1 += " %s/%s_realn.bam" % (outDir, self.id)
            cmd2 = "ln -s %s.bai" % pathWoExt
            cmd2 += " %s/%s_realn.bai" % (outDir, self.id)
        else:
            cmd1 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
            cmd1 += " -T RealignerTargetCreator"
            cmd1 += " -R %s.fa" % pathToPrefixRefGenome
            for i in range(len(self.lRealignedSampleBamFiles)):
                cmd1 += " -I %s" % self.lRealignedSampleBamFiles[i]
            if knownIndelsFile:
                cmd1 += " --known %s" % knownIndelsFile
            cmd1 += " -o %s/%s.intervals" % (outDir, self.id)
            cmd2 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
            cmd2 += " -T IndelRealigner"
            cmd2 += " -R %s.fa" % pathToPrefixRefGenome
            for i in range(len(self.lRealignedSampleBamFiles)):
                cmd2 += " -I %s" % self.lRealignedSampleBamFiles[i]
            if knownIndelsFile:
                cmd2 += " --known %s" % knownIndelsFile
            cmd2 += " -targetIntervals %s/%s.intervals" % (outDir, self.id)
            cmd2 += " -o %s/%s_realn.bam" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = None
        if useBashScript:
            bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd1
            txt += "\n%s" % cmd2
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(groupId=iJobGroup.id, name=jobName, bashFile=bashFile,
                       dir=outDir)
        else:
            cmd = "%s; %s" % (cmd1, cmd2)
            iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
        iJobGroup.insert(iJob)
        
    def setRealignedIndividualBamFile(self, pathToDir):
        self.realignedIndBamFile = "%s/%s_realn.bam" % (pathToDir, self.id)
        
    def variantCalling(self, memJvm, pathToPrefixRefGenome, knownFile, outDir,
                       iJobGroup, useBashScript=True):
        """
        https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
        http://gatkforums.broadinstitute.org/discussion/comment/14337/#Comment_14337
        """
        cmd = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % memJvm
        cmd += " -T HaplotypeCaller"
        cmd += " -R %s.fa" % pathToPrefixRefGenome
        cmd += " -I %s" % self.realignedIndBamFile
        # listBamsFile = "%s/bams_%s.list" % (outDir, self.id)
        # listBamsHandle = open(listBamsFile, "w")
        # for f in self.lPreprocessedBamFiles:
        #     listBamsHandle.write("%s\n" % f)
        # listBamsHandle.close()
        # cmd += " -I %s" % listBamsFile
        cmd += " --emitRefConfidence GVCF"
        cmd += " --variant_index_type LINEAR"
        cmd += " --variant_index_parameter 128000"
        if knownFile:
            cmd += " --known %s" % knownFile
        cmd += " --activeRegionOut %s/%s_active-regions.tab" % (outDir, self.id)
        cmd += " --activityProfileOut %s/%s_activity-profiles.tab" % (outDir, self.id)
        cmd += " -o %s/%s.g.vcf.gz" % (outDir, self.id)
        cmd += " --genotyping_mode DISCOVERY"
        cmd += " --heterozygosity 0.001" # 1 het site in 100 bp (across all samples, for humans)
        # cmd += " --output_mode " ?
        cmd += " --sample_name %s" % self.id
        cmd += " --annotateNDA"
        cmd += " --activeRegionMaxSize 300" # to be increased?
        # cmd += " --input_prior " # to be used for bi-parental crosses?
        cmd += " --maxNumHaplotypesInPopulation 128"  # change? no, see link above
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        iJob = None
        if useBashScript:
            bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
            bashHandle = open(bashFile, "w")
            txt = "#!/usr/bin/env bash"
            txt += "\ndate"
            txt += "\n%s" % cmd
            txt += "\ndate"
            bashHandle.write("%s\n" % txt)
            bashHandle.close()
            os.chmod(bashFile, stat.S_IREAD | stat.S_IEXEC)
            iJob = Job(groupId=iJobGroup.id, name=jobName, bashFile=bashFile,
                       dir=outDir)
        else:
            cmd = "%s; %s" % (cmd1, cmd2)
            iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=outDir)
        iJobGroup.insert(iJob)
        
        
class Gbs(object):
    
    def __init__(self):
        self.verbose = 1
        self.lSteps = []
        self.forceRerunSteps = False
        self.samplesFile = None
        self.scheduler = "SGE"
        self.queue = "normal.q"
        self.projectId = None
        self.enzyme = "ApeKI"
        self.jobManager = None
        self.samplesCol2idx = {"individual": None,
                               "flowcell": None,
                               "lane": None,
                               "species": None,
                               "library": None,
                               "barcode": None,
                               "seq_center": None,
                               "seq_platform": None,
                               "seq_platform_model": None,
                               "date": None,
                               "fastq_file_R1": None,
                               "fastq_file_R2": None}
        self.dSamples = {}
        self.dLanes = {}
        self.dInds = {}
        self.allLanesDir = None # depends on projectId
        self.allIndsDir = None # depends on projectId
        self.lDirSteps = ["init-quality",
                          "demultiplex",
                          "clean-reads",
                          "align-reads",
                          "realign-reads-sample",
                          "realign-reads-ind",
                          "call-variants-ind",
                          "joint-genotyping",
                          "variant-recalib"]
        self.adpFile = None
        self.adapters = {}
        self.pathToPrefixRefGenome = None
        self.dictFile = None
        self.tmpDir = "."
        self.memJvm = 4 # in Gb
        self.knownIndelsFile = None
        self.knownFile = None
        
        
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
        msg += "      --step\tstep(s) to perform (1/2/3/4/..., can be 1-2-...)\n"
        msg += "\t\t1: raw read quality per lane (with FastQC v >= 0.11.2)\n"
        msg += "\t\t2: demultiplexing per lane (with demultiplex.py v >= 1.9.0\n"
        msg += "\t\t3: cleaning per sample (with CutAdapt v >= 1.8)\n"
        msg += "\t\t4: alignment per sample (with BWA MEM v >= 0.7.12 and Samtools v >= 1.1)\n"
        msg += "\t\t5: local realignment per sample (with GATK v >= 3.3)\n"
        msg += "\t\t6: local realignment per individual (with GATK v >= 3.3)\n"
        msg += "\t\t7: variant and genotype calling per individual (with GATK HaplotypeCaller v >= 3.3)\n"
        msg += "\t\t8: variant and genotype calling jointly across individuals (with GATK GenotypeGVCFs v >= 3.3)\n"
        # msg += "\t\t9: variant recalibration (with GATK v >= 3.3)\n"
        # msg += "\t\t10: genotype refinement (with GATK v >= 3.3)\n"
        # msg += "\t\t11: base quality scores recalibration (with GATK v >= 3.3)\n"
        msg += "      --samples\tpath to the 'samples' file\n"
        msg += "\t\tthe file should be encoded in ASCII\n"
        msg += "\t\tthe first row should be a header with column names\n"
        msg += "\t\teach 'sample' (see details below) should have one and only one row\n"
        msg += "\t\tany two columns should be separated with one tabulation\n"
        msg += "\t\tcolumns can be in any order\n"
        msg += "\t\trows starting by '#' are skipped\n"
        msg += "\t\t12 columns are compulsory (but there can be more):\n"
        msg += "\t\t individual (see details below, e.g. 'Col-0', but no underscore '_')\n"
        msg += "\t\t species (e.g. 'Arabidopsis thaliana')\n"
        msg += "\t\t library (e.g. can be the same as 'individual')\n"
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
        msg += "      --tmpd\tpath to a temporary directory on child nodes (default=.)\n"
        msg += "\t\te.g. it can be /tmp or /scratch\n"
        msg += "      --jvm\tmemory given to the Java Virtual Machine (default=4, in Gb)\n"
        msg += "\t\tused in steps 4, 5 and 6 for Picard and GATK\n"
        msg += "      --knowni\tpath to a VCF file with known indels (for local realignment)\n"
        msg += "      --known\tpath to a VCF file with known sites (e.g. from dbSNP)\n"
        msg += "      --force\tforce to re-run step(s)\n"
        msg += "\t\tthis removes without warning the step directory if it exists\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --step 1 --samples samples.txt\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Details:\n"
        msg += "This program aims at genotyping a set of 'individuals' using data from\n"
        msg += "a restriction-assisted DNA sequencing (RAD-seq) experiment, also known\n"
        msg += "as a genotyping-by-sequencing (GBS) experiment.\n"
        msg += "Here, by 'individual', we mean the entity which is the focus of the\n"
        msg += "study. For instance, it can be a plant variety (or a human being), or\n"
        msg += "the specific clone of a given plant variety (or a specific tumor of a\n"
        msg += "given human being), etc.\n"
        msg += "Importantly, note that the content of the 'individual' variable will\n"
        msg += "be used to set the 'SM' (sample) tag of the 'RG' (read group) header\n"
        msg += "record type of the SAM format (see http://www.htslib.org/). However,\n"
        msg += "internal to this program, the term 'sample' corresponds to the unique\n"
        msg += "triplet (individual,flowcell,lane).\n"
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
                                         "proj=", "step=", "samples=", "dict=",
                                         "schdlr=", "queue=", "enz=", "adp=",
                                         "ref=", "tmpd=", "jvm=", "knowni=",
                                         "known=", "force"])
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
                self.lSteps = sorted(a.split("-"))
            elif o == "--samples":
                 self.samplesFile = a
            elif o == "--enz":
                self.enzyme = a
            elif o == "--adp":
                self.adpFile = a
            elif o == "--ref":
                self.pathToPrefixRefGenome = a
            elif o == "--dict":
                 self.dictFile = a
            elif o == "--tmpd":
                self.tmpDir = a
            elif o == "--jvm":
                self.memJvm = int(a)
            elif o == "--knowni":
                self.knownIndelsFile = a
            elif o == "--known":
                self.knownFile = a
            elif o == "--force":
                self.forceRerunSteps = True
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if not self.projectId:
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
        self.allLanesDir = "%s/%s_all-lanes" % (os.getcwd(),
                                                self.projectId)
        if not os.path.exists(self.allLanesDir):
            os.mkdir(self.allLanesDir)
        self.allIndsDir = "%s/%s_all-individuals" % (os.getcwd(),
                                                     self.projectId)
        if not os.path.exists(self.allIndsDir):
            os.mkdir(self.allIndsDir)
        if not self.samplesFile:
            msg = "ERROR: missing compulsory option --samples"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.samplesFile):
            msg = "ERROR: can't find file %s" % self.samplesFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not self.scheduler:
            msg = "ERROR: missing compulsory option --schdlr"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.scheduler == "OGE":
            self.scheduler = "SGE"
        if not self.queue:
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
        if "1" in self.lSteps:
            if not Utils.isProgramInPath("fastqc"):
                msg = "ERROR: can't find 'fastqc' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        if "2" in self.lSteps:
            if not Utils.isProgramInPath("demultiplex.py"):
                msg = "ERROR: can't find 'demultiplex.py' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        if "3" in self.lSteps:
            if not Utils.isProgramInPath("cutadapt"):
                msg = "ERROR: can't find 'cutadapt' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not self.adpFile:
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
            if not Utils.isProgramInPath("bwa"):
                msg = "ERROR: can't find 'bwa' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not Utils.isProgramInPath("samtools"):
                msg = "ERROR: can't find 'samtools' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not Utils.isProgramInPath("picard.jar"):
                msg = "ERROR: can't find 'picard.jar' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if not self.dictFile:
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
        if "5" in self.lSteps:
            if not Utils.isProgramInPath("GenomeAnalysisTK.jar"):
                msg = "ERROR: can't find 'GenomeAnalysisTK.jar' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
            if self.knownIndelsFile and not os.path.exists(self.knownIndelsFile):
                msg = "ERROR: can't find file %s" % self.knownIndelsFile
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        if "6" in self.lSteps:
            if not Utils.isProgramInPath("GenomeAnalysisTK.jar"):
                msg = "ERROR: can't find 'GenomeAnalysisTK.jar' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        if "4" in self.lSteps or "5" in self.lSteps or "6" in self.lSteps or \
           "7" in self.lSteps or "8" in self.lSteps or "9" in self.lSteps:
            if not self.pathToPrefixRefGenome:
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
            if os.path.dirname(self.pathToPrefixRefGenome) == "":
                self.pathToPrefixRefGenome = "%s/%s" % (os.getcwd(),
                                                        self.pathToPrefixRefGenome)
                
                
    def loadHeaderSamplesFile(self, line):
        """
        Set values to self.samplesCol2idx.
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
            if tok in self.samplesCol2idx:
                self.samplesCol2idx[tok] = idx
        for samplesCol,idx in self.samplesCol2idx.items():
            if idx is None:
                msg = "column '%s' not found in samples file" % samplesCol
                raise ValueError(msg)
            
            
    def loadContentSamplesFile(self, lines):
        for line in lines:
            tokens = line.rstrip("\n").split("\t")
            
            # create and fill a "GbsSample" object
            for samplesCol in ["individual", "flowcell", "lane"]:
                if "_" in tokens[self.samplesCol2idx[samplesCol]]:
                    msg = "underscore in %s '%s'" \
                          % (samplesCol, tokens[self.samplesCol2idx[samplesCol]])
                    raise ValueError(msg)
            ind = tokens[self.samplesCol2idx["individual"]]
            flowcell = tokens[self.samplesCol2idx["flowcell"]]
            laneNum = int(tokens[self.samplesCol2idx["lane"]])
            iSample = GbsSample(ind, flowcell, laneNum)
            iSample.species = tokens[self.samplesCol2idx["species"]]
            iSample.library = tokens[self.samplesCol2idx["library"]]
            iSample.barcode = tokens[self.samplesCol2idx["barcode"]]
            iSample.seqCenter = tokens[self.samplesCol2idx["seq_center"]]
            iSample.seqPlatform = tokens[self.samplesCol2idx["seq_platform"]]
            iSample.seqPlatformModel = tokens[self.samplesCol2idx["seq_platform_model"]]
            iSample.date = tokens[self.samplesCol2idx["date"]]
            iSample.initFastqFile1 = tokens[self.samplesCol2idx["fastq_file_R1"]]
            if tokens[self.samplesCol2idx["fastq_file_R2"]] != "":
                iSample.initFastqFile2 = tokens[self.samplesCol2idx["fastq_file_R2"]]
            self.dSamples[iSample.id] = iSample
            
            laneId = "%s_%i" % (flowcell, laneNum)
            if laneId not in self.dLanes:
                self.dLanes[laneId] = GbsLane(laneId, flowcell, laneNum)
            self.dLanes[laneId].insert(iSample)
            
            if ind not in self.dInds:
                self.dInds[ind] = GbsInd(ind, flowcell, laneNum)
            self.dInds[ind].insert(iSample)
            
            
    def loadSamplesFile(self):
        if self.verbose > 0:
            msg = "load samples file '%s'..." % self.samplesFile
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        samplesHandle = open(self.samplesFile)
        
        lines = samplesHandle.readlines()
        i = 0
        while lines[i].startswith("#"):
            i += 1
        self.loadHeaderSamplesFile(lines[i])
        self.loadContentSamplesFile(lines[(i+1):])
        
        samplesHandle.close()
        
        if self.verbose > 0:
            msg = "nb of samples: %i" % len(self.dSamples)
            msg += "\nnb of lanes: %i" % len(self.dLanes)
            msg += "\nnb of individuals: %i" % len(self.dInds)
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def setupLaneDirectories(self):
        for laneId,iLane in self.dLanes.items():
            dirLane = "%s/lane_%s" % (self.allLanesDir, laneId)
            iLane.dir = dirLane
            if not os.path.exists(dirLane):
                os.mkdir(dirLane)
                
                
    def setupIndividualDirectories(self):
        for indId,iInd in self.dInds.items():
            dirInd = "%s/ind_%s" % (self.allIndsDir, indId)
            iInd.dir = dirInd
            if not os.path.exists(dirInd):
                os.mkdir(dirInd)
                
                
    def beginStep(self, stepNum):
        if self.verbose > 0:
            msg = "perform step %i..." % stepNum
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        cwd = os.getcwd()
        startTime = time.time()
        return cwd, startTime
    
    
    def makeGroupJobId(self, prefix):
        return "%s_%s" % (prefix, Utils.uniq_alphanum(5))
    
    
    def makeDir(self, dirName):
        if os.path.isdir(dirName):
            if self.forceRerunSteps:
                shutil.rmtree(dirName)
            else:
                msg = "directory '%s' already exists" % dirName
                warnings.warn(msg, UserWarning)
                msg = "Do you want to remove it? [y/n] "
                wantRmvDir = Utils.user_input(msg)
                if wantRmvDir == "y":
                    shutil.rmtree(dirName)
                else:
                    raise OSError("can't continue")
        os.mkdir(dirName)
        
        
    def endStep(self, stepNum, cwd, startTime):
        os.chdir(cwd)
        if self.verbose > 0:
            endTime = time.time()
            stepDuration = datetime.timedelta(seconds=
                                              math.floor(endTime - startTime))
            msg = "step %i is done (%s, see '%s/')" \
                  % (stepNum,
                     stepDuration,
                     os.path.basename(self.lDirSteps[stepNum - 1]))
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def launchFastqcOnInputFastqFiles(self):
        if self.verbose > 0:
            msg = "assess quality per lane..."
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
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[0])
            self.makeDir(stepDir)
            iLane.setInitFastqFiles()
            iLane.makeInitFastqFileSymlinks()
            iLane.initQc(stepDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all quality jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def combineFastqcResults(self):
        # TODO
        # if self.verbose > 0:
        #     msg = "combine FastQC results..."
        #     sys.stdout.write("%s\n" % msg)
        #     sys.stdout.flush()
        pass
    
    
    def step1(self):
        cwd, startTime = self.beginStep(1)
        self.launchFastqcOnInputFastqFiles()
        self.combineFastqcResults()
        self.endStep(1, cwd, startTime)
        
        
    def loadAdpFile(self):
        if self.verbose > 0:
            msg = "load adapter file..."
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
        if "adpR1" not in self.adapters:
            msg = "'adpR1' missing from adapter file '%s'" % \
                  self.adpFile
            raise ValueError(msg)
        firstSampleId = self.dSamples.keys()[0]
        if self.dSamples[firstSampleId].initFastqFile2 and \
           "adpR2" not in self.adapters:
            msg = "'adpR2' missing from adapter file '%s'" % \
                  self.adpFile
            raise ValueError(msg)
        
        if self.verbose > 0:
            msg = "nb of adapters: %i" % len(self.adapters)
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def demultiplexInputFastqFiles(self):
        if self.verbose > 0:
            msg = "demultiplex samples per lane..."
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
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[1])
            self.makeDir(stepDir)
            iLane.setInitFastqFiles()
            iLane.makeInitFastqFileSymlinks()
            iLane.demultiplex(stepDir, self.enzyme, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all demultiplexing jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step2(self):
        cwd, startTime = self.beginStep(2)
        self.demultiplexInputFastqFiles()
        self.endStep(2, cwd, startTime)
        
        
    def cleanDemultiplexedFiles(self):
        if self.verbose > 0:
            msg = "clean reads per sample..."
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
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[2])
            self.makeDir(stepDir)
            iLane.setDemultiplexedFastqFiles("%s/%s" % (iLane.dir,
                                                        self.lDirSteps[1]))
            iLane.clean(self.adapters["adpR1"], self.adapters["adpR2"],
                        stepDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all cleaning jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def launchFastqcOnCleanFastqFiles(self):
        if self.verbose > 0:
            msg = "assess quality per sample..."
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
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[2])
            iLane.setCleanedFastqFiles(stepDir)
            iLane.cleanQc(stepDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all quality jobs finished"
            sys.stdout.write("%s\n" % msg)


    def saveNbReadsFromFastqc(self):
        if self.verbose > 0:
            msg = "save nb of clean reads per sample..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        outFile = "%s/%s_clean-reads-per-sample.txt.gz" % (self.allLanesDir,
                                                           self.projectId)
        if os.path.exists(outFile):
            if self.forceRerunSteps:
                os.remove(outFile)
        outHandle = gzip.open(outFile, "w")
        
        txt = "ind"
        txt += "\tflowcell"
        txt += "\tlane"
        txt += "\tmate" # R1 or R2
        txt += "\tnb.clean.reads"
        outHandle.write("%s\n" % txt)
        
        lLanes = self.dLanes.keys()
        lLanes.sort()
        for laneId in lLanes:
            iLane = self.dLanes[laneId]
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[2])
            iLane.saveNbReadsFromFastqc(stepDir, outHandle)
            
        outHandle.close()
        
        
    def step3(self):
        self.loadAdpFile()
        cwd, startTime = self.beginStep(3)
        self.cleanDemultiplexedFiles()
        self.launchFastqcOnCleanFastqFiles()
        self.saveNbReadsFromFastqc()
        self.endStep(3, cwd, startTime)
        
        
    def alignCleanedReads(self):
        if self.verbose > 0:
            msg = "align reads per sample..."
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
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[3])
            self.makeDir(stepDir)
            iLane.setCleanedFastqFiles("%s/%s" % (iLane.dir, self.lDirSteps[2]))
            iLane.align(self.pathToPrefixRefGenome, self.tmpDir, self.dictFile,
                        stepDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all mapping jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def gatherSamplesPerLane(self):
        if self.verbose > 0:
            msg = "gather samples per lane..."
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
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[3])
            iLane.setInitialBamFiles(stepDir)
            iLane.gather(self.memJvm, stepDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all gathering jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def saveBasicAlignStats(self):
        if self.verbose > 0:
            msg = "save basic alignment statistics..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        outFile = "%s/%s_basic-align-stats.txt.gz" % (self.allLanesDir,
                                                      self.projectId)
        if os.path.exists(outFile):
            if self.forceRerunSteps:
                os.remove(outFile)
        outHandle = gzip.open(outFile, "w")
        
        txt = "ind"
        txt += "\tflowcell"
        txt += "\tlane"
        txt += "\t%s" % SamtoolsFlagstat.header2str()
        outHandle.write("%s\n" % txt)
        
        lLanes = self.dLanes.keys()
        lLanes.sort()
        for laneId in lLanes:
            iLane = self.dLanes[laneId]
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[3])
            iLane.saveSamtoolsFlagstat(stepDir, outHandle)
            
        outHandle.close()
        
        outFile = "%s/%s_reads-proppaired-primaln-q5.txt" % (self.allLanesDir,
                                                             self.projectId)
        if os.path.exists(outFile):
            if self.forceRerunSteps:
                os.remove(outFile)
        if os.path.exists("%s.gz" % outFile):
            if self.forceRerunSteps:
                os.remove("%s.gz" % outFile)
        outHandle = open(outFile, "w")
        txt = "ind"
        txt += "\tflowcell"
        txt += "\tlane"
        txt += "\tnb.pairs.proppaired.primaln.q5"
        outHandle.write("%s\n" % txt)
        outHandle.close()
        cmd = "cat %s/lane_*/%s/reads-proppaired-primaln-q5_*.txt" % \
              (self.allLanesDir, self.lDirSteps[3])
        # cmd += " | awk '{a[$3]+=$4} END{for(x in a){print x,a[x]}}'"
        # cmd += " | sort -k2,2n"
        cmd += " >> %s" % outFile
        cmd += "; gzip %s" % outFile
        p = subprocess.Popen(cmd, shell=True,
                             stdout=subprocess.PIPE).communicate()
        
        
    def step4(self):
        cwd, startTime = self.beginStep(4)
        self.alignCleanedReads()
        self.gatherSamplesPerLane()
        self.saveBasicAlignStats()
        self.endStep(4, cwd, startTime)
        
        
    def localRealignmentSamples(self):
        if self.verbose > 0:
            msg = "locally realign reads per sample..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step5-realign-samples" %
                                         self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[4])
            self.makeDir(stepDir)
            iLane.setInitialBamFiles("%s/%s" % (iLane.dir, self.lDirSteps[3]))
            iLane.localRealignSamples(self.memJvm, self.pathToPrefixRefGenome,
                                      self.knownIndelsFile, stepDir,
                                      iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all locally realignment samples jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step5(self):
        cwd, startTime = self.beginStep(5)
        self.localRealignmentSamples()
        self.endStep(5, cwd, startTime)
        
        
    def localRealignmentIndividuals(self):
        if self.verbose > 0:
            msg = "locally realign reads per individual..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step6-realign-ind" %
                                         self.projectId)
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
            stepDir = "%s/%s" % (iInd.dir, self.lDirSteps[5])
            self.makeDir(stepDir)
            iInd.setRealignedSampleBamFiles(self.allLanesDir,
                                            self.lDirSteps[4])
            iInd.localRealign(self.memJvm,
                              self.pathToPrefixRefGenome,
                              self.knownIndelsFile, stepDir,
                              iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all locally realignment individuals jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step6(self):
        cwd, startTime = self.beginStep(6)
        self.localRealignmentIndividuals()
        self.endStep(6, cwd, startTime)
        
        
    def variantCallingPerIndividual(self):
        if self.verbose > 0:
            msg = "call variants per individual..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step7-varcall-ind" % self.projectId)
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
            stepDir = "%s/%s" % (iInd.dir, self.lDirSteps[6])
            self.makeDir(stepDir)
            iInd.setRealignedIndividualBamFile("%s/%s" % (iInd.dir,
                                                          self.lDirSteps[5]))
            iInd.variantCalling(self.memJvm, self.pathToPrefixRefGenome,
                                self.knownFile, stepDir, iJobGroup)
            
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all variant calling per ind jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step7(self):
        cwd, startTime = self.beginStep(7)
        self.variantCallingPerIndividual()
        self.endStep(7, cwd, startTime)
        
        
    def jointGenotyping(self):
        if self.verbose > 0:
            msg = "combine variant calls across individuals..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step8-joincall" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        stepDir = "%s/%s" % (self.allIndsDir, self.lDirSteps[7])
        self.makeDir(stepDir)
        cmd = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % self.memJvm
        cmd += " -T GenotypeGVCFs"
        cmd += " -R %s.fa" % self.pathToPrefixRefGenome
        lIndIds = self.dInds.keys()
        lIndIds.sort()
        for i,indId in enumerate(lIndIds):
            iInd = self.dInds[indId]
            prevStepDir = "%s/%s" % (iInd.dir, self.lDirSteps[6])
            cmd += " --variant %s/%s.g.vcf.gz" % (prevStepDir, indId)
        cmd += " -o %s/%s_raw.vcf.gz" % (stepDir, self.projectId)
        if self.verbose > 1:
            print(cmd)
        jobName = "stdout_%s" % (iJobGroup.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=stepDir)
        iJobGroup.insert(iJob)
        
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "joint calling job finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step8(self):
        cwd, startTime = self.beginStep(8)
        self.jointGenotyping()
        self.endStep(8, cwd, startTime)
        
        
    def variantRecalibration(self):
        if self.verbose > 0:
            msg = "recalibrate variant calls across individuals..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step9-recalibcall" % self.projectId)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.scheduler, self.queue)
        self.jobManager.insert(iJobGroup)
        
        stepDir = "%s/%s" % (self.allIndsDir, self.lDirSteps[8])
        self.makeDir(stepDir)
        cmd1 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % self.memJvm
        cmd1 += " -T VariantRecalibrator"
        cmd1 += " -R %s.fa" % self.pathToPrefixRefGenome
        prevStepDir = "%s/%s" % (self.allIndsDir, self.lDirSteps[8])
        cmd1 += " -input %s/%s_raw.vcf.gz" % (prevStepDir, self.projectId)
        # cmd1 += " -resource:%s" % TODO
        cmd1 += " -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff"
        cmd1 += " -mode BOTH"
        cmd1 += " -recalFile %s/%s.recal" % (stepDir, self.projectId)
        cmd1 += " -tranchesFile %s/%s.tranches" % (stepDir, self.projectId)
        cmd1 += " -rscriptFile %s/plots_%s.R" % (stepDir, self.projectId)
        if self.verbose > 1:
            print(cmd1)
        cmd2 = "java -Xmx%ig -jar `which GenomeAnalysisTK.jar`" % self.memJvm
        cmd2 += " -T ApplyRecalibration"
        cmd2 += " -R %s.fa" % self.pathToPrefixRefGenome
        cmd2 += " -input %s/%s_raw.vcf.gz" % (prevStepDir, self.projectId)
        cmd2 += " -tranchesFile %s/%s.tranches" % (stepDir, self.projectId)
        cmd2 += " -recalFile %s/%s.recal" % (stepDir, self.projectId)
        cmd2 += " -mode BOTH"
        cmd2 += " -o %s/%s.vcf.gz" % (stepDir, self.projectId)
        if self.verbose > 1:
            print(cmd2)
        cmd = "%s; %s" % (cmd1, cmd2)
        jobName = "stdout_%s" % (iJobGroup.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd, dir=stepDir)
        iJobGroup.insert(iJob)
        
        self.jobManager[iJobGroup.id].submit()
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "variant recalibration job finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step9(self):
        cwd, startTime = self.beginStep(9)
        self.variantRecalibration()
        self.endStep(9, cwd, startTime)
        
        
    def step10(self):
        cwd, startTime = self.beginStep(10)
        # self.genotypeRefinement()
        self.endStep(10, cwd, startTime)
        
        
    def baseQualityRecalibration(self):
        if self.verbose > 0:
            msg = "recalibrate base qualities per sample..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step11-recalib" % self.projectId)
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
        self.jobManager[iJobGroup.id].wait(self.verbose)
        
        if self.verbose > 0:
            msg = "all recalibration jobs finished"
            sys.stdout.write("%s\n" % msg)
            
            
    def step11(self):
        cwd, startTime = self.beginStep(11)
        self.baseQualityRecalibration()
        self.endStep(11, cwd, startTime)
        
        
    def run(self):
        self.loadSamplesFile()
        self.setupLaneDirectories()
        self.setupIndividualDirectories()
        
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
            
        if "7" in self.lSteps:
            self.step7()
            
        if "8" in self.lSteps:
            self.step8()
            
        if "9" in self.lSteps:
            self.step9()
            
        if "10" in self.lSteps:
            self.step10()
            
            
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
            p = subprocess.Popen(["grep", "VmHWM", "/proc/%s/status" % os.getpid()],
                      shell=False, stdout=subprocess.PIPE).communicate()
            maxMem = p[0].split()[1]
            msg += "; %s kB)" % maxMem
        else:
            msg += ")"
        print(msg); sys.stdout.flush()
