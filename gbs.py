#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: perform the computational aspects of genotyping-by-sequencing
# Copyright (C) 2015-2017 Institut National de la Recherche Agronomique
# License: AGPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

# TODO:
# - refactor variantFiltering()
# - use env var for Java jar (e.g. for Picard)
# - see where to introduce BQSR
# - catch exception if step dir already exists, and skip step for the given lane(s)
# - add step to launch PLINK --mendel after using GATK VariantsToBinaryPed
# - check version of all external programs
# - check that dates in samples file agree with SAM specification
# - add option to ignore R2 files
# - add option to give VCF of known indels (for local realign)
# - try drmaa-python? https://github.com/pygridtools/drmaa-python
# - try sgeparse? https://pypi.python.org/pypi/sgeparse
# - try logging? https://docs.python.org/2/library/logging.html

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
import pip
import pkg_resources
# import xml.dom.minidom
# import sqlite3

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n" % msg)
        sys.exit(1)

## check dependencies
installed_packages = pip.get_installed_distributions()
versions = {package.key: package.version for package in installed_packages}
deps = {"biopython": "1.64",
        "pyutilstimflutre": "0.7"}
for depN,depV in deps.items():
    msg = "ERROR: %s depends on %s in version %s or higher" % \
          (os.path.basename(sys.argv[0]), depN, depV)
    if depN not in versions:
        sys.stderr.write("%s\n" % msg)
        sys.exit(1)
    if pkg_resources.parse_version(versions[depN]) \
       < pkg_resources.parse_version(depV):
        sys.stderr.write("%s\n" % msg)
        sys.exit(1)
        
## http://biopython.org/
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Data.IUPACData import ambiguous_dna_values

## https://github.com/timflutre/pyutilstimflutre
from pyutilstimflutre import Utils, ProgVersion, Job, JobGroup, JobManager, \
    Fastqc, SamtoolsFlagstat
        
progVersion = "0.11.5" # http://semver.org/


class GbsSample(object):
    """
    A GbsSample corresponds to a unique quadruplet (genotype,flowcell,lane,barcode),
    given that the 'genotype' is the focus of the analysis.
    Caution, in some rare cases, the same genotype can be present with multiple barcodes in the same lane: see the method below to get the sample identifier, with (before demultiplexing) or without barcode (after demultiplexing).
    """
    
    def __init__(self, geno, flowcell, laneNum, barcode, demult="before"):
        self.genotype = geno
        self.flowcell = flowcell
        self.lane = laneNum
        self.barcode = barcode
        self.id = "%s_%s_%s" % (self.genotype, self.flowcell, self.lane)
        if demult == "before":
            self.id += "_%s" % self.barcode
        self.refGenome = ""
        self.library = ""
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
        txt += ";genotype=%s" % self.genotype
        txt += ";flowcell=%s" % self.flowcell
        txt += ";lane=%s" % self.lane
        txt += ";refGenome=%s" % self.refGenome
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
        lFilesR1 = glob.glob("%s/*_%s_R1.fastq.gz" % (pathToDir,
                                                      self.genotype))
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
        lFilesR2 = glob.glob("%s/*_%s_R2.fastq.gz" % (pathToDir,
                                                       self.genotype))
        if len(lFilesR2) > 1:
            msg = "%i demultiplexed R2 files found for sample '%s' in '%s'" % \
                  (len(lFilesR2), self.id, pathToDir)
            raise ValueError(msg)
        if len(lFilesR2) == 1:
            self.dDemultiplexedFastqFiles["R2"] = lFilesR2[0]
            
    def clean(self, adpR1, adpR2, errTol, minOvl, minReadLen, minQual,
              maxNPerc, outDir, iJobGroup):
        """
        https://cutadapt.readthedocs.org/en/stable/
        """
        cmd = "# generated by gbs.py %s" % progVersion
        cmd += "\ntime cutadapt"
        cmd += " -a %s" % reverse_complement(str(adpR2)) # to be removed from R1 reads
        if "R2" in self.dDemultiplexedFastqFiles:
            cmd += " -A %s" % reverse_complement(str(adpR1) + str(self.barcode)) # idem from R2
        cmd += " -o %s/%s_clean_R1.fastq.gz" % (outDir, self.id)
        if "R2" in self.dDemultiplexedFastqFiles:
            cmd += " -p %s/%s_clean_R2.fastq.gz" % (outDir, self.id)
        cmd += " -e %f" % errTol # error tolerance
        cmd += " -O %i" % minOvl # min overlap len btw reads and seq passed to -a/-A
        cmd += " -m %i" % minReadLen # min read length
        # cmd += " -U 3" # fixed nb of bases removed from starts of R2 reads
        cmd += " -q %i,%i" % (minQual, minQual) # quality trimming
        cmd += " --max-n %f" % maxNPerc
        # cmd += " --maximum-length 150"
        cmd += " %s" % self.dDemultiplexedFastqFiles["R1"]
        if "R2" in self.dDemultiplexedFastqFiles:
            cmd += " %s" % self.dDemultiplexedFastqFiles["R2"]
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
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
            if Ri in self.dCleanedFastqFiles:
                cmd = "# generated by gbs.py %s" % progVersion
                cmd += "\ntime fastqc -o %s %s" % (outDir,
                                                  self.dCleanedFastqFiles[Ri])
                jobName = "stdout_%s_%s_%s" % (iJobGroup.id, self.id, Ri)
                bashFile = "%s/job_%s_%s_%s.bash" % (outDir, iJobGroup.id,
                                                     self.id, Ri)
                iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                           bashFile=bashFile, dir=outDir)
                iJobGroup.insert(iJob)
                
    def saveNbReadsFromFastqc(self, inDir, outHandle):
        inFile = "%s/%s_clean_R1_fastqc.zip" % (inDir, self.id)
        iFqc = Fastqc(inFile)
        txt = "%s" % self.genotype
        txt += "\t%s" % self.flowcell
        txt += "\t%s" % self.lane
        txt += "\tR1"
        txt += "\t%s" % iFqc.lStats[1]["content"][3]["value"]
        outHandle.write("%s\n" % txt)
        if"R2" in self.dCleanedFastqFiles:
            inFile = "%s/%s_clean_R2_fastqc.zip" % (inDir, self.id)
            iFqc = Fastqc(inFile)
            txt = "%s" % self.genotype
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
        cmd = "# generated by gbs.py %s" % progVersion
        cmd += "\noutDir=\"%s\"" % outDir
        
        cmd += "\n\necho \"align, fixmate and sort...\""
        cmd += "\nbwa mem"
        cmd += " -R \'@RG"
        cmd += "\\tID:%s" % self.id
        cmd += "\\tCN:%s" % self.seqCenter
        cmd += "\\tDT:%s" % self.date
        cmd += "\\tLB:%s" % self.library
        cmd += "\\tPL:%s" % self.seqPlatform
        cmd += "\\tPM:%s" % self.seqPlatformModel
        cmd += "\\tPU:%s-%s.%s" % (self.flowcell, self.barcode, self.lane)
        cmd += "\\tSM:%s" % self.genotype # see GATK's FAQ for what a sample is
        cmd += "\'"
        cmd += " -M %s" % pathToPrefixRefGenome
        cmd += " %s" % self.dCleanedFastqFiles["R1"]
        if "R2" in self.dCleanedFastqFiles:
            cmd += " %s" % self.dCleanedFastqFiles["R2"]
            
        tmpBamFile = "tmp_%s.bam" % self.id
        cmd += " | samtools fixmate"
        cmd += " -O bam"
        cmd += " -" # stdin
        cmd += " -" # stdout
        cmd += " | samtools sort"
        cmd += " -o ${outDir}/%s" % tmpBamFile
        cmd += " -O bam"
        cmd += " -T %s/tmp%s_%s" % (tmpDir, Utils.uniq_alphanum(5), self.id)
        cmd += " -" # stdin
        
        # update the header with @SQ from dictFile
        tmpHeaderFile = "tmp_%s_header.sam" % self.id
        cmd += "\n\necho \"update header...\""
        cmd += "\ncat"
        cmd += " <(samtools view -H ${outDir}/tmp_%s.bam" % self.id
        cmd += " | grep -v '@SQ')"
        cmd += " <(grep '@SQ' %s)" % dictFile
        cmd += " > ${outDir}/%s" % tmpHeaderFile
        cmd += "\ntime samtools reheader"
        cmd += " ${outDir}/%s" % tmpHeaderFile
        cmd += " ${outDir}/%s" % tmpBamFile
        cmd += " > ${outDir}/%s.bam" % self.id
        cmd += "\nrm ${outDir}/%s ${outDir}/%s" % (tmpBamFile, tmpHeaderFile)
        
        # index
        cmd += "\n\necho \"index\""
        cmd += "\nsamtools index"
        cmd += " ${outDir}/%s.bam" % self.id
        
        # basic stats
        cmd += "\n\necho \"flagstat\""
        cmd += "\ntime samtools flagstat"
        cmd += " ${outDir}/%s.bam" % self.id
        cmd += " >& ${outDir}/flagstat_%s.txt" % self.id
        if "R2" in self.dCleanedFastqFiles:
            cmd += "\necho \"paired-end reads in proper pairs and primary alignments\""
        else:
            cmd += "\necho \"single reads in primary alignments\""
        cmd += "\necho -e \"%s" % self.genotype
        cmd += "\t%s" % self.flowcell
        cmd += "\t%s" % self.lane
        if "R2" in self.dCleanedFastqFiles:
            cmd += "\t\"$(samtools view -f 0x0002 -F 0x0100 -q 5"
            cmd += " ${outDir}/%s.bam" % self.id
            cmd += " | cut -f1 | sort | uniq | wc -l)"
        else:
            cmd += "\t\"$(samtools view -F 0x0004 -F 0x0100 -F 0x0800 -q 5"
            cmd += " ${outDir}/%s.bam" % self.id
            cmd += " | cut -f1 | sort | uniq | wc -l)"
        cmd += " > ${outDir}/reads_count-ok_%s.txt" % self.id
        
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
        iJobGroup.insert(iJob)
        
    def setInitialBamFile(self, dirName):
        self.initialBamFile = "%s/%s/%s.bam" % (self.dir, dirName, self.id)
        
    def saveSamtoolsFlagstat(self, outHandle):
        inDir = os.path.dirname(self.initialBamFile)
        inFile = "%s/flagstat_%s.txt" % (inDir, self.id)
        iSf = SamtoolsFlagstat(inFile)
        txt = "%s" % self.genotype
        txt += "\t%s" % self.flowcell
        txt += "\t%s" % self.lane
        txt += "\t%s" % iSf.getTxtToWrite()
        outHandle.write("%s\n" % txt)
        
    def localRealign(self, jvmXms, jvmXmx, pathToPrefixRefGenome, knownIndelsFile,
                     outDir, iJobGroup):
        cmd = "# generated by gbs.py %s" % progVersion
        cmd += "\njava"
        cmd += " -Xms%s" % jvmXms
        cmd += " -Xmx%s" % jvmXmx
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T RealignerTargetCreator"
        cmd += " -R %s.fa" % pathToPrefixRefGenome
        cmd += " -I %s" % self.initialBamFile
        if knownIndelsFile:
            cmd += " --known %s" % knownIndelsFile
        cmd += " -o %s/%s.intervals" % (outDir, self.id)
        cmd += "\njava"
        cmd += " -Xms%s" % jvmXms
        cmd += " -Xmx%s" % jvmXmx
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T IndelRealigner"
        cmd += " -R %s.fa" % pathToPrefixRefGenome
        cmd += " -I %s" % self.initialBamFile
        if knownIndelsFile:
            cmd += " --known %s" % knownIndelsFile
        cmd += " -targetIntervals %s/%s.intervals" % (outDir, self.id)
        cmd += " -o %s/%s_realn.bam" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
        iJobGroup.insert(iJob)
        
    def baseQualityRecalibrate(self, jvmXms, jvmXmx, pathToPrefixRefGenome,
                               knownFile, outDir, iJobGroup):
        cmd = "# generated by gbs.py %s" % progVersion
        cmd += "\njava"
        cmd += " -Xms%s" % jvmXms
        cmd += " -Xmx%s" % jvmXmx
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T BaseRecalibrator"
        cmd += " -R %s.fa" % pathToPrefixRefGenome
        cmd += " -I %s/%s_realn.bam" % (outDir, self.id)
        if knownFile != "":
            cmd += " --known %s" % knownFile
        else:
            cmd += " --run_without_dbsnp_potentially_ruining_quality"
        cmd += " -o %s/%s_recal.table" % (outDir, self.id)
        cmd += "\njava"
        cmd += " -Xms%s" % jvmXms
        cmd += " -Xmx%s" % jvmXmx
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T PrintReads"
        cmd += " -R %s.fa" % pathToPrefixRefGenome
        cmd += " -I %s/%s_realn.bam" % (outDir, self.id)
        cmd += " --BQSR %s/%s_recal.table" % (outDir, self.id)
        cmd += " -o %s/%s_recal.bam" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
        iJobGroup.insert(iJob)
        
        
class GbsLane(object):
    
    def __init__(self, laneId, flowcell, number):
        self.id = laneId # flowcell identifier + "_" + lane number
        self.flowcell = flowcell
        self.number = number
        self.dir = None
        self.dSamples = {} # keys: sample id ; values: GbsSample object
        self.dInitFastqFiles = {} # key(s): R1 (and R2, optional)
                                  # values: [init, symlink]
        
    def insert(self, iSample):
        if iSample.id not in self.dSamples:
            self.dSamples[iSample.id] = iSample
            
    def nbGenos(self):
        genos = set()
        for sampleId,iSample in self.dSamples.items():
            genos.add(iSample.genotype)
        return len(genos)
    
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
            cmd = "echo \"commands generated by gbs.py %s\"" % progVersion
            cmd += "\ntime fastqc -o %s %s" % (outDir, lFiles[1])
            jobName = "stdout_%s_%s_%s" % (iJobGroup.id, self.id, Ri)
            bashFile = "%s/job_%s_%s_%s.bash" % (outDir, iJobGroup.id,
                                                 self.id, Ri)
            iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                       bashFile=bashFile, dir=outDir)
            iJobGroup.insert(iJob)
            
    def saveBarcodeFile(self, outDir, format="table"):
        fileName = "%s/barcodes_%s" % (outDir, self.id)
        if format == "fasta":
            fileName += ".fa"
        elif format == "table":
            fileName += ".tsv"
        else:
            msg = "format for demultiplex.py should be 'table' or 'fasta'"
            raise ValueError(msg)
        
        fileHandle = open(fileName, "w")
        
        if format == "table":
            fileHandle.write("id\ttag\n")
            
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sample in lSamples:
            iSample = self.dSamples[sample]
            if format == "fasta":
                fileHandle.write(">%s\n%s\n" % (iSample.genotype,
                                                iSample.barcode))
            else:
                fileHandle.write("%s\t%s\n" % (iSample.genotype,
                                               iSample.barcode))
        fileHandle.close()
        return fileName
    
    def demultiplex(self, outDir, enzyme, method, nbSubstitutionsAllowed,
                    enforceSubst, iJobGroup):
        cmd = "echo \"commands generated by gbs.py %s\"" % progVersion
        cmd += "\ndemultiplex.py"
        cmd += " --idir %s" % os.path.dirname(self.dInitFastqFiles["R1"][1])
        cmd += " --ifq1 %s" % os.path.basename(self.dInitFastqFiles["R1"][1])
        if "R2" in self.dInitFastqFiles:
            cmd += " --ifq2 %s" % os.path.basename(self.dInitFastqFiles["R2"][1])
        cmd += " --it %s" % self.saveBarcodeFile(outDir, format="table")
        cmd += " --ofqp %s/%s" % (outDir, self.id)
        cmd += " --met %s" % method
        cmd += " --subst %i" % nbSubstitutionsAllowed
        cmd += " --ensubst %s" % enforceSubst
        cmd += " --dist %i" % 0
        cmd += " --re %s" % enzyme
        cmd += " --chim 1"
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
        iJobGroup.insert(iJob)
        
    def setDemultiplexedFastqFiles(self, pathToDir):
        for sampleId,iSample in self.dSamples.items():
            iSample.setDemultiplexedFastqFiles(pathToDir)
                
    def clean(self, adpR1, adpR2, errTol, minOvl, minReadLen, minQual,
              maxNPerc, outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.clean(adpR1, adpR2, errTol, minOvl, minReadLen, minQual,
                          maxNPerc, outDir, iJobGroup)
            
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
            
    def align(self, pathToPrefixRefGenome, tmpDir, dictFile, outDirName,
              iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            outDir = "%s/%s" % (iSample.dir, outDirName)
            if os.path.isdir(outDir):
                shutil.rmtree(outDir)
            os.mkdir(outDir)
            iSample.align(pathToPrefixRefGenome, tmpDir, dictFile, outDir,
                          iJobGroup)
            
    def setInitialBamFiles(self, dirName):
        for sampleId,iSample in self.dSamples.items():
            iSample.setInitialBamFile(dirName)
            
    def gather(self, jvmXms, jvmXmx, tmpDir, outDir, iJobGroup):
        lSamples = self.dSamples.keys()
        cmd = "echo \"commands generated by gbs.py %s\"" % progVersion
        cmd += "\necho \"merge, sort and index BAMs from %s (%i samples)...\"" \
               % (self.id, len(lSamples))
        cmd += "\ntime samtools merge -f"
        cmd += " %s/%s_unsorted.bam" % (outDir, self.id)
        lSamples.sort()
        for sampleId in lSamples:
            cmd += " %s" % self.dSamples[sampleId].initialBamFile
        # https://www.biostars.org/p/251721/
        cmd += "\nsamtools sort"
        cmd += " -o %s/%s.bam" % (outDir, self.id)
        cmd += " -T %s/tmp%s_%s" % (tmpDir, Utils.uniq_alphanum(5), self.id)
        cmd += " %s/%s_unsorted.bam" % (outDir, self.id)
        cmd += "\nrm -f %s/%s_unsorted.bam*" % (outDir, self.id)
        cmd += "\nsamtools index"
        cmd += " %s/%s.bam" % (outDir, self.id)
        cmd += "\n\necho \"collect insert sizes...\""
        cmd += "\njava"
        cmd += " -Xms%s" % jvmXms
        cmd += " -Xmx%s" % jvmXmx
        cmd += " -jar `which picard.jar`"
        cmd += " CollectInsertSizeMetrics"
        cmd += " HISTOGRAM_FILE=%s/hist_insert-sizes_picard_%s.pdf" \
               % (outDir, self.id)
        cmd += " INPUT=%s/%s.bam" % (outDir, self.id)
        cmd += " OUTPUT=%s/insert-sizes_picard_%s.txt" % (outDir, self.id)
        cmd += " VALIDATION_STRINGENCY=%s" % "STRICT" #LENIENT"
        cmd += "\n\necho \"clean...\""
        cmd += "\nrm -f %s/%s.bam*" % (outDir, self.id) # to save space
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
        iJobGroup.insert(iJob)
        
    def saveSamtoolsFlagstat(self, outHandle):
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            iSample = self.dSamples[sampleId]
            iSample.saveSamtoolsFlagstat(outHandle)
            
    def localRealignSamples(self, jvmXms, jvmXmx, pathToPrefixRefGenome,
                            knownIndelsFile, outDirName,
                            iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            outDir = "%s/%s" % (iSample.dir, outDirName)
            if os.path.isdir(outDir):
                shutil.rmtree(outDir)
            os.mkdir(outDir)
            iSample.localRealign(jvmXms, jvmXmx, pathToPrefixRefGenome,
                                 knownIndelsFile, outDir, iJobGroup)
            
    def baseQualityRecalibrate(self, jvmXms, jvmXmx, pathToPrefixRefGenome,
                               knownFile, outDir, iJobGroup):
        for sampleId,iSample in self.dSamples.items():
            iSample.baseQualityRecalibrate(jvmXms, jvmXmx,
                                           pathToPrefixRefGenome, knownFile,
                                           outDir, iJobGroup)
            
            
class GbsGeno(object):
    
    def __init__(self, geno):
        self.id = geno
        self.dir = None
        self.dSamples = {}
        self.lRealignedSampleBamFiles = []
        self.realignedGenoBamFile = ""
        
    def insert(self, iSample):
        if iSample.id not in self.dSamples:
            self.dSamples[iSample.id] = iSample
            
    def setRealignedSampleBamFiles(self, stepDir):
        lSamples = self.dSamples.keys()
        lSamples.sort()
        for sampleId in lSamples:
            iSample = self.dSamples[sampleId]
            self.lRealignedSampleBamFiles.append("%s/%s/%s_realn.bam" \
                                                 % (iSample.dir,
                                                    stepDir,
                                                    iSample.id))
            
    def localRealign(self, jvmXms, jvmXmx, pathToPrefixRefGenome,
                     knownIndelsFile, outDir, iJobGroup):
        cmd = "echo \"commands generated by gbs.py %s\"" % progVersion
        if len(self.lRealignedSampleBamFiles) == 1:
            pathWoExt = os.path.splitext(self.lRealignedSampleBamFiles[0])[0]
            cmd += "\nln -s %s.bam" % pathWoExt
            cmd += " %s/%s_realn.bam" % (outDir, self.id)
            cmd += "\nln -s %s.bai" % pathWoExt
            cmd += " %s/%s_realn.bai" % (outDir, self.id)
        else:
            cmd += "\njava"
            cmd += " -Xms%s" % jvmXms
            cmd += " -Xmx%s" % jvmXmx
            cmd += " -jar `which GenomeAnalysisTK.jar`"
            cmd += " -T RealignerTargetCreator"
            cmd += " -R %s.fa" % pathToPrefixRefGenome
            for i in range(len(self.lRealignedSampleBamFiles)):
                cmd += " -I %s" % self.lRealignedSampleBamFiles[i]
            if knownIndelsFile:
                cmd += " --known %s" % knownIndelsFile
            cmd += " -o %s/%s.intervals" % (outDir, self.id)
            cmd += "\njava"
            cmd += " -Xms%s" % jvmXms
            cmd += " -Xmx%s" % jvmXmx
            cmd += " -jar `which GenomeAnalysisTK.jar`"
            cmd += " -T IndelRealigner"
            cmd += " -R %s.fa" % pathToPrefixRefGenome
            for i in range(len(self.lRealignedSampleBamFiles)):
                cmd += " -I %s" % self.lRealignedSampleBamFiles[i]
            if knownIndelsFile:
                cmd += " --known %s" % knownIndelsFile
            cmd += " -targetIntervals %s/%s.intervals" % (outDir, self.id)
            cmd += " -o %s/%s_realn.bam" % (outDir, self.id)
        jobName = "stdout_%s_%s" % (iJobGroup.id, self.id)
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
        iJobGroup.insert(iJob)
        
    def setRealignedGenotypeBamFile(self, pathToDir):
        self.realignedGenoBamFile = "%s/%s_realn.bam" % (pathToDir, self.id)
        
    def variantCalling(self, jvmXms, jvmXmx, tmpDir, pathToPrefixRefGenome,
                       knownFile, outDir, iJobGroup, saveActiveRegions=False,
                       saveActivityProfiles=False):
        """
        https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
        http://gatkforums.broadinstitute.org/discussion/comment/14337/#Comment_14337
        """
        cmd = "echo \"commands generated by gbs.py %s\"" % progVersion
        cmd += "\njava"
        cmd += " -Xms%s" % jvmXms
        cmd += " -Xmx%s" % jvmXmx
        if tmpDir:
            cmd += " -Djava.io.tmpdir=%s" % tmpDir
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T HaplotypeCaller"
        cmd += " -R %s.fa" % pathToPrefixRefGenome
        cmd += " -I %s" % self.realignedGenoBamFile
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
        if saveActiveRegions:
            cmd += " --activeRegionOut %s/%s_active-regions.tab" % (outDir, self.id)
        if saveActivityProfiles:
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
        bashFile = "%s/job_%s_%s.bash" % (outDir, iJobGroup.id, self.id)
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=outDir)
        iJobGroup.insert(iJob)
        
        
class Gbs(object):
    
    def __init__(self):
        self.verbose = 1
        self.project1Id = None
        self.project2Id = None
        self.scheduler = "SGE"
        self.queue = "normal.q"
        self.lResources = None
        self.rmvBash = False
        self.lSteps = []
        self.forceRerunSteps = False
        self.samplesFile = None
        self.fclnToKeep = None
        self.pathToInReadsDir = ""
        self.enzyme = "ApeKI"
        self.dmxMethod = "4c"
        self.nbSubstsAllowedDemult = 2
        self.enforceSubst = "lenient"
        self.jobManager = None # instantiated in run()
        self.samplesCol2idx = {"genotype": None,
                               "flowcell": None,
                               "lane": None,
                               "ref_genome": None,
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
        self.dFlowcells = {}
        self.dGenos = {}
        self.allLanesDir = None # depends on cwd and project1Id
        self.allSamplesDir = None # depends on cwd and project2Id
        self.allGenosDir = None # depends on cwd and project2Id
        self.lDirSteps = ["init-quality",
                          "demultiplex",
                          "clean-reads",
                          "align-reads",
                          "realign-reads-sample",
                          "realign-reads-geno",
                          "call-variants-geno",
                          "joint-genotyping", # to be completed by jointGenoId
                          "variant-geno-filter"]
        self.adpFile = None
        self.adapters = {}
        self.errTol = 0.2
        self.minOvl = 3
        self.minReadLen = 35
        self.minQual = 20
        self.maxNPerc = 0.2
        self.pathToPrefixRefGenome = None
        self.dictFile = None
        self.jointGenoId = None
        self.restrictAllelesTo = "ALL"
        self.minDp = None
        self.minGq = None
        self.maxNbFilterGenos = None
        self.maxFracFilterGenos = None
        self.maxNbNocallGenos = None
        self.maxFracNocallGenos = None
        self.famFile = None
        self.mendelianViolationQualThreshold = 0
        self.excludeSampleFile = None
        self.tmpDir = "."
        self.jvmXms = "512m"
        self.jvmXmx = "4g"
        self.queue2 = "bigmem.q"
        self.knownIndelsFile = None
        self.knownFile = None
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' performs the computational aspects of genotyping-by-sequencing.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "      --proj1\tname of the project used for steps 1 to 4\n"
        msg += "\t\tmention a reference genome only if all samples belong to\n"
        msg += "\t\t the same species, and will be mapped to the same ref genome\n"
        msg += "      --proj2\tname of the project used for steps 4 to 8\n"
        msg += "\t\tcan be the same as --proj1, or can be different\n"
        msg +="\t\t notably when samples come from different species\n"
        msg += "\t\t or if one wants to align reads to different ref genomes\n"
        msg += "      --schdlr\tname of the cluster scheduler (default=SGE)\n"
        msg += "      --queue\tname of the cluster queue (default=normal.q)\n"
        msg += "      --resou\tcluster resources (e.g. 'test' for 'qsub -l test')\n"
        msg += "      --rmvb\tremove bash scripts for jobs launched in parallel\n"
        msg += "      --step\tstep to perform (1/2/3/.../9)\n"
        msg += "\t\t1: raw read quality per lane (with FastQC v >= 0.11.2)\n"
        msg += "\t\t2: demultiplexing per lane (with demultiplex.py v >= 1.14.0)\n"
        msg += "\t\t3: cleaning per sample (with CutAdapt v >= 1.8)\n"
        msg += "\t\t4: alignment per sample (with BWA MEM v >= 0.7.12, Samtools v >= 1.3, Picard and R v >= 3)\n"
        msg += "\t\t5: local realignment per sample (with GATK v >= 3.5)\n"
        msg += "\t\t6: local realignment per genotype (with GATK v >= 3.5)\n"
        msg += "\t\t7: variant and genotype calling per genotype (with GATK HaplotypeCaller v >= 3.5)\n"
        msg += "\t\t8: variant and genotype calling jointly across genotypes (with GATK GenotypeGVCFs v >= 3.5)\n"
        msg += "\t\t9: variant and genotype filtering (with GATK v >= 3.5)\n"
        msg += "      --samples\tpath to the 'samples' file\n"
        msg += "\t\tcompulsory for all steps, but can differ between steps\n"
        msg += "\t\t e.g. if samples come from different species or are aligned\n"
        msg += "\t\t on different ref genomes, different samples file should\n"
        msg += "\t\t be used for steps 4-9, representing different subsets of\n"
        msg += "\t\t the file used for steps 1-3\n"
        msg += "\t\tthe file should be encoded in ASCII\n"
        msg += "\t\tthe first row should be a header with column names\n"
        msg += "\t\teach 'sample' (see details below) should have one and only one row\n"
        msg += "\t\tany two columns should be separated with one tabulation\n"
        msg += "\t\tcolumns can be in any order\n"
        msg += "\t\trows starting by '#' are skipped\n"
        msg += "\t\t12 columns are compulsory (but there can be more):\n"
        msg += "\t\t genotype (see details below, e.g. 'Col-0', but use neither underscore '_' nor space ' ' nor dot '.', use dash '-' instead)\n"
        msg += "\t\t ref_genome (identifier of the reference genome used for alignment, e.g. 'Atha_v2', but use neither space ' ' nor dot '.'; the full species name, e.g. 'Arabidopsis thaliana', will be present in the file given to --dict)\n"
        msg += "\t\t library (e.g. can be the same as 'genotype')\n"
        msg += "\t\t barcode (e.g. 'ATGG')\n"
        msg += "\t\t seq_center (e.g. 'Broad Institute', 'GenoToul', etc)\n"
        msg += "\t\t seq_platform (e.g. 'ILLUMINA', see SAM format specification)\n"
        msg += "\t\t seq_platform_model (e.g. 'HiSeq 2000')\n"
        msg += "\t\t flowcell (e.g. 'C5YMDACXX')\n"
        msg += "\t\t lane (e.g. '3', can be '31' if a first demultiplexing was done per index)\n"
        msg += "\t\t date (e.g. '2015-01-15', see SAM format specification)\n"
        msg += "\t\t fastq_file_R1 (filename, one per lane, gzip-compressed)\n"
        msg += "\t\t fastq_file_R2 (filename, one per lane, gzip-compressed)\n"
        msg += "      --fcln\tidentifier of a flowcell and lane number\n"
        msg += "\t\tformat as <flowcell>_<lane-number>, e.g. 'C5YMDACXX_1'\n"
        msg += "\t\tif set, only the samples from this lane will be analyzed\n"
        msg += "      --pird\tpath to the input reads directory\n"
        msg += "\t\tcompulsory for steps 1 and 2\n"
        msg += "\t\twill be added to the columns 'fastq_file_R*' from the sample file\n"
        msg += "\t\tif not set, input read files should be in current directory\n"
        msg += "      --enz\tname of the restriction enzyme\n"
        msg += "\t\tcompulsory for step 2\n"
        msg += "\t\tdefault=ApeKI\n"
        msg += "      --dmxmet\tmethod used to demultiplex\n"
        msg += "\t\tcompulsory for step 2\n"
        msg += "\t\tdefault=4c (see the help of demultiplex.py to know more)\n"
        msg += "      --subst\tnumber of substitutions allowed during demultiplexing\n"
        msg += "\t\tcompulsory for step 2\n"
        msg += "\t\tdefault=2\n"
        msg += "      --ensubst\tenforce the nb of substitutions allowed\n"
        msg += "\t\tcompulsory for step 2\n"
        msg += "\t\tdefault=lenient/strict\n"
        msg += "      --adp\tpath to the file containing the adapters\n"
        msg += "\t\tcompulsory for step 3\n"
        msg += "\t\tsame format as FastQC: name<tab>sequence\n"
        msg += "\t\tname: at least 'adpR1' (also 'adpR2' if paired-end)\n"
        msg += "\t\tsequence: from 5' (left) to 3' (right)\n"
        msg += "       --errtol\terror tolerance to find adapters\n"
        msg += "\t\tcompulsory for step 3\n"
        msg += "\t\tdefault=0.2\n"
        msg += "       --minovl\tminimum overlap length between reads and adapters\n"
        msg += "\t\tcompulsory for step 3\n"
        msg += "\t\tdefault=3 (in bases)\n"
        msg += "       --minrl\tminimum length to keep a read\n"
        msg += "\t\tcompulsory for step 3\n"
        msg += "\t\tdefault=35 (in bases)\n"
        msg += "       --minq\tminimum quality to trim a read\n"
        msg += "\t\tcompulsory for step 3\n"
        msg += "\t\tdefault=20 (used for both reads if paired-end)\n"
        msg += "       --maxNp\tmaximum percentage of N to keep a read\n"
        msg += "\t\tcompulsory for step 3\n"
        msg += "\t\tdefault=0.2\n"
        msg += "      --ref\tpath to the prefix of files for the reference genome\n"
        msg += "\t\tcompulsory for steps 4, 5, 6, 7, 8, 9\n"
        msg += "\t\tshould correspond to the 'ref_genome' column in --samples\n"
        msg += "\t\te.g. '/data/Atha_v2' for '/data/Atha_v2.fa', '/data/Atha_v2.bwt', etc\n"
        msg += "\t\tthese files are produced via 'bwa index ...'\n"
        msg += "      --dict\tpath to the 'dict' file (SAM header with @SQ tags)\n"
        msg += "\t\tcompulsory for step 4\n"
        msg += "\t\tsee 'CreateSequenceDictionary' in the Picard software\n"
        msg += "      --jgid\tcohort identifier to use for joint genotyping\n"
        msg += "\t\tcompulsory for steps 8, 9\n"
        msg += "\t\tuseful to launch several, different cohorts in parallel\n"
        msg += "      --rat\trestrict alleles to be of a particular allelicity\n"
        msg += "\t\tused in step 9\n"
        msg += "\t\tdefault=ALL/BIALLELIC/MULTIALLELIC\n"
        msg += "\t\tsee '--restrictAllelesTo' in GATK's SelectVariant\n"
        msg += "      --mdp\tminimum value for DP (read depth; e.g. 10)\n"
        msg += "\t\tused in step 9\n"
        msg += "\t\tsee GATK's VariantFiltration\n"
        msg += "      --mgq\tminimum value for GQ (genotype quality; e.g. 20)\n"
        msg += "\t\tused in step 9\n"
        msg += "\t\tsee GATK's VariantFiltration\n"
        msg += "      --mnfg\tmaximum number of filtered genotypes to keep a variant\n"
        msg += "\t\tused in step 9\n"
        msg += "\t\tsee '--maxFilteredGenotypes' in GATK's SelectVariants\n"
        msg += "      --mffg\tmaximum fraction of filtered genotypes to keep a variant\n"
        msg += "\t\tused in step 9\n"
        msg += "\t\tsee '--maxFractionFilteredGenotypes' in GATK's SelectVariants\n"
        msg += "      --mnnc\tmaximum number of not-called genotypes to keep a variant\n"
        msg += "\t\tused in step 9\n"
        msg += "      --mfnc\tmaximum fraction of not-called genotypes to keep a variant\n"
        msg += "\t\tused in step 9\n"
        msg += "\t\tsee '--maxNOCALLfraction' in GATK's SelectVariants\n"
        msg += "      --fam\tpath to the file containing pedigree information\n"
        msg += "\t\tused in step 9\n"
        msg += "\t\tdiscard variants with Mendelian violations (see Semler et al, 2012)\n"
        msg += "\t\tshould be in the 'fam' format specified by PLINK\n"
        msg += "\t\tvalidation strictness (GATK '-pedValidationType') is set at 'SILENT'\n"
        msg += "\t\t allowing some samples to be absent from the pedigree\n"
        msg += "      --mvq\tminimum GQ for each trio member to accept a variant as a Mendelian violation\n"
        msg += "\t\tused in step 9 if '--fam' is specified\n"
        msg += "\t\tdefault=0\n"
        msg += "      --xlssf\tpath to the file with genotypes to exclude\n"
        msg += "\t\tused in step 9 (can be especially useful if '--fam' is specified)\n"
        msg += "      --tmpd\tpath to a temporary directory on child nodes (default=.)\n"
        msg += "\t\te.g. it can be /tmp or /scratch\n"
        msg += "\t\tused in step 4 for 'samtools sort'\n"
        msg += "\t\tused in step 7 for 'GATK HaplotypeCaller'\n"
        msg += "      --jvmXms\tinitial memory allocated to the Java Virtual Machine\n"
        msg += "\t\tdefault=512m (can also be specified as 1024k, 1g, etc)\n"
        msg += "\t\tused in steps 4, 5, 6, 7 and 8 for Picard and GATK\n"
        msg += "      --jvmXmx\tmaximum memory allocated to the Java Virtual Machine\n"
        msg += "\t\tdefault=4g\n"
        msg += "\t\tused in steps 4, 5, 6, 7 and 8 for Picard and GATK\n"
        msg += "      --queue2\tname of the second cluster queue (default=bigmem.q)\n"
        msg += "\t\tused in step 4 for Picard to collect insert sizes\n"
        msg += "      --knowni\tpath to a VCF file with known indels (for local realignment)\n"
        msg += "      --known\tpath to a VCF file with known variants (e.g. from dbSNP)\n"
        msg += "      --force\tforce to re-run step(s)\n"
        msg += "\t\tthis removes without warning the step directory if it exists\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --step 1 --samples samples.txt\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Details:\n"
        msg += "This program aims at genotyping a set of 'genotypes' using data from\n"
        msg += "a restriction-assisted DNA sequencing (RAD-seq) experiment, also known\n"
        msg += "as a genotyping-by-sequencing (GBS) experiment.\n"
        msg += "Here, by 'genotype', we mean the entity which is the focus of the\n"
        msg += "study. For instance, it can be a plant variety (or a human being), or\n"
        msg += "the specific clone of a given plant variety (or a specific tumor of a\n"
        msg += "given human being), etc.\n"
        msg += "Importantly, note that the content of the 'genotype' column will\n"
        msg += "be used to set the 'SM' (sample) tag of the 'RG' (read group) header\n"
        msg += "record type of the SAM format (see http://www.htslib.org/). However,\n"
        msg += "internal to this program, the term 'sample' corresponds to the unique\n"
        msg += "quadruplet (genotype,flowcell,lane,barcode) for steps 1 and 2, and to\n"
        msg += "the unique triplet (genotype,flowcell,lane) for the others.\n"
        msg += "Jobs are executed in parallel (--schdlr). Their return status is\n"
        msg += "recorded in a SQLite database which is removed at the end. If a job\n"
        msg += "fails, the whole script stops with an error.\n"
        msg += "\n"
        msg += "Dependencies:\n"
        msg += "Python >= 2.7; Biopython; pyutilstimflutre >= 0.5\n"
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        
        The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012).
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2015-2017 Institut National de la Recherche Agronomique.\n"
        msg += "License AGPLv3+: GNU AGPL version 3 or later <http://gnu.org/licenses/agpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt(sys.argv[1:],
                                       "hVv:i:",
                                       ["help", "version", "verbose=",
                                        "proj1=",
                                        "proj2=",
                                        "step=",
                                        "samples=",
                                        "fcln=",
                                        "dict=",
                                        "schdlr=",
                                        "queue=",
                                        "enz=",
                                        "dmxmet=",
                                        "subst=",
                                        "ensubst=",
                                        "adp=",
                                        "errtol=",
                                        "minovl=",
                                        "minrl=",
                                        "minq=",
                                        "maxNp=",
                                        "ref=",
                                        "jgid=",
                                        "rat=",
                                        "mdp=",
                                        "mgq=",
                                        "mnfg=",
                                        "mffg=",
                                        "mnnc=",
                                        "mfnc=",
                                        "fam=",
                                        "mvq=",
                                        "xlssf=",
                                        "tmpd=",
                                        "jvmXms=",
                                        "jvmXmx=",
                                        "queue2=",
                                        "knowni=",
                                        "known=",
                                        "force",
                                        "pird=",
                                        "resou=",
                                        "rmvb"])
        except getopt.GetoptError as err:
            sys.stderr.write("%s\n\n" % str(err))
            # self.help()
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
            elif o == "--proj1":
                self.project1Id = a
            elif o == "--proj2":
                self.project2Id = a
            elif o == "--schdlr":
                self.scheduler = a
            elif o == "--queue":
                self.queue = a
            elif o == "--resou":
                self.lResources = a.split()
            elif o == "--rmvb":
                self.rmvBash = True
            elif o == "--step":
                self.lSteps = [a]
            elif o == "--samples":
                 self.samplesFile = a
            elif o == "--fcln":
                self.fclnToKeep = a
            elif o == "--pird":
                self.pathToInReadsDir = a
            elif o == "--enz":
                self.enzyme = a
            elif o == "--dmxmet":
                self.dmxMethod = a
            elif o == "--subst":
                self.nbSubstsAllowedDemult = int(a)
            elif o == "--ensubst":
                self.enforceSubst = a
            elif o == "--adp":
                self.adpFile = a
            elif o == "--errtol":
                self.errTol = float(a)
            elif o == "--minovl":
                self.minOvl = int(a)
            elif o == "--minrl":
                self.minReadLen = int(a)
            elif o == "--minq":
                self.minQual = int(a)
            elif o == "--maxNp":
                self.maxNPerc = float(a)
            elif o == "--ref":
                self.pathToPrefixRefGenome = a
            elif o == "--dict":
                 self.dictFile = a
            elif o == "--jgid":
                 self.jointGenoId = a
            elif o == "--rat":
                self.restrictAllelesTo = a
            elif o == "--mdp":
                self.minDp = int(a)
            elif o == "--mgq":
                self.minGq = int(a)
            elif o == "--mnfg":
                self.maxNbFilterGenos = int(a)
            elif o == "--mffg":
                self.maxFracFilterGenos = float(a)
            elif o == "--mnnc":
                self.maxNbNocallGenos = int(a)
            elif o == "--mfnc":
                self.maxFracNocallGenos = float(a)
            elif o == "--fam":
                self.famFile = a
            elif o == "--mvq":
                self.mendelianViolationQualThreshold = int(a)
            elif o == "--xlssf":
                self.excludeSampleFile = a
            elif o == "--tmpd":
                self.tmpDir = a
            elif o == "--jvmXms":
                self.jvmXms = a
            elif o == "--jvmXmx":
                self.jvmXmx = a
            elif o == "--queue2":
                self.queue2 = a
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
        if len(self.lSteps) == 0:
            msg = "ERROR: missing compulsory option --step"
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if len(self.lSteps) > 1:
            msg = "ERROR: --step takes a single step"
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if self.lSteps[0] not in ["1","2","3","4","5","6","7","8","9"]:
            msg = "ERROR: unknown --step %s" % self.lSteps[0]
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if "1" in self.lSteps or "2" in self.lSteps or "3" in self.lSteps \
           or "4" in self.lSteps:
            if not self.project1Id:
                msg = "ERROR: missing compulsory option --proj1"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if "4" in self.lSteps or "5" in self.lSteps or "6" in self.lSteps \
           or "7" in self.lSteps or "8" in self.lSteps:
            if not self.project2Id:
                msg = "ERROR: missing compulsory option --proj2"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if self.project1Id and "_" in self.project1Id:
            msg = "ERROR: forbidden underscore '_' in project identifier '%s'" \
                  % self.project1Id
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if self.project2Id and "_" in self.project2Id:
            msg = "ERROR: forbidden underscore '_' in project identifier '%s'" \
                  % self.project2Id
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if not self.samplesFile:
            msg = "ERROR: missing compulsory option --samples"
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if not os.path.exists(self.samplesFile):
            msg = "ERROR: can't find file %s" % self.samplesFile
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if not self.scheduler:
            msg = "ERROR: missing compulsory option --schdlr"
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if self.scheduler == "OGE":
            self.scheduler = "SGE"
        if not self.queue:
            msg = "ERROR: missing compulsory option --queue"
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if self.lSteps == []:
            msg = "ERROR: missing compulsory option --step"
            sys.stderr.write("%s\n\n" % msg)
            # self.help()
            sys.exit(1)
        if "1" in self.lSteps:
            if not Utils.isProgramInPath("fastqc"):
                msg = "ERROR: can't find 'fastqc' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if "2" in self.lSteps:
            if not Utils.isProgramInPath("demultiplex.py"):
                msg = "ERROR: can't find 'demultiplex.py' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            obsMajVer, obsMinVer = ProgVersion.getVersion("demultiplex.py")
            if not (obsMajVer == 1 and obsMinVer >= 14):
                msg = "ERROR: 'demultiplex.py' is in version %s.%s" % \
                      (obsMajVer, obsMinVer)
                msg += " instead of >= 1.14.0"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if "3" in self.lSteps:
            if not Utils.isProgramInPath("cutadapt"):
                msg = "ERROR: can't find 'cutadapt' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not self.adpFile:
                msg = "ERROR: missing compulsory option --adp"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not os.path.exists(self.adpFile):
                msg = "ERROR: can't find file %s" % self.adpFile
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if self.maxNPerc < 0 or self.maxNPerc > 1:
                msg = "ERROR: --maxNp %f should be between 0 and 1" \
                      % self.maxNPerc
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if "4" in self.lSteps:
            if not Utils.isProgramInPath("bwa"):
                msg = "ERROR: can't find 'bwa' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not Utils.isProgramInPath("samtools"):
                msg = "ERROR: can't find 'samtools' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not Utils.isProgramInPath("picard.jar"):
                msg = "ERROR: can't find 'picard.jar' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not self.dictFile:
                msg = "ERROR: missing compulsory option --dict"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not os.path.exists(self.dictFile):
                msg = "ERROR: can't find file %s" % self.dictFile
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if os.path.dirname(self.dictFile) == '':
                self.dictFile = "%s/%s" % (os.getcwd(), self.dictFile)
            if not self.queue2:
                msg = "ERROR: missing compulsory option --queue2"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if "5" in self.lSteps:
            if not Utils.isProgramInPath("GenomeAnalysisTK.jar"):
                msg = "ERROR: can't find 'GenomeAnalysisTK.jar' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            obsMajVer, obsMinVer = ProgVersion.getVersionGatk()
            expMajVer = 3
            expMinVer = 5
            if not (obsMajVer == expMajVer and obsMinVer >= expMinVer):
                msg = "ERROR: 'GATK' is in version %s.%s" % \
                      (obsMajVer, obsMinVer)
                msg += " instead of >= %i.%i" % (expMajVer, expMinVer)
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if self.knownIndelsFile and not os.path.exists(self.knownIndelsFile):
                msg = "ERROR: can't find file %s" % self.knownIndelsFile
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if "6" in self.lSteps or "7" in self.lSteps or "8" in self.lSteps or \
           "9" in self.lSteps:
            if not Utils.isProgramInPath("GenomeAnalysisTK.jar"):
                msg = "ERROR: can't find 'GenomeAnalysisTK.jar' in PATH"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            obsMajVer, obsMinVer = ProgVersion.getVersionGatk()
            if not (obsMajVer == 3 and obsMinVer >= 5):
                msg = "ERROR: 'GATK' is in version %s.%s" % \
                      (obsMajVer, obsMinVer)
                msg += " instead of >= 3.5"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
        if "4" in self.lSteps or "5" in self.lSteps or "6" in self.lSteps or \
           "7" in self.lSteps or "8" in self.lSteps or "9" in self.lSteps:
            if not self.pathToPrefixRefGenome:
                msg = "ERROR: missing compulsory option --ref"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not os.path.exists("%s.bwt" % self.pathToPrefixRefGenome):
                msg = "ERROR: can't find file %s.bwt" % self.pathToPrefixRefGenome
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if not os.path.exists("%s.fa.fai" % self.pathToPrefixRefGenome):
                msg = "ERROR: can't find file %s.fa.fai" % self.pathToPrefixRefGenome
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if os.path.dirname(self.pathToPrefixRefGenome) == "":
                self.pathToPrefixRefGenome = "%s/%s" % (os.getcwd(),
                                                        self.pathToPrefixRefGenome)
        if "8" in self.lSteps or "9" in self.lSteps:
            if not self.jointGenoId:
                msg = "ERROR: missing compulsory option --jgid"
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
                
        if "9" in self.lSteps:
            if self.restrictAllelesTo not in ["ALL", "BIALLELIC",
                                              "MULTIALLELIC"]:
                msg = "ERROR: unknown option --rat %s" % self.restrictAllelesTo
                sys.stderr.write("%s\n\n" % msg)
                # self.help()
                sys.exit(1)
            if self.famFile:
                if not os.path.exists(self.famFile):
                    msg = "ERROR: can't find file %s" % self.famFile
                    sys.stderr.write("%s\n\n" % msg)
                    # self.help()
                    sys.exit(1)
            if self.excludeSampleFile:
                if not os.path.exists(self.excludeSampleFile):
                    msg = "ERROR: can't find file %s" % self.excludeSampleFile
                    sys.stderr.write("%s\n\n" % msg)
                    # self.help()
                    sys.exit(1)
                    
                    
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
        """
        Fill self.dSamples, self.dLanes, self.dFlowcells and self.dGenos.
        """
        refgenome = set()
        
        for line in lines:
            if line.startswith("#"):
                continue
            
            tokens = line.rstrip("\n").split("\t")
            
            # create and fill a "GbsSample" object
            for samplesCol in ["genotype", "flowcell", "lane"]:
                if "_" in tokens[self.samplesCol2idx[samplesCol]]:
                    msg = "underscore in %s '%s', replace by dash '-'" \
                          % (samplesCol, tokens[self.samplesCol2idx[samplesCol]])
                    raise ValueError(msg)
                if " " in tokens[self.samplesCol2idx[samplesCol]]:
                    msg = "space in %s '%s', replace by dash '-'" \
                          % (samplesCol, tokens[self.samplesCol2idx[samplesCol]])
                    raise ValueError(msg)
                if "." in tokens[self.samplesCol2idx[samplesCol]]:
                    msg = "dot in %s '%s', replace by dash '-'" \
                          % (samplesCol, tokens[self.samplesCol2idx[samplesCol]])
                    raise ValueError(msg)
            geno = tokens[self.samplesCol2idx["genotype"]]
            flowcell = tokens[self.samplesCol2idx["flowcell"]]
            laneNum = int(tokens[self.samplesCol2idx["lane"]])
            barcode = tokens[self.samplesCol2idx["barcode"]]
            if self.fclnToKeep is not None and \
               "%s_%i" % (flowcell, laneNum) != self.fclnToKeep:
                continue
            iSample = GbsSample(geno, flowcell, laneNum, barcode,
                                "before" if int(self.lSteps[0]) < 3 \
                                else "after")
            ## at this stage: iSample.id = genotype_flowcell_lane_barcode
            iSample.refGenome = tokens[self.samplesCol2idx["ref_genome"]]
            iSample.library = tokens[self.samplesCol2idx["library"]]
            iSample.seqCenter = tokens[self.samplesCol2idx["seq_center"]]
            iSample.seqPlatform = tokens[self.samplesCol2idx["seq_platform"]]
            iSample.seqPlatformModel = tokens[self.samplesCol2idx["seq_platform_model"]]
            iSample.date = tokens[self.samplesCol2idx["date"]]
            iSample.initFastqFile1 \
                = "%s/%s" % (self.pathToInReadsDir,
                             tokens[self.samplesCol2idx["fastq_file_R1"]])
            if tokens[self.samplesCol2idx["fastq_file_R2"]] != "":
                iSample.initFastqFile2 \
                    = "%s/%s" % (self.pathToInReadsDir,
                                 tokens[self.samplesCol2idx["fastq_file_R2"]])
            self.dSamples[iSample.id] = iSample
            refgenome.add(iSample.refGenome)

            if flowcell not in self.dFlowcells:
                self.dFlowcells[flowcell] = []
            if laneNum not in self.dFlowcells[flowcell]:
                self.dFlowcells[flowcell].append(laneNum)

            # create and fill a "GbsLane" object
            laneId = "%s_%i" % (flowcell, laneNum)
            if laneId not in self.dLanes:
                self.dLanes[laneId] = GbsLane(laneId, flowcell, laneNum)
            self.dLanes[laneId].insert(iSample)

            # create and fill a "GbsGeno" object
            if geno not in self.dGenos:
                self.dGenos[geno] = GbsGeno(geno)
            self.dGenos[geno].insert(iSample)
            
        if ("4" in self.lSteps or "5" in self.lSteps or "6" in self.lSteps \
            or "7" in self.lSteps or "8" in self.lSteps) \
            and len(refgenome) > 1:
            print(refgenome)
            msg = "samples file contains more than one reference genome"
            raise ValueError(msg)
        
        
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
            msg = "nb of genotypes: %i" % len(self.dGenos)
            msg += "\nnb of samples: %i" % len(self.dSamples)
            msg += "\nnb of flowcell%s: %i" \
                   % ("s" if len(self.dFlowcells) > 1 else "",
                      len(self.dFlowcells))
            msg += "\nnb of lane%s: %i" \
                   % ("s" if len(self.dLanes) > 1 else "",
                      len(self.dLanes))
            for flowcell in self.dFlowcells:
                for laneNum in self.dFlowcells[flowcell]:
                    laneId = "%s_%i" % (flowcell, laneNum)
                    msg += "\nflowcell %s, lane %i: %i samples" \
                           % (flowcell, laneNum,
                              len(self.dLanes[laneId].dSamples))
                    if self.lSteps[0] in ["1","2"]:
                        msg += ", %i genotypes" % self.dLanes[laneId].nbGenos()
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def setupLaneDirectories(self):
        """
        One sub-directory per lane.
        """
        self.allLanesDir = "%s/%s_all-lanes" % (os.getcwd(), self.project1Id)
        if not os.path.exists(self.allLanesDir):
            os.mkdir(self.allLanesDir)
        for laneId,iLane in self.dLanes.items():
            dirLane = "%s/%s" % (self.allLanesDir, laneId)
            iLane.dir = dirLane
            if not os.path.exists(dirLane):
                os.mkdir(dirLane)
        if self.verbose > 0:
            msg = "lane directories: %s" % self.allLanesDir
            print(msg); sys.stdout.flush()


    def setupSampleDirectories(self):
        """
        One sub-directory per sample.
        """
        self.allSamplesDir = "%s/%s_all-samples" % (os.getcwd(),
                                                    self.project2Id)
        if not os.path.exists(self.allSamplesDir):
            os.mkdir(self.allSamplesDir)
        for sampleId,iSample in self.dSamples.items():
            dirSample = "%s/%s" % (self.allSamplesDir, sampleId)
            iSample.dir = dirSample
            if not os.path.exists(dirSample):
                os.mkdir(dirSample)
        if self.verbose > 0:
            msg = "sample directories: %s" % self.allSamplesDir
            print(msg); sys.stdout.flush()


    def setupGenotypeDirectories(self):
        """
        One sub-directory per genotype.
        """
        self.allGenosDir = "%s/%s_all-genotypes" % (os.getcwd(),
                                                    self.project2Id)
        if not os.path.exists(self.allGenosDir):
            os.mkdir(self.allGenosDir)
        for genoId,iGeno in self.dGenos.items():
            dirGeno = "%s/%s" % (self.allGenosDir, genoId)
            iGeno.dir = dirGeno
            if not os.path.exists(dirGeno):
                os.mkdir(dirGeno)
        if self.verbose > 0:
            msg = "genotype directories: %s" % self.allGenosDir
            print(msg); sys.stdout.flush()


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
                
                # http://stackoverflow.com/a/24862213/597069
                try:
                    if os.getpgrp() == os.tcgetpgrp(sys.stdout.fileno()):
                        msg = "Do you want to remove it? [y/n] "
                        wantRmvDir = Utils.user_input(msg)
                        if wantRmvDir == "y":
                            shutil.rmtree(dirName)
                        else:
                            raise OSError("can't continue")
                    else: # running in background
                        raise OSError("use --force, or remove it yourself and re-launch me.")
                except OSError:
                    raise OSError("use --force, or remove it yourself and re-launch me.")
                    
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
                     self.lDirSteps[stepNum - 1])
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
            
    def launchFastqcOnInputFastqFiles(self):
        if self.verbose > 0:
            msg = "assess quality per lane..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step1-fastqc" % self.project1Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[0])
            self.makeDir(stepDir)
            iLane.setInitFastqFiles()
            iLane.makeInitFastqFileSymlinks()
            iLane.initQc(stepDir, iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
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
        
        
    def demultiplexInputFastqFiles(self):
        if self.verbose > 0:
            msg = "demultiplex samples per lane..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step2-demultiplex" % self.project1Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[1])
            self.makeDir(stepDir)
            iLane.setInitFastqFiles()
            iLane.makeInitFastqFileSymlinks()
            iLane.demultiplex(stepDir, self.enzyme, self.dmxMethod,
                              self.nbSubstsAllowedDemult, self.enforceSubst,
                              iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step2(self):
        cwd, startTime = self.beginStep(2)
        self.demultiplexInputFastqFiles()
        self.endStep(2, cwd, startTime)
        
        
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
            
            
    def cleanDemultiplexedFiles(self):
        if self.verbose > 0:
            msg = "clean reads per sample..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step3-clean" % self.project1Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[2])
            self.makeDir(stepDir)
            iLane.setDemultiplexedFastqFiles("%s/%s" % (iLane.dir,
                                                        self.lDirSteps[1]))
            iLane.clean(self.adapters["adpR1"], self.adapters["adpR2"],
                        self.errTol, self.minOvl, self.minReadLen,
                        self.minQual, self.maxNPerc, stepDir, iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def launchFastqcOnCleanFastqFiles(self):
        if self.verbose > 0:
            msg = "assess quality per sample..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step3-fastqc" % self.project1Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            stepDir = "%s/%s" % (iLane.dir, self.lDirSteps[2])
            iLane.setCleanedFastqFiles(stepDir)
            iLane.cleanQc(stepDir, iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def saveNbReadsFromFastqc(self):
        if self.verbose > 0:
            msg = "save nb of clean reads per sample..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        outFile = "%s/%s_clean-reads-per-sample.txt.gz" % (self.allLanesDir,
                                                           self.project1Id)
        if os.path.exists(outFile):
            if self.forceRerunSteps:
                os.remove(outFile)
        outHandle = gzip.open(outFile, "w")
        
        txt = "geno"
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
            
        groupJobId = self.makeGroupJobId("%s_gbs-step4-align" % self.project2Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            
            # create this directory here, even if it is used only later,
            # in gatherSamplesPerLane, so that the user is immediately
            # warned if the directory already exists, before launching jobs
            outDir = "%s/%s" % (self.allSamplesDir, iLane.id)
            self.makeDir(outDir)
            
            iLane.setCleanedFastqFiles("%s/%s" % (iLane.dir, self.lDirSteps[2]))
            iLane.align(self.pathToPrefixRefGenome, self.tmpDir, self.dictFile,
                        self.lDirSteps[3], iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def saveBasicAlignStats(self):
        """
        via samtools flagstat per sample
        """
        if self.verbose > 0:
            msg = "save basic alignment statistics..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        outFile = "%s/%s_basic-align-stats.txt.gz" % \
                  (self.allSamplesDir, self.project2Id)
        if os.path.exists(outFile):
            if self.forceRerunSteps:
                os.remove(outFile)
        outHandle = gzip.open(outFile, "w")
        
        txt = "geno"
        txt += "\tflowcell"
        txt += "\tlane"
        txt += "\t%s" % SamtoolsFlagstat.header2str()
        outHandle.write("%s\n" % txt)
        
        lLanes = self.dLanes.keys()
        lLanes.sort()
        for laneId in lLanes:
            iLane = self.dLanes[laneId]
            iLane.setInitialBamFiles(self.lDirSteps[3])
            iLane.saveSamtoolsFlagstat(outHandle)
            
        outHandle.close()
        
        outFile = "%s/%s_reads-count-ok.txt" % \
                  (self.allSamplesDir, self.project2Id)
        if os.path.exists(outFile):
            os.remove(outFile)
        if os.path.exists("%s.gz" % outFile):
            os.remove("%s.gz" % outFile)
        outHandle = open(outFile, "w")
        txt = "geno"
        txt += "\tflowcell"
        txt += "\tlane"
        txt += "\tnb.reads.align.ok" # nb of reads if single, of pairs if paired-end
        outHandle.write("%s\n" % txt)
        for laneId in self.dLanes:
            for sampleId in self.dLanes[laneId].dSamples:
                tmpFile = "%s/%s/align-reads/reads_count-ok_%s.txt" \
                          % (self.allSamplesDir, sampleId, sampleId)
                tmpHandle = open(tmpFile, "rb")
                outHandle.write(tmpHandle.read())
                tmpHandle.close()
        outHandle.close()
        with open(outFile, "rb") as f_in, \
             gzip.open("%s.gz" % outFile, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        if os.path.exists(outFile):
            os.remove(outFile)
            
            
    def gatherSamplesPerLane(self):
        """
        merge, sort and index BAM files, and collect insert sizes
        """
        if self.verbose > 0:
            msg = "gather samples per lane..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step4-gather" % self.project2Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue2, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.allSamplesDir, iLane.id)
            iLane.gather(self.jvmXms, self.jvmXmx, self.tmpDir,
                         outDir, iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step4(self):
        cwd, startTime = self.beginStep(4)
        self.alignCleanedReads()
        self.saveBasicAlignStats()
        self.gatherSamplesPerLane()
        self.endStep(4, cwd, startTime)
        
        
    def localRealignmentSamples(self):
        if self.verbose > 0:
            msg = "locally realign reads per sample..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step5-realign-samples" %
                                         self.project2Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            iLane.setInitialBamFiles(self.lDirSteps[3]) # from previous step
            iLane.localRealignSamples(self.jvmXms, self.jvmXmx,
                                      self.pathToPrefixRefGenome,
                                      self.knownIndelsFile, self.lDirSteps[4],
                                      iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step5(self):
        cwd, startTime = self.beginStep(5)
        self.localRealignmentSamples()
        self.endStep(5, cwd, startTime)
        
        
    def localRealignmentGenotypes(self):
        if self.verbose > 0:
            msg = "locally realign reads per genotype..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step6-realign-geno" %
                                         self.project2Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        lGenoIds = self.dGenos.keys()
        lGenoIds.sort()
        for i,genoId in enumerate(lGenoIds):
            iGeno = self.dGenos[genoId]
            stepDir = "%s/%s" % (iGeno.dir, self.lDirSteps[5])
            self.makeDir(stepDir)
            iGeno.setRealignedSampleBamFiles(self.lDirSteps[4])
            iGeno.localRealign(self.jvmXms, self.jvmXmx,
                               self.pathToPrefixRefGenome,
                               self.knownIndelsFile, stepDir,
                               iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step6(self):
        cwd, startTime = self.beginStep(6)
        self.localRealignmentGenotypes()
        self.endStep(6, cwd, startTime)
        
        
    def variantCallingPerGenotype(self):
        if self.verbose > 0:
            msg = "call variants per genotype..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_gbs-step7-varcall-geno" \
                                         % self.project2Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        lGenoIds = self.dGenos.keys()
        lGenoIds.sort()
        for i,genoId in enumerate(lGenoIds):
            iGeno = self.dGenos[genoId]
            stepDir = "%s/%s" % (iGeno.dir, self.lDirSteps[6])
            self.makeDir(stepDir)
            iGeno.setRealignedGenotypeBamFile("%s/%s" % (iGeno.dir,
                                                         self.lDirSteps[5]))
            iGeno.variantCalling(self.jvmXms, self.jvmXmx, self.tmpDir,
                                 self.pathToPrefixRefGenome, self.knownFile,
                                 stepDir, iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step7(self):
        cwd, startTime = self.beginStep(7)
        self.variantCallingPerGenotype()
        self.endStep(7, cwd, startTime)
        
        
    def jointGenotyping(self):
        if self.verbose > 0:
            msg = "combine variant calls across genotypes..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_%s_gbs-step8-joincall" \
                                         % (self.project2Id,
                                            self.jointGenoId))
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        stepDir = "%s/%s_%s" % (self.allGenosDir, self.lDirSteps[7],
                                self.jointGenoId)
        self.makeDir(stepDir)
        if self.verbose > 0:
            msg = "results will be in '%s'" % stepDir
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        cmd = "echo \"commands generated by gbs.py %s\"" % progVersion
        cmd += "\njava"
        cmd += " -Xms%s" % self.jvmXms
        cmd += " -Xmx%s" % self.jvmXmx
        # if self.tmpDir:
        #     cmd += " -Djava.io.tmpdir=%s" % self.tmpDir
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T GenotypeGVCFs"
        cmd += " -R %s.fa" % self.pathToPrefixRefGenome
        lGenoIds = self.dGenos.keys()
        lGenoIds.sort()
        for i,genoId in enumerate(lGenoIds):
            iGeno = self.dGenos[genoId]
            prevStepDir = "%s/%s" % (iGeno.dir, self.lDirSteps[6])
            cmd += " --variant %s/%s.g.vcf.gz" % (prevStepDir, genoId)
        cmd += " -o %s/%s_%s_raw.vcf.gz" % (stepDir, self.project2Id,
                                            self.jointGenoId)
        cmd += " --heterozygosity 0.001" # 1 het site in 100 bp (across all samples, for humans)
        if self.verbose > 1:
            print(cmd)
        jobName = "stdout_%s" % (iJobGroup.id)
        bashFile = "%s/job_%s_%s.bash" % (stepDir, iJobGroup.id, "jointcalling")
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=stepDir)
        iJobGroup.insert(iJob)
        
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step8(self):
        cwd, startTime = self.beginStep(8)
        self.jointGenotyping()
        self.endStep(8, cwd, startTime)
        
        
    def variantFiltering(self):
        """
        http://www.nature.com/articles/sdata201511
        """
        if self.verbose > 0:
            msg = "filter variant calls across genotypes..."
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
            
        groupJobId = self.makeGroupJobId("%s_%s_gbs-step9-filter" \
                                         % (self.project2Id,
                                            self.jointGenoId))
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        stepDir = "%s/%s_%s" % (self.allGenosDir, self.lDirSteps[7],
                                self.jointGenoId)
        if not os.path.exists(stepDir):
            msg = "dir %s doesn't exist\n" % stepDir
            msg += "\nlaunch step 8 first"
            raise OSError(msg)
        if self.verbose > 0:
            msg = "results will be in '%s'" % stepDir
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()

        # 1) SelectVariants -selectType SNP
        cmd = "echo \"commands generated by gbs.py %s\"" % progVersion
        cmd += "\njava"
        cmd += " -Xms%s" % self.jvmXms
        cmd += " -Xmx%s" % self.jvmXmx
        # if self.tmpDir:
        #     cmd += " -Djava.io.tmpdir=%s" % self.tmpDir
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T SelectVariants"
        cmd += " -R %s.fa" % self.pathToPrefixRefGenome
        cmd += " -V %s/%s_%s_raw.vcf.gz" % (stepDir, self.project2Id,
                                            self.jointGenoId)
        cmd += " -selectType SNP"
        cmd += " -o %s/%s_%s_raw-snps.vcf.gz" % (stepDir, self.project2Id,
                                                 self.jointGenoId)
        
        # 2) VariantFiltration --filterExpression --genotypeFilterExpression
        cmd += "\n"
        cmd += "java"
        cmd += " -Xms%s" % self.jvmXms
        cmd += " -Xmx%s" % self.jvmXmx
        # if self.tmpDir:
        #     cmd += " -Djava.io.tmpdir=%s" % self.tmpDir
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T VariantFiltration"
        cmd += " -R %s.fa" % self.pathToPrefixRefGenome
        cmd += " -V %s/%s_%s_raw-snps.vcf.gz" % (stepDir, self.project2Id,
                                                 self.jointGenoId)
        cmd += " --clusterSize 3" # default value from GATK
        cmd += " --clusterWindowSize 0" # default value from GATK
        cmd += " --filterName \"high_DP\""
        cmd += " --filterExpression \"DP > 100000\""
        cmd += " --filterName \"low_QD\""
        cmd += " --filterExpression \"QD < 2.0\""
        cmd += " --filterName \"high_FS\""
        cmd += " --filterExpression \"FS > 60.0\""
        cmd += " --filterName \"low_MQ\""
        cmd += " --filterExpression \"MQ < 40.0\""
        cmd += " --filterName \"low_MQRankSum\""
        cmd += " --filterExpression \"MQRankSum < -12.5\""
        cmd += " --filterName \"low_RPRS\""
        cmd += " --filterExpression \"ReadPosRankSum < -8.0\""
        cmd += " --filterName \"high_HS\""
        cmd += " --filterExpression \"HaplotypeScore > 13.0\""
        cmd += " --filterName \"high_SOR\""
        cmd += " --filterExpression \"SOR > 4.0\""
        if self.minDp:
            cmd += " --genotypeFilterName \"low_DP\""
            cmd += " --genotypeFilterExpression \"DP < %i\"" % self.minDp
        if self.minGq:
            cmd += " --genotypeFilterName \"low_GQ\""
            cmd += " --genotypeFilterExpression \"GQ < %i\"" % self.minGq
        cmd += " -o %s/%s_%s_raw-snps_tmp1.vcf.gz" % (stepDir, self.project2Id,
                                                     self.jointGenoId)
        
        # 3) SelectVariants <exclude sites> --setFilteredGtToNocall -ped ...
        cmd += "\n"
        cmd += "java"
        cmd += " -Xms%s" % self.jvmXms
        cmd += " -Xmx%s" % self.jvmXmx
        # if self.tmpDir:
        #     cmd += " -Djava.io.tmpdir=%s" % self.tmpDir
        cmd += " -jar `which GenomeAnalysisTK.jar`"
        cmd += " -T SelectVariants"
        cmd += " -R %s.fa" % self.pathToPrefixRefGenome
        cmd += " -V %s/%s_%s_raw-snps_tmp1.vcf.gz" % (stepDir, self.project2Id,
                                                     self.jointGenoId)
        cmd += " --restrictAllelesTo %s" % self.restrictAllelesTo
        cmd += " --excludeFiltered"
        cmd += " --excludeNonVariants"
        cmd += " --setFilteredGtToNocall"
        if self.famFile:
            cmd += " -ped %s" % self.famFile
            cmd += " -mv -mvq %i -invMv" % self.mendelianViolationQualThreshold
            cmd += " --pedigreeValidationType SILENT"
        if self.excludeSampleFile:
            cmd += " --exclude_sample_file %s" % self.excludeSampleFile
        cmd += " -o %s/%s_%s_raw-snps_tmp2.vcf.gz" % (stepDir,
                                                      self.project2Id,
                                                      self.jointGenoId)
        
        # 4) SelectVariants <FilteredGenotypes> <NOCALL>
        # http://gatkforums.broadinstitute.org/gatk/discussion/9795
        if self.maxNbFilterGenos or self.maxFracFilterGenos or \
           self.maxNbNocallGenos or self.maxFracNocallGenos:
            cmd += "\n"
            cmd += "java"
            cmd += " -Xms%s" % self.jvmXms
            cmd += " -Xmx%s" % self.jvmXmx
            # if self.tmpDir:
            #     cmd += " -Djava.io.tmpdir=%s" % self.tmpDir
            cmd += " -jar `which GenomeAnalysisTK.jar`"
            cmd += " -T SelectVariants"
            cmd += " -R %s.fa" % self.pathToPrefixRefGenome
            cmd += " -V %s/%s_%s_raw-snps_tmp2.vcf.gz" % (stepDir, self.project2Id,
                                                          self.jointGenoId)
            if self.maxNbFilterGenos:
                cmd += " --maxFilteredGenotypes %i" % self.maxNbFilterGenos
            if self.maxFracFilterGenos:
                cmd += " --maxFractionFilteredGenotypes %f" % self.maxFracFilterGenos
            if self.maxNbNocallGenos:
                cmd += " --maxNOCALLnumber %i" % self.maxNbNocallGenos
            if self.maxFracNocallGenos:
                cmd += " --maxNOCALLfraction %f" % self.maxFracNocallGenos
            cmd += " -o %s/%s_%s_raw-snps_gatk-filter.vcf.gz" % (stepDir,
                                                                 self.project2Id,
                                                                 self.jointGenoId)
            cmd += "\n"
            cmd += "rm %s/%s_%s_raw-snps_tmp2.vcf.gz*" % (stepDir, self.project2Id,
                                                          self.jointGenoId)
        else:
            cmd += "\n"
            cmd += "mv %s/%s_%s_raw-snps_tmp2.vcf.gz" % (stepDir,
                                                         self.project2Id,
                                                         self.jointGenoId)
            cmd += " %s/%s_%s_raw-snps_gatk-filter.vcf.gz" % (stepDir,
                                                              self.project2Id,
                                                              self.jointGenoId)
            cmd += "\n"
            cmd += "mv %s/%s_%s_raw-snps_tmp2.vcf.gz.tbi" % (stepDir,
                                                             self.project2Id,
                                                             self.jointGenoId)
            cmd += " %s/%s_%s_raw-snps_gatk-filter.vcf.gz.tbi" % (stepDir,
                                                                  self.project2Id,
                                                                  self.jointGenoId)
        cmd += "\n"
        cmd += "rm %s/%s_%s_raw-snps_tmp1.vcf.gz*" % (stepDir, self.project2Id,
                                                      self.jointGenoId)
        
        if self.verbose > 1:
            print(cmd)
        jobName = "stdout_%s" % (iJobGroup.id)
        bashFile = "%s/job_%s_%s.bash" % (stepDir, iJobGroup.id,
                                          "select-filter-snps")
        iJob = Job(groupId=iJobGroup.id, name=jobName, cmd=cmd,
                   bashFile=bashFile, dir=stepDir)
        iJobGroup.insert(iJob)
        
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step9(self):
        cwd, startTime = self.beginStep(9)
        self.variantFiltering()
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
            
        groupJobId = self.makeGroupJobId("%s_gbs-step11-recalib" \
                                         % self.project2Id)
        if self.verbose > 0:
            msg = "groupJobId=%s" % groupJobId
            sys.stdout.write("%s\n" % msg)
            sys.stdout.flush()
        iJobGroup = JobGroup(groupJobId, self.queue, self.lResources)
        self.jobManager.insert(iJobGroup)
        
        for laneId,iLane in self.dLanes.items():
            outDir = "%s/%s" % (self.lDirSteps[4], laneId)
            iLane.baseQualityRecalibrate(self.jvmXms, self.jvmXmx,
                                         self.pathToPrefixRefGenome,
                                         self.knownFile, outDir, iJobGroup)
            
        self.jobManager.submit(iJobGroup.id)
        self.jobManager.wait(iJobGroup.id, self.rmvBash, self.verbose)
        
        
    def step11(self):
        cwd, startTime = self.beginStep(11)
        self.baseQualityRecalibration()
        self.endStep(11, cwd, startTime)
        
        
    def run(self):
        self.loadSamplesFile()
        
        # set up directories
        if "1" in self.lSteps or "2" in self.lSteps or "3" in self.lSteps \
           or "4" in self.lSteps:
            self.setupLaneDirectories()
        if "4" in self.lSteps or "5" in self.lSteps or "6" in self.lSteps:
            self.setupSampleDirectories()
        if "6" in self.lSteps or "7" in self.lSteps or "8" in self.lSteps or \
           "9" in self.lSteps:
            self.setupGenotypeDirectories()
            
        # set up job manager
        if "1" in self.lSteps or "2" in self.lSteps or "3" in self.lSteps:
            self.jobManager = JobManager(self.scheduler, self.project1Id)
        if "4" in self.lSteps or "5" in self.lSteps or "6" in self.lSteps \
           or "7" in self.lSteps or "8" in self.lSteps or "9" in self.lSteps:
            self.jobManager = JobManager(self.scheduler, self.project2Id)
            
        # execute the step(s)
        if "1" in self.lSteps: # init read quality
            self.step1()
            
        if "2" in self.lSteps: # demultiplex reads
            self.step2()
            
        if "3" in self.lSteps: # clean reads
            self.step3()
            
        if "4" in self.lSteps: # align per sample
            self.step4()
            
        if "5" in self.lSteps: # realign per sample
            self.step5()
            
        if "6" in self.lSteps: # realign per geno
            self.step6()
            
        if "7" in self.lSteps: # var and geno call per geno
            self.step7()
            
        if "8" in self.lSteps: # var and geno call jointly
            self.step8()
            
        if "9" in self.lSteps: # var filter
            self.step9()
            
        self.jobManager.close()
        
        
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
