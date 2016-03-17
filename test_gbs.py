#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: test gbs.py
# Copyright (C) 2015-2016 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

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
import tempfile
import shutil
import itertools
import random
import hashlib

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
        
progVersion = "0.3.0" # http://semver.org/


class TestGbs(object):
    
    def __init__(self):
        self.verbose = 0
        self.pathToProg = ""
        self.testsToRun = ["1"]
        self.queue = "normal.q"
        self.lResources = None
        self.clean = True
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' tests gbs.py.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (default=0/1/2/3)\n"
        msg += "  -p, --p2p\tfull path to the program to be tested\n"
        msg += "  -t, --test\tidentifiers of test(s) to run (default=1)\n"
        msg += "  -q, --queue\tname of the cluster queue (default=normal.q)\n"
        msg += "  -r, --resou\tcluster resources (e.g. 'test' for 'qsub -l test')\n"
        msg += "  -n, --noclean\tkeep temporary directory with all files\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -p ~/src/gbs.py\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        
        The person roles complies with R's guidelines (The R Journal Vol. 4/1, June 2012).
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2015-2016 Institut National de la Recherche Agronomique (INRA).\n"
        msg += "License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:p:t:q:r:n",
                                        ["help", "version", "verbose=",
                                         "p2p=", "test=", "queue=", "resou=",
                                         "noclean"])
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
            elif o == "-p" or o == "--p2p":
                 self.pathToProg = a
            elif o == "-t" or o == "--test":
                self.testsToRun = a.split("-")
            elif o == "-q" or o == "--queue":
                self.queue = a
            elif o == "-r" or o == "--resou":
                self.lResources = a.split()
            elif o == "-n" or o == "--noclean":
                self.clean = False
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if self.pathToProg == "":
            msg = "ERROR: missing compulsory option --p2p"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.pathToProg):
            msg = "ERROR: can't find file %s" % self.pathToProg
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        for t in self.testsToRun:
            if t not in ["1"]:
                msg = "ERROR: unknown --test %s" % t
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
                
                
    #==========================================================================
    
    
    def beforeTest(self):
        """
        Create a temporary folder with a unique ID
        """
        cwd = os.getcwd()
        testDir = tempfile.mkdtemp(dir=cwd, prefix="tmp_test_")
        os.chdir(testDir)
        if self.verbose > 0:
            print("temp dir: %s" % os.getcwd()); sys.stdout.flush()
        return cwd, testDir
        
        
    def makeSamplesFile(self):
        """
        4 individuals: ind1 (mother), ind2 (father), ind3 (child), ind4 (child)
        2 lanes, 6 samples (3 individuals per lane, parents present in both)
        2 extractions, and thus 2 libraries, for ind1
        """
        samplesFile = "samples.txt"
        lane2samples = {}
        
        samplesHandle = open(samplesFile, "w")
        
        lane2samples["1"] = {"R1": "init_reads_lane1_R1.fastq.gz",
                          "R2": "init_reads_lane1_R2.fastq.gz",
                          "inds": {"ind1": {"gen":0, "lib":"ind1-A", "tag":"AAAA"},
                                   "ind2": {"gen":0, "lib":"ind2", "tag":"GGGG"},
                                   "ind3": {"gen":1, "lib":"ind3", "tag":"TTTT"}}}
        lane2samples["2"] = {"R1": "init_reads_lane2_R1.fastq.gz",
                          "R2": "init_reads_lane2_R2.fastq.gz",
                          "inds": {"ind1": {"gen":0, "lib":"ind1-B", "tag":"CCCC"},
                                   "ind2": {"gen":0, "lib":"ind2", "tag":"TTTT"},
                                   "ind4": {"gen":1, "lib":"ind4", "tag":"AAAA"}}}
        
        # header
        txt = "individual"
        txt += "\tspecies"
        txt += "\tlibrary"
        txt += "\tbarcode"
        txt += "\tseq_center"
        txt += "\tseq_platform"
        txt += "\tseq_platform_model"
        txt += "\tflowcell"
        txt += "\tlane"
        txt += "\tdate"
        txt += "\tfastq_file_R1"
        txt += "\tfastq_file_R2"
        
        # ind1 in lane 1
        txt += "\n"
        txt += "ind1"
        txt += "\tTest example"
        txt += "\t%s" % lane2samples["1"]["inds"]["ind1"]["lib"]
        txt += "\t%s" % lane2samples["1"]["inds"]["ind1"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 1
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["1"]["R1"]
        txt += "\t%s" % lane2samples["1"]["R2"]
        
        # ind2 in lane 1
        txt += "\n"
        txt += "ind2"
        txt += "\tTest example"
        txt += "\t%s" % lane2samples["1"]["inds"]["ind2"]["lib"]
        txt += "\t%s" % lane2samples["1"]["inds"]["ind2"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%s" % 1
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["1"]["R1"]
        txt += "\t%s" % lane2samples["1"]["R2"]
        
        # ind3 in lane 1
        txt += "\n"
        txt += "ind3"
        txt += "\tTest example"
        txt += "\t%s" % lane2samples["1"]["inds"]["ind3"]["lib"]
        txt += "\t%s" % lane2samples["1"]["inds"]["ind3"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%s" % 1
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["1"]["R1"]
        txt += "\t%s" % lane2samples["1"]["R2"]
        
        # ind1 in lane 2
        txt += "\n"
        txt += "ind1"
        txt += "\tTest example"
        txt += "\t%s" % lane2samples["2"]["inds"]["ind1"]["lib"]
        txt += "\t%s" % lane2samples["2"]["inds"]["ind1"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 2
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["2"]["R1"]
        txt += "\t%s" % lane2samples["2"]["R2"]
        
        # ind2 in lane 2
        txt += "\n"
        txt += "ind2"
        txt += "\tTest example"
        txt += "\t%s" % lane2samples["2"]["inds"]["ind2"]["lib"]
        txt += "\t%s" % lane2samples["2"]["inds"]["ind2"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 2
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["2"]["R1"]
        txt += "\t%s" % lane2samples["2"]["R2"]
        
        # ind4 in lane 2
        txt += "\n"
        txt += "ind4"
        txt += "\tTest example"
        txt += "\t%s" % lane2samples["2"]["inds"]["ind4"]["lib"]
        txt += "\t%s" % lane2samples["2"]["inds"]["ind4"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 2
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["2"]["R1"]
        txt += "\t%s" % lane2samples["2"]["R2"]
        
        samplesHandle.write("%s\n" % txt)
        samplesHandle.close()
        return samplesFile, lane2samples
    
    
    def simulRefGenomeWoMotif(self, nbChrs, lenChr):
        """
        >>> i = TestGbs()
        >>> obs = i.simulRefGenomeWoMotif(2, 100)
        >>> type(obs) == type([])
        True
        >>> len(obs)
        2
        >>> [len(x) for x in obs]
        [100, 100]
        """
        refGen = []
        for chrNum in range(nbChrs):
            refGen.append(SeqRecord(Seq("".join([random.choice(["A","T","G","C"]) \
                                                 for _ in xrange(lenChr)]),
                                        generic_dna),
                                    id="chr%i" % (chrNum + 1),
                                    name="chr%i" % (chrNum + 1),
                                    description="Texample|v1|chr%i" % (chrNum + 1)))
        return refGen
        
        
    def chooseFragCoords(self, nbChrs, lenChr, nbFragsPerChr, minFragLen,
                         maxFragLen):
        """
        >>> i = TestGbs()
        >>> obs = i.chooseFragCoords(2, 1000, 4, 5, 20)
        >>> type(obs) == type([])
        True
        >>> len(obs)
        2
        >>> len(obs[0])
        4
        >>> [type(x) == type([]) for x in obs]
        [True, True]
        >>> [len(x) for x in obs[0]]
        [2, 2, 2, 2]
        """
        lFragCoords = []
        for chrNum in range(nbChrs):
            lStarts = range(maxFragLen + 1, lenChr - (maxFragLen + 1),
                            int(math.ceil((lenChr - 2 * (maxFragLen + 1)) \
                                          / float(nbFragsPerChr))))
            lSizes = [random.choice(xrange(minFragLen, maxFragLen)) \
                      for _ in xrange(len(lStarts))]
            lCoords = [[i[0], i[0]+i[1]] for i in zip(lStarts,
                                                       lSizes)]
            lFragCoords.append(lCoords)
        return lFragCoords
        
        
    def insertMotifsPerFrag(self, refGen, lFragCoords, motif):
        """
        >>> i = TestGbs()
        >>> refGen = [SeqRecord(Seq("AAAAAAAAAAAAAAAAAAAA"), id="test")]
        >>> lFragCoords = [ [ [2,10] ] ]
        >>> obs = i.insertMotifsPerFrag(refGen, lFragCoords, "TT")
        >>> str(obs[0].seq)
        'ATTAAAAATTAAAAAAAAAA'
        """
        nbChrs = len(refGen)
        for c in range(nbChrs):
            tmp = str(refGen[c].seq)
            for fragCoord in lFragCoords[c]:
                tmp = tmp[0:(fragCoord[0]-1)] \
                      + motif \
                      + tmp[(fragCoord[0]-1+len(motif)):(fragCoord[1]-len(motif))] \
                      + motif \
                      + tmp[fragCoord[1]:]
            refGen[c].seq = Seq(tmp, generic_dna)
        return refGen
        
        
    def simulRefGenome(self, nbChrs, lenChr, nbFragsPerChr, minFragLen,
                       maxFragLen, motif):
        """
        the reference genome is haploid
        motif -- string with what remains of the restriction site after being cut
        """
        refGen = self.simulRefGenomeWoMotif(nbChrs, lenChr)
        lFragCoords = self.chooseFragCoords(nbChrs, lenChr, nbFragsPerChr,
                                            minFragLen, maxFragLen)
        refGen = self.insertMotifsPerFrag(refGen, lFragCoords, motif)
        return refGen, lFragCoords
        
        
    def simulParentGenome(self, refGen):
        """
        diploid
        TODO: add SNPs to ref genome (no indels)
        """
        return refGen
        
        
    def simulChildGenome(self, genMother, genFather):
        """
        diploid
        TODO: add crossing-overs to parental genomes
        """
        return genMother
        
        
    def simulIndGenomes(self, lane2samples, refGen):
        ind2genome = {}
        
        parents = []
        for lane in lane2samples:
            for ind in lane2samples[lane]["inds"]:
                if ind not in ind2genome and \
                   lane2samples[lane]["inds"][ind]["gen"] == 0:
                    ind2genome[ind] = self.simulParentGenome(refGen)
                    parents.append(ind)
        assert len(parents) == 2
        
        for lane in lane2samples:
            for ind in lane2samples[lane]["inds"]:
                if ind not in ind2genome and \
                   lane2samples[lane]["inds"][ind]["gen"] == 1:
                    ind2genome[ind] \
                        = self.simulChildGenome(ind2genome[parents[0]],
                                                ind2genome[parents[1]])
                    
        return ind2genome
        
        
    def simulInserts(self, lane2samples, lFragCoords):
        """
        TODO: some fragments could not be sequenced, in all or some individuals
        """
        lane2inserts = {}
        for lane in lane2samples:
            lane2inserts[lane] = {}
            for ind in lane2samples[lane]["inds"]:
                lane2inserts[lane][ind] = lFragCoords
        return lane2inserts
        
        
    def extractReadSequences(self, indGenChr, lInsertsChr, barcode, lenRead):
        """
        >>> i = TestGbs()
        >>> indGenChr = SeqRecord(Seq("AATTTAGGGA"), id="chr1")
        >>> i.extractReadSequences(indGenChr, [[3,9]], "C", 3)
        [['CTT', 'CCC']]
        """
        return [[str(barcode) \
                 + str(indGenChr.seq[(i[0]-1):(i[0]-1+lenRead-len(barcode))]),
                 reverse_complement(str(indGenChr.seq[(i[1]-lenRead):i[1]]))]
                for i in lInsertsChr]
        
        
    def simulPairedReadsPerSample(self, lane, barcode, indGen, lInserts,
                                  nbReads, lenRead, readCounter):
        """
        barcode -- string
        indGen -- [SeqRecChr1, SeqRecChr2, ...]
        lInserts -- [listChr1[], listChr2[], ...]
        TODO: use nbReads
        """
        lR1 = []
        lR2 = []
        nbChrs = len(indGen)
        for c in range(nbChrs):
            lInsertPairs = self.extractReadSequences(indGen[c], lInserts[c],
                                                     barcode, lenRead)
            lR1 = lR1 + [SeqRecord(Seq(lInsertPair[0], generic_dna),
                                   id="lane%s:read%i" % (lane,
                                                         readCounter+i+1),
                                   description="1:N:0") \
                         for i,lInsertPair in enumerate(lInsertPairs)]
            lR2 = lR2 + [SeqRecord(Seq(lInsertPair[1], generic_dna),
                                   id="lane%s:read%i" % (lane,
                                                         readCounter+i+1),
                                   description="2:N:0") \
                         for i,lInsertPair in enumerate(lInsertPairs)]
            readCounter += len(lInsertPairs)
        return lR1, lR2, readCounter
        
        
    def simulPairedReadsPerLane(self, lane, samples, ind2gen, lInserts, nbReads,
                                lenRead):
        """
        lane -- string
        samples -- dict[ind]{gen:..., tag:...}
        ind2gen -- dict[ind][SeqRecChr1, SeqRecChr2, ...]
        lInserts -- dict{ind}[ ... ]
        """
        lR1 = []
        lR2 = []
        readCounter = 0
        for ind in samples:
            lSampleR1, lSampleR2, readCounter \
                = self.simulPairedReadsPerSample(lane,
                                                 samples[ind]["tag"],
                                                 ind2gen[ind],
                                                 lInserts[ind],
                                                 nbReads,
                                                 lenRead,
                                                 readCounter)
            lR1 = lR1 + lSampleR1
            lR2 = lR2 + lSampleR2
        return lR1, lR2
        
        
    def simulReads(self, lane2samples, ind2gen, lane2inserts, nbReads, lenRead):
        """
        lane2samples -- dict[lane]{R1:file, R2:file, samples={ind:{gen, tag}, ...}}
        ind2gen -- dict[ind][SeqRecChr1, SeqRecChr2, ...]
        lane2inserts: dict{lane}{ind}[listChr1[], listChr2[], ...]
        """
        lane2reads = {}
        for lane in lane2samples:
            lane2reads[lane] = {"R1": [], "R2": []}
            lane2reads[lane]["R1"], lane2reads[lane]["R2"] \
                = self.simulPairedReadsPerLane(lane,
                                               lane2samples[lane]["inds"],
                                               ind2gen,
                                               lane2inserts[lane],
                                               nbReads,
                                               lenRead)
        return lane2reads
        
        
    def saveRefGenome(self, refGen):
        refFile = "refgenome.fa"
        refHandle = open(refFile, "w")
        for chrRec in refGen:
            refHandle.write(chrRec.format("fasta"))
        refHandle.close()
        return refFile
        
        
    def saveReads(self, lane2samples, lane2reads, paired):
        for lane in lane2samples:
            fastqHandle = gzip.open(lane2samples[lane]["R1"], "w")
            for rec in lane2reads[lane]["R1"]:
                rec.letter_annotations["phred_quality"] = [40] * len(rec)
                fastqHandle.write(rec.format("fastq"))
            fastqHandle.close()
            if paired:
                fastqHandle = gzip.open(lane2samples[lane]["R2"], "w")
                for rec in lane2reads[lane]["R2"]:
                    rec.letter_annotations["phred_quality"] = [40] * len(rec)
                    fastqHandle.write(rec.format("fastq"))
                fastqHandle.close()
                
                
    def makeInputSequences(self, lane2samples, nbChrs=2, lenChr=100000, nbFragsPerChr=50,
                           minFragLen=30, maxFragLen=600, motif="CAGC",
                           nbReads=1000, lenRead=100, paired=True):
        refGen, lFragCoords \
            = self.simulRefGenome(nbChrs, lenChr, nbFragsPerChr, minFragLen,
                                  maxFragLen, motif)
        
        ind2gen = self.simulIndGenomes(lane2samples, refGen)
        
        lane2inserts = self.simulInserts(lane2samples, lFragCoords)
        
        lane2reads = self.simulReads(lane2samples, ind2gen, lane2inserts, nbReads, lenRead)
        
        refFile = self.saveRefGenome(refGen)
        self.saveReads(lane2samples, lane2reads, paired)
        return refGen, refFile
    
    
    def makeAdapterFile(self):
        """
        Illumina's universal adapters, written from 5' (left) to 3' (right).
        """
        adpFile = "adapters.txt"
        dAdps = {}
        
        dAdps["adpR1"] = "CTCTTCCGATCT"
        dAdps["adpR2"] = dAdps["adpR1"]
        
        adpHandle = open(adpFile, "w")
        for adp,seq in dAdps.items():
            adpHandle.write("%s\t%s\n" % (adp, seq))
        adpHandle.close()
        
        return adpFile, dAdps
    
    
    def makeDictFile(self, refGen, refFile):
        dictFile = "refgenome.dict"
        
        dictHandle = open(dictFile, "w")
        dictHandle.write("@HD\tVN:1.4\tSO:unsorted\n")
        for c in range(len(refGen)):
            txt = "@SQ"
            txt += "\tSN:%s" % refGen[c].id
            txt += "\tLN:%i" % len(refGen[c])
            txt += "\tSP:%s" % "Test example"
            txt += "\tAS:%s" % "v1"
            txt += "\tM5:%s" % hashlib.md5(str(refGen[c].seq)).hexdigest()
            txt += "\tUR:file:%s" % refFile
            dictHandle.write("%s\n" % txt)
        dictHandle.close()
        
        return dictFile
    
    
    def makeBwaIndexFiles(self, refFile):
        args = ["bwa", "index", "-p", refFile.split(".fa")[0], refFile]
        stdoutFilehandle = open("bwa_index.log", "w")
        subprocess.check_call(args, stdout=stdoutFilehandle,
                              stderr=subprocess.STDOUT)
        stdoutFilehandle.close()
        
        
    def makeFaIndexFile(self, refFile):
        args = ["samtools", "faidx", refFile]
        stdoutFilehandle = open("samtools_faidx.log", "w")
        subprocess.check_call(args, stdout=stdoutFilehandle,
                              stderr=subprocess.STDOUT)
        stdoutFilehandle.close()
        
        
    def launchProg(self, options):
        """
        Launch gbs.py.
        """
        args = [self.pathToProg] \
               + options \
               + ["-v", str(self.verbose - 1)]
        if self.verbose > 0:
            print(" ".join(args))
        msgs = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE).communicate()
        if self.verbose > 1:
            if msgs[0] != "":
                print("stdout:")
                print(msgs[0][:-1])
            if msgs[1] != "":
                print("sterr:")
                print(msgs[1][:-1])
        return msgs
        
        
    def afterTest(self, cwd, testDir):
        os.chdir(cwd)
        if self.clean:
            shutil.rmtree(testDir)
            
            
    #==========================================================================
    
    
    def test_step1_prepare(self):
        options = []
        
        options.append("--proj")
        options.append("testGbs")
        options.append("--queue")
        options.append(self.queue)
        options.append("--resou")
        options.append(" ".join(self.lResources))
        options.append("--step")
        options.append("1")
        
        samplesFile, lane2samples = self.makeSamplesFile()
        nbChrs = 2
        lenChr = 100000
        nbFragsPerChr = 50
        minFragLen = 30
        maxFragLen = 300
        motif = "CAGC"
        nbReads = 1000 # for a given sample, over all its chrs
        lenRead = 100
        refGen, refFile = self.makeInputSequences(lane2samples, nbChrs, lenChr,
                                                  nbFragsPerChr, minFragLen,
                                                  maxFragLen, motif,
                                                  nbReads, lenRead)
        adpFile, dAdps = self.makeAdapterFile()
        dictFile = self.makeDictFile(refGen, refFile)
        self.makeBwaIndexFiles(refFile)
        self.makeFaIndexFile(refFile)
        options.append("--samples")
        options.append(samplesFile)
        options.append("--pird")
        options.append(os.getcwd())
        
        return options
        
        
    def test_step1_comp(self, msgs):
        pass
        
        
    def test_step1(self):
        if self.verbose > 0:
            print("launch test step 1 ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        options = self.test_step1_prepare()
        msgs = self.launchProg(options)
        self.test_step1_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================
    
    
    def run(self):
        if "1" in self.testsToRun:
            self.test_step1()
            
            
if __name__ == "__main__":
    i = TestGbs()
    
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
            p = subprocess.Popen(["grep", "VmHWM", "/proc/%s/status" % \
                                  os.getpid()],
                                 shell=False,
                                 stdout=subprocess.PIPE).communicate()
            maxMem = p[0].split()[1]
            msg += "; %s kB)" % maxMem
        else:
            msg += ")"
        print(msg); sys.stdout.flush()
