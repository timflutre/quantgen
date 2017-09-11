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
        
progVersion = "0.4.1" # http://semver.org/


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
        msg += "  ../gbs.py --proj1 testGbs --step 1 --samples samples.txt --pird ~/src/...\n"
        msg += "  ../gbs.py --proj1 testGbs --step 2 --samples samples.txt --pird ~/src/...\n"
        msg += "  ../gbs.py --proj1 testGbs --step 3 --samples samples.txt --pird ~/src/... --adp adapters.txt\n"
        msg += "  ../gbs.py --proj1 testGbs --proj2 testGbs-Atha --step 4 --samples samples_Atha_v2.txt --ref ~/src/...refgenome_Atha_v2 --dict ~/src/...refgenome_Atha_v2.dict\n"
        msg += "  ../gbs.py --proj2 testGbs-Atha --step 5 --samples samples_Atha_v2.txt --ref ~/src/...refgenome_Atha_v2\n"
        msg += "  ../gbs.py --proj2 testGbs-Atha --step 6 --samples samples_Atha_v2.txt --ref ~/src/...refgenome_Atha_v2\n"
        msg += "  ../gbs.py --proj2 testGbs-Atha --step 7 --samples samples_Atha_v2.txt --ref ~/src/...refgenome_Atha_v2\n"
        msg += "  ../gbs.py --proj2 testGbs-Atha --step 8 --samples samples_Atha_v2.txt --ref ~/src/...refgenome_Atha_v2 --jgid cohort-Atha\n"
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
        args = ["which", "bwa"]
        try:
            p = subprocess.check_output(args, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, e:
            msg = "can't find 'bwa' in PATH"
            raise ValueError(msg)
        args = ["which", "samtools"]
        try:
            p = subprocess.check_output(args, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, e:
            msg = "can't find 'samtools' in PATH"
            raise ValueError(msg)
        
        
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
        refgenome 1, 4 genotypes: geno1 (mother), geno2 (father), geno3 (child), geno4 (child)
        refgenome 2, 1 genotype: geno5
        2 lanes, 6 samples (parents present in both lanes)
        2 extractions, and thus 2 libraries, for geno1
        """
        samplesFile = "samples.txt"
        samplesAthaFile = "samples_Atha_v2.txt"
        samplesVvinFile = "samples_Vvin_v1.txt"
        lane2samples = {}
        refgenomeId2species = {"Atha_v2": "Arabidopsis thaliana",
                               "Vvin_v1": "Vitis vinifera"}
        
        samplesHandle = open(samplesFile, "w")
        samplesAthaHandle = open(samplesAthaFile, "w")
        samplesVvinHandle = open(samplesVvinFile, "w")
        
        lane2samples["1"] = {"R1": "init_reads_lane1_R1.fastq.gz",
                             "R2": "init_reads_lane1_R2.fastq.gz",
                             "genos": {"geno1": {"refgenomeId":"Atha_v2", "gen":0,
                                                 "lib":"geno1-A", "tag":"AAAA"},
                                       "geno2": {"refgenomeId":"Atha_v2", "gen":0,
                                                 "lib":"geno2", "tag":"GGGG"},
                                       "geno3": {"refgenomeId":"Atha_v2", "gen":1,
                                                 "lib":"geno3", "tag":"TTTT"}}}
        lane2samples["2"] = {"R1": "init_reads_lane2_R1.fastq.gz",
                             "R2": "init_reads_lane2_R2.fastq.gz",
                             "genos": {"geno1": {"refgenomeId":"Atha_v2", "gen":0,
                                                 "lib":"geno1-B", "tag":"CCCC"},
                                       "geno2": {"refgenomeId":"Atha_v2", "gen":0,
                                                 "lib":"geno2", "tag":"TTTT"},
                                       "geno4": {"refgenomeId":"Atha_v2", "gen":1,
                                                 "lib":"geno4", "tag":"AAAA"},
                                       "geno5": {"refgenomeId":"Vvin_v1", "gen":0,
                                                 "lib":"geno5", "tag":"AATTCC"}}}
        
        # header
        txt = "genotype"
        txt += "\tref_genome"
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
        samplesHandle.write("%s\n" % txt)
        samplesAthaHandle.write("%s\n" % txt)
        samplesVvinHandle.write("%s\n" % txt)
        
        # geno1 in lane 1
        txt = "geno1"
        txt += "\t%s" % lane2samples["1"]["genos"]["geno1"]["refgenomeId"]
        txt += "\t%s" % lane2samples["1"]["genos"]["geno1"]["lib"]
        txt += "\t%s" % lane2samples["1"]["genos"]["geno1"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 1
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["1"]["R1"]
        txt += "\t%s" % lane2samples["1"]["R2"]
        samplesHandle.write("%s\n" % txt)
        samplesAthaHandle.write("%s\n" % txt)
        
        # geno2 in lane 1
        txt = "geno2"
        txt += "\t%s" % lane2samples["1"]["genos"]["geno2"]["refgenomeId"]
        txt += "\t%s" % lane2samples["1"]["genos"]["geno2"]["lib"]
        txt += "\t%s" % lane2samples["1"]["genos"]["geno2"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%s" % 1
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["1"]["R1"]
        txt += "\t%s" % lane2samples["1"]["R2"]
        samplesHandle.write("%s\n" % txt)
        samplesAthaHandle.write("%s\n" % txt)
        
        # geno3 in lane 1
        txt = "geno3"
        txt += "\t%s" % lane2samples["1"]["genos"]["geno3"]["refgenomeId"]
        txt += "\t%s" % lane2samples["1"]["genos"]["geno3"]["lib"]
        txt += "\t%s" % lane2samples["1"]["genos"]["geno3"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%s" % 1
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["1"]["R1"]
        txt += "\t%s" % lane2samples["1"]["R2"]
        samplesHandle.write("%s\n" % txt)
        samplesAthaHandle.write("%s\n" % txt)
        
        # geno1 in lane 2
        txt = "geno1"
        txt += "\t%s" % lane2samples["2"]["genos"]["geno1"]["refgenomeId"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno1"]["lib"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno1"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 2
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["2"]["R1"]
        txt += "\t%s" % lane2samples["2"]["R2"]
        samplesHandle.write("%s\n" % txt)
        samplesAthaHandle.write("%s\n" % txt)
        
        # geno2 in lane 2
        txt = "geno2"
        txt += "\t%s" % lane2samples["2"]["genos"]["geno2"]["refgenomeId"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno2"]["lib"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno2"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 2
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["2"]["R1"]
        txt += "\t%s" % lane2samples["2"]["R2"]
        samplesHandle.write("%s\n" % txt)
        samplesAthaHandle.write("%s\n" % txt)
        
        # geno4 in lane 2
        txt = "geno4"
        txt += "\t%s" % lane2samples["2"]["genos"]["geno4"]["refgenomeId"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno4"]["lib"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno4"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 2
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["2"]["R1"]
        txt += "\t%s" % lane2samples["2"]["R2"]
        samplesHandle.write("%s\n" % txt)
        samplesAthaHandle.write("%s\n" % txt)
        
        # geno5 in lane 2
        txt = "geno5"
        txt += "\t%s" % lane2samples["2"]["genos"]["geno5"]["refgenomeId"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno5"]["lib"]
        txt += "\t%s" % lane2samples["2"]["genos"]["geno5"]["tag"]
        txt += "\tGenoToul"
        txt += "\tILLUMINA"
        txt += "\tHiSeq 2000"
        txt += "\tF6YMDACDT"
        txt += "\t%i" % 2
        txt += "\t2015-01-15"
        txt += "\t%s" % lane2samples["2"]["R1"]
        txt += "\t%s" % lane2samples["2"]["R2"]
        samplesHandle.write("%s\n" % txt)
        samplesVvinHandle.write("%s\n" % txt)
        
        samplesHandle.close()
        samplesAthaHandle.close()
        samplesVvinHandle.close()
        
        return samplesFile, lane2samples, refgenomeId2species
    
    
    def simulRefGenomeWoMotif(self, nbChrs, lenChr, speciesName,
                              refgenomeVersion):
        """
        >>> i = TestGbs()
        >>> obs = i.simulRefGenomeWoMotif(2, 100, "Arabidopsis thaliana", "2")
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
                                    description="%s|v%s|chr%i" % \
                                    (speciesName, refgenomeVersion,
                                     chrNum + 1)))
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
        
        
    def simulRefGenomes(self, refgenomeId2species, nbChrs, lenChr,
                        nbFragsPerChr, minFragLen, maxFragLen, motif):
        """
        the reference genome is haploid
        motif -- string with what remains of the restriction site after being cut
        """
        dRefGens = {}
        dFragCoords = {}
        
        for refgenomeId in refgenomeId2species:
            speciesName = refgenomeId2species[refgenomeId]
            refgenomeVersion = refgenomeId.split("_")[1].replace("v", "")
            refGen = self.simulRefGenomeWoMotif(nbChrs, lenChr, speciesName,
                                                refgenomeVersion)
            lFragCoords = self.chooseFragCoords(nbChrs, lenChr, nbFragsPerChr,
                                                minFragLen, maxFragLen)
            refGen = self.insertMotifsPerFrag(refGen, lFragCoords, motif)
            dRefGens[refgenomeId] = refGen
            dFragCoords[refgenomeId] = lFragCoords
            
        return dRefGens, dFragCoords
    
    
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
    
    
    def simulGenoGenomes(self, lane2samples, dRefGens):
        geno2genome = {}
        
        parents = []
        for lane in lane2samples:
            for geno in lane2samples[lane]["genos"]:
                if geno not in geno2genome and \
                   lane2samples[lane]["genos"][geno]["gen"] == 0:
                    geno2genome[geno] = self.simulParentGenome(
                        dRefGens[lane2samples[lane]["genos"][geno]["refgenomeId"]])
                    parents.append(geno)
                    
        for lane in lane2samples:
            for geno in lane2samples[lane]["genos"]:
                if geno not in geno2genome and \
                   lane2samples[lane]["genos"][geno]["gen"] == 1:
                    geno2genome[geno] \
                        = self.simulChildGenome(geno2genome[parents[0]],
                                                geno2genome[parents[1]])
                    
        return geno2genome
    
    
    def simulInserts(self, lane2samples, dFragCoords):
        """
        lane2samples -- dict of dicts ...
        dFragCoords -- dict of lists
        TODO: some fragments could not be sequenced, in all or some genotypes
        """
        lane2inserts = {}
        for lane in lane2samples:
            lane2inserts[lane] = {}
            for geno in lane2samples[lane]["genos"]:
                lane2inserts[lane][geno] = dFragCoords[
                    lane2samples[lane]["genos"][geno]["refgenomeId"]]
        return lane2inserts
    
    
    def extractReadSequences(self, genoGenChr, lInsertsChr, barcode, lenRead):
        """
        >>> i = TestGbs()
        >>> genoGenChr = SeqRecord(Seq("AATTTAGGGA"), id="chr1")
        >>> i.extractReadSequences(genoGenChr, [[3,9]], "C", 3)
        [['CTT', 'CCC']]
        """
        return [[str(barcode) \
                 + str(genoGenChr.seq[(i[0]-1):(i[0]-1+lenRead-len(barcode))]),
                 reverse_complement(str(genoGenChr.seq[(i[1]-lenRead):i[1]]))]
                for i in lInsertsChr]
    
    
    def simulPairedReadsPerSample(self, lane, barcode, genoGen, lInserts,
                                  nbReads, lenRead, readCounter):
        """
        barcode -- string
        genoGen -- [SeqRecChr1, SeqRecChr2, ...]
        lInserts -- [listChr1[], listChr2[], ...]
        TODO: use nbReads
        """
        lR1 = []
        lR2 = []
        nbChrs = len(genoGen)
        for c in range(nbChrs):
            lInsertPairs = self.extractReadSequences(genoGen[c], lInserts[c],
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
    
    
    def simulPairedReadsPerLane(self, lane, samples, geno2gen, lInserts, nbReads,
                                lenRead):
        """
        lane -- string
        samples -- dict[geno]{gen:..., tag:...}
        geno2gen -- dict[geno][SeqRecChr1, SeqRecChr2, ...]
        lInserts -- dict{geno}[ ... ]
        """
        lR1 = []
        lR2 = []
        readCounter = 0
        for geno in samples:
            lSampleR1, lSampleR2, readCounter \
                = self.simulPairedReadsPerSample(lane,
                                                 samples[geno]["tag"],
                                                 geno2gen[geno],
                                                 lInserts[geno],
                                                 nbReads,
                                                 lenRead,
                                                 readCounter)
            lR1 = lR1 + lSampleR1
            lR2 = lR2 + lSampleR2
        return lR1, lR2
    
    
    def simulReads(self, lane2samples, geno2gen, lane2inserts, nbReads, lenRead):
        """
        lane2samples -- dict[lane]{R1:file, R2:file, samples={geno:{gen, tag}, ...}}
        geno2gen -- dict[geno][SeqRecChr1, SeqRecChr2, ...]
        lane2inserts: dict{lane}{geno}[listChr1[], listChr2[], ...]
        """
        lane2reads = {}
        for lane in lane2samples:
            lane2reads[lane] = {"R1": [], "R2": []}
            lane2reads[lane]["R1"], lane2reads[lane]["R2"] \
                = self.simulPairedReadsPerLane(lane,
                                               lane2samples[lane]["genos"],
                                               geno2gen,
                                               lane2inserts[lane],
                                               nbReads,
                                               lenRead)
        return lane2reads
    
    
    def saveRefGenomes(self, dRefGens):
        dRefFiles = {}
        for refgenomeId in dRefGens:
            refFile = "refgenome_%s.fa" % refgenomeId
            refHandle = open(refFile, "w")
            for chrRec in dRefGens[refgenomeId]:
                refHandle.write(chrRec.format("fasta"))
            refHandle.close()
            dRefFiles[refgenomeId] = refFile
        return dRefFiles
    
    
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
                
                
    def makeInputSequences(self, refgenomeId2species, lane2samples,
                           nbChrs=2, lenChr=100000, nbFragsPerChr=50,
                           minFragLen=30, maxFragLen=600, motif="CAGC",
                           nbReads=1000, lenRead=100, paired=True):
        dRefGens, dFragCoords \
            = self.simulRefGenomes(refgenomeId2species, nbChrs, lenChr,
                                   nbFragsPerChr, minFragLen,
                                   maxFragLen, motif)
        
        geno2gen = self.simulGenoGenomes(lane2samples, dRefGens)
        
        lane2inserts = self.simulInserts(lane2samples, dFragCoords)
        
        lane2reads = self.simulReads(lane2samples, geno2gen, lane2inserts, nbReads, lenRead)
        
        dRefFiles = self.saveRefGenomes(dRefGens)
        self.saveReads(lane2samples, lane2reads, paired)
        return dRefGens, dRefFiles
    
    
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
    
    
    def makeDictFiles(self, refgenomeId2species, dRefGens, dRefFiles):
        dDictFiles = {}
        
        for refgenomeId in refgenomeId2species:
            dictFile = "refgenome_%s.dict" % refgenomeId
            
            dictHandle = open(dictFile, "w")
            dictHandle.write("@HD\tVN:1.4\tSO:unsorted\n")
            refGen = dRefGens[refgenomeId]
            for c in range(len(refGen)):
                txt = "@SQ"
                txt += "\tSN:%s" % refGen[c].id
                txt += "\tLN:%i" % len(refGen[c])
                txt += "\tSP:%s" % refgenomeId2species[refgenomeId]
                txt += "\tAS:%s" % refgenomeId.split("_")[1].replace("v", "")
                txt += "\tM5:%s" % hashlib.md5(str(refGen[c].seq)).hexdigest()
                txt += "\tUR:file:%s" % dRefFiles[refgenomeId]
                dictHandle.write("%s\n" % txt)
            dictHandle.close()
            dDictFiles[refgenomeId] = dictFile
            
        return dDictFiles
    
    
    def makeBwaIndexFiles(self, dRefFiles):
        for refgenomeId in dRefFiles:
            refFile = dRefFiles[refgenomeId]
            args = ["bwa", "index", "-p", refFile.split(".fa")[0], refFile]
            stdoutFilehandle = open("bwa_index.log", "w")
            subprocess.check_call(args, stdout=stdoutFilehandle,
                                  stderr=subprocess.STDOUT)
            stdoutFilehandle.close()
            
            
    def makeFaIndexFiles(self, dRefFiles):
        for refgenomeId in dRefFiles:
            refFile = dRefFiles[refgenomeId]
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
        
        options.append("--proj1")
        options.append("testGbs")
        options.append("--queue")
        options.append(self.queue)
        if self.lResources:
            options.append("--resou")
            options.append(" ".join(self.lResources))
        options.append("--step")
        options.append("1")
        
        samplesFile, lane2samples, refgenomeId2species = self.makeSamplesFile()
        nbChrs = 2
        lenChr = 100000
        nbFragsPerChr = 50
        minFragLen = 30
        maxFragLen = 300
        motif = "CAGC"
        nbReads = 1000 # for a given sample, over all its chrs
        lenRead = 100
        dRefGens, dRefFiles = self.makeInputSequences(
            refgenomeId2species, lane2samples,
            nbChrs, lenChr, nbFragsPerChr, minFragLen, maxFragLen, motif,
            nbReads, lenRead)
        adpFile, dAdps = self.makeAdapterFile()
        dDictFiles = self.makeDictFiles(refgenomeId2species, dRefGens,
                                        dRefFiles)
        self.makeBwaIndexFiles(dRefFiles)
        self.makeFaIndexFiles(dRefFiles)
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
