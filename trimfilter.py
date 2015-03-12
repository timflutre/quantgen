#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: trim and filter reads in fastq files
# Copyright (C) 2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

# Inspired from:
# http://bcb.io/2009/08/09/trimming-adaptors-from-short-read-sequences/ by Brad Chapman
# http://news.open-bio.org/news/2009/09/biopython-fast-fastq/ by Peter Cock

# TODO:
# filter reads on quality
# implement Trimmomatic's palindrom mode for adapter read-through

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
import re
import itertools

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Restriction
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio import pairwise2
from Bio import Align
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "0.2.0" # http://semver.org/


class Trimfilter(object):
    
    def __init__(self):
        self.verbose = 1
        self.inDir = "."
        self.inFqFile1 = ""
        self.inFqFile2 = ""
        self.outPrefix = ""
        self.adpFile = ""
        self.maxNbErrors = 0 # i.e. exact matches
        self.dAdps = {}   # adapters to be searched in both reads of a pair
        self.dAdps1 = {}  # adapters to be searched in the R1 read of a pair
        self.dAdps2 = {}  # adapters to be searched in the R2 read of a pair
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' trims and filters reads in fastq files.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "      --idir\tpath to the input directory with the fastq files (default=.)\n"
        msg += "      --ifq1\tpath to the first input fastq file\n"
        msg += "\t\tcan be compressed with gzip\n"
        msg += "      --ifq2\tpath to the second input fastq file, optional\n"
        msg += "\t\tcan be compressed with gzip\n"
        msg += "\t\terror raised if reads not in same order as --ifq1\n"
        msg += "      --adp\tpath to file with adapters in fasta format\n"
        msg += "\t\tsame format as Trimmomatic, but only 'simple' mode is implemented\n"
        msg += "      --err\tmax number of errors (default=0, i.e. exact matches)\n"
        msg += "      --op\tprefix for the output files (2 paired, 2 unpaired, 1 log)\n"
        msg += "\t\twill be compressed with gzip\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --ifq1 reads1.fastq.gz --ifq2 reads2.fastq.gz --adp adapters.fa --op test\n" % os.path.basename(sys.argv[0])
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
        msg += "Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA).\n"
        msg += "License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt(sys.argv[1:], "hVv:",
                                       ["help", "version", "verbose=",
                                        "idir=", "ifq1=", "ifq2=", "adp=",
                                        "err=", "op=",])
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
            elif o == "--idir":
                self.inDir = a
            elif o == "--ifq1":
                 self.inFqFile1 = a
            elif o == "--ifq2":
                 self.inFqFile2 = a
            elif o == "--adp":
                self.adpFile = a
            elif o == "--err":
                self.maxNbErrors = int(a)
            elif o == "--op":
                 self.outPrefix = a
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if self.inDir == "":
            msg = "ERROR: missing compulsory option --idir"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.inDir):
            msg = "ERROR: can't find dir %s" % self.inDir
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.inFqFile1 == "":
            msg = "ERROR: missing compulsory option --ifq1"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        self.inFqFile1 = "%s/%s" % (self.inDir, self.inFqFile1)
        if not os.path.exists(self.inFqFile1):
            msg = "ERROR: can't find file %s" % self.inFqFile1
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.inFqFile2 != "":
            self.inFqFile2 = "%s/%s" % (self.inDir, self.inFqFile2)
            if not os.path.exists(self.inFqFile2):
                msg = "ERROR: can't find file %s" % self.inFqFile2
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
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
        if self.maxNbErrors < 0:
            self.maxNbErrors = 0
        if self.outPrefix == "":
            msg = "ERROR: missing compulsory option --op"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
            
    def loadAdapters(self):
        if self.verbose > 0:
            msg = "load adapter file..."
            print(msg); sys.stdout.flush()
            
        tmp = SeqIO.to_dict(SeqIO.parse(self.adpFile, "fasta",
                                        alphabet=IUPAC.unambiguous_dna))
        for adp in tmp:
            if "/1" in tmp[adp].id:
                self.dAdps1[tmp[adp].id] = str(tmp[adp].seq)
            elif "/2" in tmp[adp].id:
                self.dAdps2[tmp[adp].id] = str(tmp[adp].seq)
            else:
                self.dAdps[tmp[adp].id] = str(tmp[adp].seq)
                
        if self.inFqFile2 == "" and len(self.dAdps) == 0:
            msg = "ERROR: no adapter was loaded"
            sys.stderr.write("%s\n\n" % msg)
            sys.exit(1)
            
        if self.verbose > 0:
            if self.inFqFile2 == "":
                msg = "nb of adapters: %i" % (len(self.dAdps))
            else:
                msg = "nb of adapters: R1=%i R2=%i both=%i" % (len(self.dAdps1),
                                                               len(self.dAdps2),
                                                               len(self.dAdps))
            print(msg); sys.stdout.flush()
            
            
    def openOutputFiles(self):
        dOutHandles = {}
        dOutHandles["log"] = gzip.open("%s.log.gz" % self.outPrefix, "w")
        if self.inFqFile2 == "":
            dOutHandles["out"] = gzip.open("%s.fq.gz" % self.outPrefix, "w")
        else:
            dOutHandles["P1"] = gzip.open("%s_%s_R1.fq.gz" %
                                          (self.outPrefix, "paired"), "w")
            dOutHandles["U1"] = gzip.open("%s_%s_R1.fq.gz" %
                                          (self.outPrefix, "unpaired"), "w")
            dOutHandles["P2"] = gzip.open("%s_%s_R2.fq.gz" %
                                          (self.outPrefix, "paired"), "w")
            dOutHandles["U1"] = gzip.open("%s_%s_R2.fq.gz" %
                                          (self.outPrefix, "unpaired"), "w")
        return dOutHandles
        
        
    def getStartExactMatches(self, seq, adpSeq):
        """
        >>> i = Trimfilter()
        >>> i.getStartExactMatches("ATGCTTT", "TTT")
        4
        >>> i.getStartExactMatches("ATGCTT", "TTT")
        -1
        """
        return seq.find(adpSeq)
        
        
    def getStartFuzzyMatches(self, seq, adpSeq):
        """
        >>> i = Trimfilter()
        >>> i.getStartFuzzyMatches("CCCCCCCCCCTT", "TTT")
        10
        >>> i.getStartFuzzyMatches("CCCCCCCCCCCTATT", "CTATTT")
        10
        >>> i.getStartFuzzyMatches("CCCCCCCCCCTGT", "TTT")
        10
        """
        # aligns = pairwise2.align.localms(seq,
        #                                  adpSeq,
        #                                  5.0, -4.0, -9.0, -0.5,
        #                                  # one_alignment_only=True,
        #                                  gap_char="-")
        # readAlignedRegion = aligns[0][0][aligns[0][3]:aligns[0][4]]
        # motifAlignedRegion = aligns[0][1][aligns[0][3]:aligns[0][4]]
        # nbMatches = sum((1 if s == motifAlignedRegion[i] else 0) 
        #                 for i, s in enumerate(readAlignedRegion))
        seq_a, adaptor_a, score, start, end = pairwise2.align.localms(
            seq, adpSeq,
            5.0, -4.0, -9.0, -0.5,
            one_alignment_only=True,
            gap_char="-")[0]
        return start
        
        
    def getAdpStarts(self, seq, dAdps, dAdp2Start):
        """
        Search for the starts of exact or fuzzy matches of adapters on the sequence.
        
        seq -- sequence to search for
        dAdps -- dict with adapter names as keys and adapters sequence as values
        dAdp2Start -- dict with adapter names as keys and match starts as values
        
        >>> i = Trimfilter()
        >>> i.maxNbErrors = 0
        >>> dAdps = {"adp1": "TTT"}
        >>> dAdp2Start = {}
        >>> i.getAdpStarts("ATGCTTT", dAdps, dAdp2Start)
        {u'adp1': 4}
        """
        if(len(dAdps) > 0):
            for (adpName,adpSeq) in dAdps.items():
                dAdp2Start[adpName] = self.getStartExactMatches(seq, adpSeq)
                if self.maxNbErrors > 0 and dAdp2Start[adpName] == -1:
                    dAdp2Start[adpName] = self.getStartFuzzyMatches(seq, adpSeq)
        return dAdp2Start
        
        
    def logAdpMatches(self, handle, read_id, dAdp2Start):
        for (adpName,adpStart) in dAdp2Start.items():
            handle.write("%s\t%s\t%i\n" % (read_id, adpName, adpStart))
            
            
    def getTrimmingCoord(self, seqLen, dAdp2Start):
        """
        >>> i = Trimfilter()
        >>> i.getTrimmingCoord(13, {"adp1": -1})
        13
        >>> i.getTrimmingCoord(13, {"adp2": 10})
        10
        >>> i.getTrimmingCoord(13, {"adp1": -1, "adp2": 10})
        10
        """
        idx = seqLen
        for (adpName,adpStart) in dAdp2Start.items():
            if adpStart != -1 and adpStart < idx:
                idx = adpStart
        return idx
        
        
    def trimfilterSingleReads(self):
        """
        Read the fastq file, trim and filter, write the output.
        """
        if self.verbose > 0:
            msg = "trim and filter single-end reads (err=%i) ..." % self.maxNbErrors
            print(msg); sys.stdout.flush()
            
        if self.inFqFile1.endswith(".gz"):
            inFqHandle = gzip.open(self.inFqFile1, "r")
        else:
            inFqHandle = open(self.inFqFile1, "r")
            
        reads = FastqGeneralIterator(inFqHandle)
        
        dOutHandles = self.openOutputFiles()
        dOutHandles["log"].write("read.id\tadp.name\tadp.start\n")
        nbReads = 0
        nbTrimmedReads = 0
        nbFilteredReads = 0
        
        for (read_id, read_seq, read_q) in reads:
            nbReads += 1
            
            dAdp2Start = {}
            dAdp2Start = self.getAdpStarts(read_seq, self.dAdps, dAdp2Start)
            
            self.logAdpMatches(dOutHandles["log"], read_id, dAdp2Start)
            idx = self.getTrimmingCoord(len(read_seq), dAdp2Start)
            dOutHandles["out"].write("@%s\n%s\n+\n%s\n" %
                                     (read_id,
                                      read_seq[:idx],
                                      read_q[:idx]))
            
            if idx != len(read_seq):
                nbTrimmedReads += 1
                
        inFqHandle.close()
        [handle.close() for (out,handle) in dOutHandles.items()]
        
        if self.verbose > 0:
            msg = "total nb of reads: %i" % nbReads
            msg += "\nnb of trimmed reads: %i" % nbTrimmedReads
            print(msg); sys.stdout.flush()
            
            
    def trimfilterPairedReads(self):
        """
        Read the fastq files, trim and filter, write the output.
        """
        if self.verbose > 0:
            msg = "trim and filter paired-end reads (err=%i) ..." % self.maxNbErrors
            print(msg); sys.stdout.flush()
            
        if self.inFqFile1.endswith(".gz"):
            inFqHandle1 = gzip.open(self.inFqFile1, "r")
        else:
            inFqHandle1 = open(self.inFqFile1, "r")
        if self.inFqFile2.endswith(".gz"):
            inFqHandle2 = gzip.open(self.inFqFile2, "r")
        else:
            inFqHandle2 = open(self.inFqFile2, "r")
            
        reads1 = FastqGeneralIterator(inFqHandle1)
        reads2 = FastqGeneralIterator(inFqHandle2)
        
        dOutHandles = self.openOutputFiles()
        dOutHandles["log"].write("read.id\tadp.name\tadp.start\n")
        nbPairs = 0
        nbTrimmedPairs = 0
        nbTrimmedRead1 = 0
        nbTrimmedRead2 = 0
        nbFilteredPairs = 0
        
        for (read1_id, read1_seq, read1_q), (read2_id, read2_seq, read2_q) \
            in itertools.izip(reads1, reads2):
            if read1_id.split(" ")[0] != read2_id.split(" ")[0]:
                msg = "ERROR: for pair %i, reads %s and %s are not paired" \
                      % (nbPairs, read1_id, read2_id)
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            nbPairs += 1
            
            dAdp2Start1 = {}
            dAdp2Start2 = {}
            dAdp2Start1 = self.getAdpStarts(read1_seq, self.dAdps1, dAdp2Start1)
            dAdp2Start2 = self.getAdpStarts(read2_seq, self.dAdps2, dAdp2Start2)
            dAdp2Start1 = self.getAdpStarts(read1_seq, self.dAdps, dAdp2Start1)
            dAdp2Start2 = self.getAdpStarts(read2_seq, self.dAdps, dAdp2Start2)
            
            self.logAdpMatches(dOutHandles["log"], read1_id, dAdp2Start1)
            idx1 = self.getTrimmingCoord(len(read1_seq), dAdp2Start1)
            dOutHandles["P1"].write("@%s\n%s\n+\n%s\n" %
                                    (read1_id,
                                     read1_seq[:idx1],
                                     read1_q[:idx1]))
            self.logAdpMatches(dOutHandles["log"], read2_id, dAdp2Start2)
            idx2 = self.getTrimmingCoord(len(read2_seq), dAdp2Start2)
            dOutHandles["P2"].write("@%s\n%s\n+\n%s\n" %
                                    (read2_id,
                                     read2_seq[:idx2],
                                     read2_q[:idx2]))
            if idx1 != len(read1_seq) and idx2 != len(read2_seq):
                nbTrimmedPairs += 1
            elif idx1 != len(read1_seq) and idx2 == len(read2_seq):
                nbTrimmedRead1 += 1
            elif idx1 == len(read1_seq) and idx2 != len(read2_seq):
                nbTrimmedRead2 += 1
                
        inFqHandle1.close()
        inFqHandle2.close()
        [handle.close() for (out,handle) in dOutHandles.items()]
        
        if self.verbose > 0:
            msg = "total nb of read pairs: %i" % nbPairs
            msg += "\nnb of trimmed read pairs: %i" % nbTrimmedPairs
            print(msg); sys.stdout.flush()
            
            
    def run(self):
        self.loadAdapters()
        if self.inFqFile2 == "":
            self.trimfilterSingleReads()
        else:
            self.trimfilterPairedReads()
            
            
if __name__ == "__main__":
    i = Trimfilter()
    
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
