#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: demultiplex samples in fastq files
# Copyright (C) 2014-2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut], Laurène Gay [ctb], Nicolas Rode [ctb]
# Versioning: https://github.com/timflutre/quantgen

# Inspired from:
# https://gist.github.com/seandavi/3015625 by Sean Devi
# http://bcb.io/2009/08/09/trimming-adaptors-from-short-read-sequences/ by Brad Chapman
# http://news.open-bio.org/news/2009/09/biopython-fast-fastq/ by Peter Cock

# TODO:
# try fuzzy matching with https://pypi.python.org/pypi/regex
# try https://github.com/faircloth-lab/splitaake/blob/master/bin/splitaake_reads_many_gz.py
# try https://humgenprojects.lumc.nl/svn/fastools/trunk/fastools/demultiplex.py
# try https://github.com/pelinakan/UBD
# try http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage
# try http://www.dnabaser.com/download/nextgen-fastq-editor/index.html

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
        
progVersion = "1.9.0" # http://semver.org/


class Demultiplex(object):
    
    def __init__(self):
        self.verbose = 1
        self.inDir = "."
        self.inFqFile1 = ""
        self.inFqFile2 = ""
        self.tagFile = ""
        self.outPrefix = ""
        self.method = "4a"
        self.dist = 20
        self.restrictEnzyme = None
        self.findChimeras = "1"
        self.clipIdx = True
        self.tags = {} # keys are individuals and values are sequences (as string)
        self.cutMotif = ""
        self.regexpCutMotif = "" # used if findChimeras != 0
        self.regexpCompilCutMotif = "" # used if findChimeras != 0
        self.lenRemainMotif = -1
        self.regexpRemainMotif = "" # remaining motif (as uncompiled regexp)
        self.patterns = {} # keys are individuals and values are tag+motif (as compiled regexp)
        self.dInd2NbAssigned = {}
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' demultiplexes samples in fastq files.\n" % os.path.basename(sys.argv[0])
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
        msg += "      --ifq2\tpath to the second input fastq file\n"
        msg += "\t\tcan be compressed with gzip\n"
        msg += "\t\terror raised if reads not in same order as --ifq1\n"
        msg += "      --it\tpath to the tag file\n"
        msg += "\t\tonly A, T, G and C (no ambiguous nucleotide allowed)\n"
        msg += "\t\tcan be in 2 formats (automatically detected)\n"
        msg += "\t\t fasta: put sample names in the fasta headers\n"
        msg += "\t\t table: 2 columns, header line should be 'id\\ttag'\n"
        msg += "      --ofqp\tprefix for the output fastq files (2 per ind, 1 for unassigned, 1 for chimeras)\n"
        msg += "\t\tthe suffix will be based on \".fastq.gz\" (not \".fq\" because of FastQC)\n"
        msg += "\t\twill be compressed with gzip\n"
        msg += "\t\tthis prefix will also be used for a gzipped text file\n"
        msg += "\t\t counting assigned read pairs per individual\n"
        msg += "      --met\tmethod to assign pairs of reads to individuals via tags (default=4a)\n"
        msg += "\t\tonly forward strand is considered\n"
        msg += "\t\t1: assign pair if both reads start with the tag\n"
        msg += "\t\t2: assign pair if at least one read starts with the tag\n"
        msg += "\t\t3: same as 2 but count if one or both reads start with the tag\n"
        msg += "\t\t4a: assign pair if first read starts with tag (ignore second read)\n"
        msg += "\t\t4b: assign pair if first read has tag in its first N bases\n"
        msg += "\t\t    (ignore second read, see --dist)\n"
        msg += "\t\t4c: assign pair only if first read has tag and remaining cut site\n"
        msg += "\t\t    in its first N bases (ignore second read, see --dist and --re)\n"
        msg += "\t\t4d: assign pair only if first and/or second read has tag and\n"
        msg += "\t\t    remaining cut site in its first N bases (see --dist and --re)\n"
        msg += "\t\t    PCR chimeras (R1 tag is different than R2 tag) are detected\n"
        msg += "\t\t    and saved in distinct files than the unassigned\n"
        msg += "      --dist\tdistance from the read start to search for the tag (in bp, default=20)\n"
        msg += "      --re\tname of the restriction enzyme (e.g. 'ApeKI')\n"
        msg += "      --chim\tsearch if full restriction site found in R1 and/or R2 (default=\"1\", see --re)\n"
        msg += "\t\t0: don't search (but some chimeras may still be detected if --met 4d)\n"
        msg += "\t\t1: if chimera, count as such, try to assign, and save in same files as others\n"
        msg += "\t\t2: if chimera, don't even try to assign and save in distinct files\n"
        msg += "      --nci\tdo not clip the tag when saving the assigned reads\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --ifq1 reads1.fastq.gz --ifq2 reads2.fastq.gz --ifat tags.fa --ofqp test --met 3\n" % os.path.basename(sys.argv[0])
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
        msg += "Copyright (C) 2014-2015 Institut National de la Recherche Agronomique (INRA).\n"
        msg += "License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut], Laurène Gay [ctb], Nicolas Rode [ctb]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt(sys.argv[1:], "hVv:",
                                       ["help", "version", "verbose=",
                                        "idir=", "ifq1=", "ifq2=", "it=",
                                        "ofqp=", "met=", "re=", "chim=",
                                        "dist=", "nci"])
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
            elif o == "--it":
                 self.tagFile = a
            elif o == "--ofqp":
                 self.outPrefix = a
            elif o == "--met":
                self.method = a
            elif o == "--dist":
                self.dist = int(a)
            elif o == "--re":
                try:
                    self.restrictEnzyme = Restriction.__dict__[a]
                except KeyError:
                    msg = "ERROR: restriction enzyme %s not recognized" % a
                    sys.stderr.write("%s\n\n" % msg)
                    sys.exit(1)
            elif o == "--chim":
                self.findChimeras = a
            elif o == "--nci":
                self.clipIdx = False
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
        if self.inFqFile2 == "":
            msg = "ERROR: missing compulsory option --ifq2"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        self.inFqFile2 = "%s/%s" % (self.inDir, self.inFqFile2)
        if not os.path.exists(self.inFqFile2):
            msg = "ERROR: can't find file %s" % self.inFqFile2
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.tagFile == "":
            msg = "ERROR: missing compulsory option --it"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.tagFile):
            msg = "ERROR: can't find file %s" % self.tagFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.outPrefix == "":
            msg = "ERROR: missing compulsory option --ofqp"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.method not in ["1","2","3","4a","4b","4c","4d","chim"]:
            msg = "ERROR: unknown option --met %s" % self.method
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if (self.method in ["4c","4d"] or self.findChimeras != "0") \
           and self.restrictEnzyme == None:
            msg = "ERROR: missing compulsory option --re"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.findChimeras not in ["0","1","2"]:
            msg = "ERROR: --chim %s is unknown" % self.findChimeras
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.dist < 0:
            self.dist = 0
            
            
    def findTagFileFormat(self):
        """
        Look for a '>' in the tag file and decide to return a 'fasta' or 'table' file format.
        """
        tagFileFormat = "table"
        tagHandle = open(self.tagFile)
        line = tagHandle.readline()
        tagHandle.close()
        tokens = line.split()
        if line[0] == ">":
            tagFileFormat = "fasta"
        elif tokens == ["id", "tag"]:
            tagFileFormat = "table"
        else:
            if len(tokens) == 2:
                msg = "WARNING: tag file should have a header line 'id\\ttag'"
                sys.stderr.write("%s\n" % msg)
                tagFileFormat = "table"
            else:
                msg = "ERROR: tag file seem to be neither in 'fasta' nor 'table' format"
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
        return tagFileFormat
        
        
    def loadTags(self):
        tagFileFormat = self.findTagFileFormat()
        if self.verbose > 0:
            msg = "load tag file (format=%s)..." % tagFileFormat
            print(msg); sys.stdout.flush()
            
        if tagFileFormat == "fasta":
            tmp = SeqIO.to_dict(SeqIO.parse(self.tagFile, "fasta",
                                            alphabet=IUPAC.unambiguous_dna))
            for tag in tmp:
                self.tags[tmp[tag].id] = str(tmp[tag].seq)
        elif tagFileFormat == "table":
            tagHandle = open(self.tagFile)
            line = tagHandle.readline()
            tokens = line.split()
            if tokens == ["id", "tag"]:
                line = tagHandle.readline()
                tokens = line.split()
            while True:
                if line == "":
                    break
                if tokens[0] in self.tags:
                    msg = "ERROR: id '%s' is present several times" % tokens[0]
                    sys.stderr.write("%s\n" % msg)
                    sys.exit(1)
                self.tags[tokens[0]] = tokens[1]
                line = tagHandle.readline()
                tokens = line.split()
            tagHandle.close()
            
        for ind in self.tags:
            self.dInd2NbAssigned[ind] = 0
            
        if self.verbose > 0:
            msg = "nb of tags: %i" % len(self.tags)
            print(msg); sys.stdout.flush()
            
            
    def retrieveRestrictionEnzyme(self):
        if self.verbose > 0:
            print("retrieve restriction enzyme from Biopython...")
            
        self.cutMotif = self.restrictEnzyme.elucidate()
        if self.verbose > 0:
            print("enzyme %s: motif=%s" % (self.restrictEnzyme, self.cutMotif))
            
        if self.findChimeras != "0":
            for nt in list(self.cutMotif):
                if nt in ["A", "T", "G", "C"]:
                    self.regexpCutMotif += nt
                elif nt in ["^", "_"]:
                    continue
                else: # ambiguous letter from IUPAC code, e.g. W
                    self.regexpCutMotif += "[" + ambiguous_dna_values[nt] + "]" # eg. [AT] for W
            if self.verbose > 0:
                print("regexp of full motif: %s" % self.regexpCutMotif)

            self.regexpCompilCutMotif = re.compile(self.regexpCutMotif)
            
            
    def prepareRemainingMotif(self):
        if self.verbose > 0:
            print("prepare remaining motif..."); sys.stdout.flush()
            
        coordCutSense = self.cutMotif.find("^")
        remainMotifAmbig = self.restrictEnzyme.site[coordCutSense:len(self.restrictEnzyme.site)]
        self.lenRemainMotif = len(remainMotifAmbig)
        for nt in list(remainMotifAmbig):
            if nt in ["A", "T", "G", "C"]:
                self.regexpRemainMotif += nt
            else: # ambiguous letter from IUPAC code, e.g. W
                self.regexpRemainMotif += "[" + ambiguous_dna_values[nt] + "]" # eg. [AT] for W
                
        if self.verbose > 0:
            print("regexp of remaining motif: %s" % self.regexpRemainMotif)
            
            
    def compilePatterns(self):
        """
        Compile patterns (each tag + remaining right of cut site as regexp) 
        for quicker searches.
        """
        if self.verbose > 0:
            print("compile patterns..."); sys.stdout.flush()
        for tagId in self.tags:
            tag = self.tags[tagId]
            self.patterns[tagId] = re.compile(tag + self.regexpRemainMotif)
            
            
    def checkDist(self):
        """
        Check that length of tag (+ remaining motif) <= self.dist.
        """
        for tagId in self.tags:
            tag = self.tags[tagId]
            if self.dist > 0 and len(tag) + self.lenRemainMotif > self.dist:
                msg = "ERROR: --dist %i is too short for tag %s" % (self.dist,
                                                                    tagId)
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)


    def identifyChimeras(self, read1_seq, read2_seq):
        chimera = False
        if self.regexpCompilCutMotif.search(read1_seq) != None \
           or self.regexpCompilCutMotif.search(read2_seq) != None:
            chimera = True
        return chimera
        
        
    def identifyIndividual_1(self, read1_seq, read2_seq):
        """
        Assign pair if both reads start with the tag.
        """
        assigned = False
        ind = None
        idx = 0
        for tagId in self.tags:
            tag = self.tags[tagId]
            if read1_seq.startswith(tag) \
               and read2_seq.startswith(tag):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx = len(tag)
                break
        return assigned, tagId, idx, idx
        
        
    def identifyIndividual_2(self, read1_seq, read2_seq):
        """
        Assign pair if at least one read starts with the tag. !!!!! IdX1=0 or IdX2=0 if the tag is only found in the alternative read
        """
        assigned = False
        ind = None
        idx1 = 0
        idx2 = 0
        for tagId in self.tags:
            tag = self.tags[tagId]
            if read1_seq.startswith(tag):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = len(tag)
                    if read2_seq.startswith(tag):
                        idx2 = len(tag)
                break
            if read2_seq.startswith(tag):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx2 = len(tag)
                break
        return assigned, ind, idx1, idx2
        
        
    def identifyIndividual_3(self, read1_seq, read2_seq):
        """
        Same as 2 but count if one or both reads start with the tag.
        """
        assigned = False
        ind = None
        idx1 = 0
        idx2 = 0
        nbAssignedPairsTwoTags = 0
        nbAssignedPairsOneTag = 0
        for tagId in self.tags:
            tag = self.tags[tagId]
            if read1_seq.startswith(tag):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = len(tag)
                if read2_seq.startswith(tag):
                    if self.clipIdx:
                        idx2 = len(tag)
                    nbAssignedPairsTwoTags += 1
                else:
                    nbAssignedPairsOneTag += 1
                break
            if read2_seq.startswith(tag):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx2 = len(tag)
                nbAssignedPairsOneTag += 1
                break
        return assigned, ind, idx1, idx2, nbAssignedPairsTwoTags, \
            nbAssignedPairsOneTag
        
        
    def identifyIndividual_4a(self, read1_seq):
        """
        Assign pair if first read starts with the tag (ignore second read).
        """
        assigned = False
        ind = None
        idx1 = 0
        for tagId in self.tags:
            tag = self.tags[tagId]
            if read1_seq.startswith(tag):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = len(tag)
                break
        return assigned, tagId, idx1, 0
        
        
    def identifyIndividual_4b(self, read1_seq):
        """
        Assign pair if first read has tag in its first N bases (ignore second
        read, use self.dist).
        """
        assigned = False
        ind = None
        idx1 = 0
        for tagId in self.tags:
            tag = self.tags[tagId]
            tmpIdx = read1_seq[:self.dist].find(tag)
            if tmpIdx != -1:
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = tmpIdx + len(tag)
                break
        return assigned, tagId, idx1, 0
        
        
    def identifyIndividual_4c(self, read1_seq):
        """
        Assign pair if first read has tag and remaining cut site in its first
        N bases (ignore second read, use self.dist and self.patterns).
        """
        assigned = False
        ind = None
        idx1 = 0
        for tagId in self.patterns:
            tmpRe = self.patterns[tagId].search(read1_seq[:self.dist])
            if tmpRe != None:
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = tmpRe.end() - self.lenRemainMotif
                break
        return assigned, tagId, idx1, 0
        
        
    def identifyIndividual_4d(self, read1_seq, read2_seq):
        """
        Assign pair if first and/or second read has tag and remaining cut site in its first
        N bases (use self.dist and self.patterns). PCR chimeras are handled.
        """
        assigned = False
        ind = None
        idx1 = 0
        idx2 = 0
        AssignedPairsTwoTags = 0
        AssignedPairsOneTag = 0
        chimera = False
        
        lTagIds = self.patterns.keys()
        lTagIds.sort()
        for i in range(0, len(lTagIds)): # for each tag
            tagId = lTagIds[i]
            tmpRe1 = self.patterns[tagId].search(read1_seq[:self.dist])
            tmpRe2 = self.patterns[tagId].search(read2_seq[:self.dist])
            
            if tmpRe1 != None: # if first read has the pattern
                ind = tagId
                idx1 = tmpRe1.end() - self.lenRemainMotif
                if tmpRe2 != None: # if second read also has the pattern
                    assigned = True
                    idx2 = tmpRe2.end() - self.lenRemainMotif
                    AssignedPairsTwoTags += 1
                else: # if only the first read has the pattern
                    for tagId2 in self.patterns:
                        if tagId2 == tagId:
                            continue
                        tmpRe2chim = self.patterns[tagId2].search(read2_seq[:self.dist])
                        if tmpRe2chim != None: # if the first and second read have different tags 
                            chimera = True
                            break
                    if not chimera:
                        assigned = True
                        idx2 = idx1
                        AssignedPairsOneTag += 1
                break
                
            elif tmpRe2 != None: # if only the second read has the pattern
                ind = tagId
                idx2 = tmpRe2.end() - self.lenRemainMotif
                for j in range(i+1, len(lTagIds)):
                    tagId2 = lTagIds[j]
                    tmpRe1chim =  self.patterns[tagId2].search(read1_seq[:self.dist])
                    if tmpRe1chim != None: # if the first and second read have different tags
                        chimera = True
                        break
                if not chimera:
                    assigned = True
                    AssignedPairsOneTag += 1
                    idx1 = idx2
                break
            
        if not self.clipIdx:
            idx1 = idx2 = 0
            
        return assigned, tagId, idx1, idx2, AssignedPairsTwoTags, AssignedPairsOneTag, chimera
    
    
    def saveStatsPerInd(self):
        outFile = "%s_stats-demultiplex.txt.gz" % self.outPrefix
        outHandle = gzip.open(outFile, "w")
        
        txt = "ind"
        txt += "\tbarcode"
        txt += "\tassigned"
        outHandle.write("%s\n" % txt)
        
        lInds = self.tags.keys()
        lInds.sort()
        for ind in lInds:
            txt = "%s" % ind
            txt += "\t%s" % self.tags[ind]
            txt += "\t%i" % self.dInd2NbAssigned[ind]
            outHandle.write("%s\n" % txt)
            
        outHandle.close()
        
        
    def demultiplexPairedReads(self):
        """
        Read the data files, launch the identifying method, writes the output.
        """
        if self.verbose > 0:
            msg = "demultiplex paired-end reads (method=%s)..." % self.method
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
        
        dOutFqHandles = {}
        nbPairs = 0
        nbAssignedPairs = 0
        nbAssignedPairsTwoTags = 0
        nbAssignedPairsOneTag = 0
        nbChimeras = 0
        nbChimerasSite = 0
        nbChimerasTags = 0
        nbUnassignedPairs = 0
        nbUnassignedPairsChimeras = 0
        meanQuals = []
        
        # iterate over all read pairs
        for (read1_id, read1_seq, read1_q), (read2_id, read2_seq, read2_q) \
            in itertools.izip(reads1, reads2):
            
            if read1_id.split(" ")[0] != read2_id.split(" ")[0]:
                msg = "ERROR: for pair %i, reads %s and %s are not paired" \
                      % (nbPairs, read1_id, read2_id)
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            nbPairs += 1
            
            assigned = False
            ind = None
            chimera = False
            chimeraSite = False
            chimeraTags = False
            if self.findChimeras != "0":
                chimeraSite = self.identifyChimeras(read1_seq, read2_seq)
            if not chimeraSite or self.findChimeras == "1":
                if self.method == "1":
                    assigned, ind, idx1, idx2 = self.identifyIndividual_1(
                        read1_seq, read2_seq)
                elif self.method == "2":
                    assigned, ind, idx1, idx2 = self.identifyIndividual_2(
                        read1_seq, read2_seq)
                elif self.method == "3":
                    assigned, ind, idx1, idx2, t2, t1 = self.identifyIndividual_3(
                        read1_seq, read2_seq)
                elif self.method == "4a":
                    assigned, ind, idx1, idx2 = self.identifyIndividual_4a(read1_seq)
                elif self.method == "4b":
                    assigned, ind, idx1, idx2 = self.identifyIndividual_4b(read1_seq)
                elif self.method == "4c":
                    assigned, ind, idx1, idx2 = self.identifyIndividual_4c(read1_seq)
                elif self.method == "4d":
                    assigned, ind, idx1, idx2, t2, t1, chimeraTags = self.identifyIndividual_4d(
                        read1_seq, read2_seq)
                    
            if chimeraSite:
                nbChimerasSite += 1 
            if chimeraTags:
                nbChimerasTags += 1
            if chimeraSite or chimeraTags:
                chimera = True
                nbChimeras += 1
                
            if assigned:
                nbAssignedPairs += 1
                if self.method in ["3","4d"]:
                    nbAssignedPairsTwoTags += t2
                    nbAssignedPairsOneTag += t1
                if ind not in dOutFqHandles:
                    dOutFqHandles[ind] = [
                        gzip.open("%s_%s_R1.fastq.gz" %
                                  (self.outPrefix, ind), "w"),
                        gzip.open("%s_%s_R2.fastq.gz" %
                                  (self.outPrefix, ind), "w")]
                dOutFqHandles[ind][0].write("@%s\n%s\n+\n%s\n" %
                                            (read1_id,
                                             read1_seq[idx1:],
                                             read1_q[idx1:]))
                dOutFqHandles[ind][1].write("@%s\n%s\n+\n%s\n" %
                                            (read2_id,
                                             read2_seq[idx2:],
                                             read2_q[idx2:]))
                self.dInd2NbAssigned[ind] += 1
            else: # chimera or unassigned
                nbUnassignedPairs += 1
                if chimera and self.findChimeras != "1":
                    nbUnassignedPairsChimeras += 1
                    if "chimeras" not in dOutFqHandles:
                        dOutFqHandles["chimeras"] = [
                            gzip.open("%s_chimeras_R1.fastq.gz" %
                                      self.outPrefix, "w"),
                            gzip.open("%s_chimeras_R2.fastq.gz" %
                                      self.outPrefix, "w")]
                    dOutFqHandles["chimeras"][0].write("@%s\n%s\n+\n%s\n" %
                                                      (read1_id,
                                                       read1_seq,
                                                       read1_q))
                    dOutFqHandles["chimeras"][1].write("@%s\n%s\n+\n%s\n" %
                                                      (read2_id,
                                                       read2_seq,
                                                       read2_q))
                else: # chimera or unassigned
                    if "unassigned" not in dOutFqHandles:
                        dOutFqHandles["unassigned"] = [
                            gzip.open("%s_unassigned_R1.fastq.gz" %
                                      self.outPrefix, "w"),
                            gzip.open("%s_unassigned_R2.fastq.gz" %
                                      self.outPrefix, "w")]
                        
                    dOutFqHandles["unassigned"][0].write("@%s\n%s\n+\n%s\n" %
                                                         (read1_id,
                                                          read1_seq,
                                                          read1_q))
                    dOutFqHandles["unassigned"][1].write("@%s\n%s\n+\n%s\n" %
                                                         (read2_id,
                                                          read2_seq,
                                                          read2_q))
                    
        inFqHandle1.close()
        inFqHandle2.close()
        [handle.close() for (ind,handles) in dOutFqHandles.items() \
         for handle in handles]
        
        self.saveStatsPerInd()
        
        if self.verbose > 0:
            msg = "total nb of read pairs: %i" % nbPairs
            msg += "\nnb of assigned read pairs: %i" % nbAssignedPairs
            if self.method in ["3","4d"]:
                msg += "; 2tags=%i 1tags=%i" % (nbAssignedPairsTwoTags, \
                                                nbAssignedPairsOneTag)
            if self.findChimeras != "0" or self.method in ["4d"]:
                msg += "\nnb of chimeric read pairs: %i (%.2f%%" % (
                    nbChimeras,
                    100 * nbChimeras / float(nbPairs))
                if self.findChimeras != "0" and self.method in ["4d"]:
                    msg += "; site=%i tags=%i" % (nbChimerasSite,
                                                  nbChimerasTags)
                msg += ")"
            msg += "\nnb of unassigned read pairs: %i (%.2f%%" % (
                nbUnassignedPairs,
                100 * nbUnassignedPairs / float(nbPairs))
            msg += ")"
            if self.findChimeras == "2" or self.method in ["4d"]:
                msg += "\nnb of unassigned read pairs (excluding chimeras): %i (%.2f%%" % (
                    nbUnassignedPairs - nbUnassignedPairsChimeras,
                    100 * (nbUnassignedPairs - nbUnassignedPairsChimeras) / float(nbPairs))
                msg += ")"
            nbInds = len(dOutFqHandles)
            if "unassigned" in dOutFqHandles:
                nbInds -= 1
            if "chimeras" in dOutFqHandles:
                nbInds -= 1
            msg += "\nnb of individuals with assigned reads: %i" % nbInds
            print(msg); sys.stdout.flush()
            
            
    def run(self):
        self.loadTags()
        if self.method in ["4c","4d"] or self.findChimeras != "0":
            self.retrieveRestrictionEnzyme()
            self.prepareRemainingMotif()
            self.compilePatterns()
        self.checkDist()
        self.demultiplexPairedReads()
        
        
if __name__ == "__main__":
    i = Demultiplex()
    
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
