#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: demultiplex individuals in two paired fastq files
# Copyright (C) 2014 Institut National de la Recherche Agronomique
# License: GPL-3+
# Author: Timothée Flutre
# Versioning: https://github.com/timflutre/quantgen

# Inspired from:
# https://gist.github.com/seandavi/3015625 by Sean Devi
# https://bcbio.wordpress.com/2009/08/09/trimming-adaptors-from-short-read-sequences/ by Brad Chapman
# http://news.open-bio.org/news/2009/12/interleaving-paired-fastq-files-with-biopython/ by Peter Cock

# TODO:
# add method using regexp (http://stackoverflow.com/a/12390748/597069) with n >= 1 nucleotide(s) between the start and the tag
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

import itertools # for izip
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
        
progVersion = "1.0.0" # http://semver.org/


# class RestrictEnzyme(Restriction.RestrictionType):
#     def __init__(self, name):
        
#         SeqRecord.__init__(self,
#                            Seq(data.split("=")[1].replace("/","").upper(),
#                                IUPAC.ambiguous_dna),
#                            id=data.split("=")[0],
#                            description="restriction enzyme")
#         self.cut = data.split("=")[1].find("/") # 1 in case of ApeKI
#         self.remainr = SeqRecord(Seq(self.seq[self.cut:len(self.seq)],
#                                      IUPAC.ambiguous_dna),
#                                  id=self.id,
#                                  description="remaining seq right of the cut")
        
        
class Demultiplex(object):
    
    def __init__(self):
        self.verbose = 1
        self.inDir = "."
        self.inFqFile1 = ""
        self.inFqFile2 = ""
        self.tagFile = ""
        self.outFqPrefix = ""
        self.method = 2
        self.restrictEnzyme = None
        self.remainingMotifs = []
        self.dist = 25
        self.clipIdx = False
        self.tags = {}
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' demultiplexes individuals in two paired fastq files.\n" % os.path.basename(sys.argv[0])
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
        msg += "\t\tonly A, T, G and C (no ambiguous nucleotide)\n"
        msg += "\t\tcan be in 2 formats (automatically detected)\n"
        msg += "\t\t fasta: put sample names in the fasta headers\n"
        msg += "\t\t table: 2 columns, header line should be 'id\\ttag'\n"
        msg += "      --ofqp\tprefix for the output fastq files (2 per ind, 1 unassigned)\n"
        msg += "\t\twill be compressed with gzip\n"
        msg += "      --met\tmethod to assign paired reads to individuals (1/2/default=3/4)\n"
        msg += "\t\t1: assign if both reads start with the tag (only fwd)\n"
        msg += "\t\t2: assign if at least one read starts with the tag (only fwd)\n"
        msg += "\t\t3: same as 2 but count if one or both reads start with the tag (only fwd)\n"
        # msg += "\t\t4: assign if at least one read contains the tag next to the cut site (only fwd)\n"
        # msg += "      --dist\tdistance from the read start to look for cut site (in bp, default=25)\n"
        # msg += "\t\twith --met 4\n"
        # msg += "      --re\trestriction enzyme with its cut site (e.g. ApeKI=G/CWGC)\n"
        # msg += "\t\twith --met 4\n"
        msg += "      --ci\tclip the tag when saving the assigned reads\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --ifq1 reads1.fastq.gz --ifq2 reads2.fastq.gz --ifat tags.fa --ofqp test --ci\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA).\n"
        msg += "License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:",
                                        ["help", "version", "verbose=",
                                         "idir=", "ifq1=", "ifq2=", "it=",
                                         "ofqp=", "met=", "re=", "dist=",
                                         "ci"])
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
                 self.outFqPrefix = a
            elif o == "--met":
                self.method = a
            elif o == "--re":
                try:
                    self.restrictEnzyme = Restriction.__dict__[a]
                except KeyError:
                    msg = "ERROR: restriction enzyme %s not recognized" % a
                    sys.stderr.write("%s\n\n" % msg)
                    sys.exit(1)
            elif o == "--dist":
                self.dist = int(a)
            elif o == "--ci":
                self.clipIdx = True
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
        if self.outFqPrefix == "":
            msg = "ERROR: missing compulsory option --ofqp"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.method not in ["1","2","3","4"]:
            msg = "ERROR: wrong value for --met"
            sys.stderr.write("%s\n" % msg)
            self.help()
            sys.exit(1)
        if self.clipIdx:
            msg = "ERROR: --ci not yet implemented"
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)



    def findTagFileFormat(self):
        """
        Return 'fasta' or 'table'.
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
                msg = "ERROR: tag file doesn't seem to be in 'fasta' or 'table' format"
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
            
        if self.verbose > 0:
            msg = "nb of tags: %i" % len(self.tags)
            print(msg); sys.stdout.flush()
            
            
    def prepareRemainingMotifs(self):
        if self.verbose > 0:
            print("prepare remaining motifs..."); sys.stdout.flush()
            
        cutMotif = self.restrictEnzyme.elucidate() # "G^CWG_C" for ApeKI
        if self.verbose > 0:
            print("enzyme %s: motif=%s" % (self.restrictEnzyme, cutMotif))
        coordCutSense = cutMotif.find("^") # position of cut in sense strand, 1 for ApeKI
        remainMotifAmbig = self.restrictEnzyme.site[coordCutSense:len(self.restrictEnzyme.site)] # "CWGC" for ApeKI
        if self.verbose > 3: # debug
            print(remainMotifAmbig)
        for i,l in enumerate(remainMotifAmbig):
            if l in ["A", "T", "G", "C"]:
                if i == 0:
                    self.remainingMotifs.append(l)
                else:
                    for j in range(len(self.remainingMotifs)):
                        self.remainingMotifs[j] += l
            else: # ambiguous letter, e.g. W
                if i == 0:
                    for nt in ambiguous_dna_values[l]:
                        self.remainingMotifs.append(nt)
                else:
                    currNbRemainMotifs = len(self.remainingMotifs)
                    nbUnambig = len(ambiguous_dna_values[l])
                    for j in range(currNbRemainMotifs):
                        self.remainingMotifs \
                            = self.remainingMotifs + \
                            [self.remainingMotifs[j]] * (nbUnambig - 1)
                        self.remainingMotifs[j] += ambiguous_dna_values[l][0]
                        for k in range(1,nbUnambig):
                            self.remainingMotifs[len(self.remainingMotifs)-k] \
                                += ambiguous_dna_values[l][k]
                            
        if self.verbose > 0:
            print("nb of remaining motif(s): %i (%s)" \
                  % (len(self.remainingMotifs),
                     ", ".join(self.remainingMotifs)))
            
            
    def identifyIndividual_1(self, read1_seq, read2_seq):
        """
        Count if both read sequences of the pair starts with the tag (only fwd).
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
                idx = len(tag)
                break
        return assigned, tagId, idx, idx
        
        
    def identifyIndividual_2(self, read1_seq, read2_seq):
        """
        Count if at least one read sequence of the pair starts with the tag (only fwd).
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
                idx1 = len(tag)
                if read2_seq.startswith(tag):
                    idx2 = len(tag)
                break
            if read2_seq.startswith(tag):
                assigned = True
                ind = tagId
                idx2 = len(tag)
                break
        return assigned, ind, idx1, idx2
        
        
    def identifyIndividual_3(self, read1_seq, read2_seq):
        """
        Count if one or both reads start with the tag (only fwd).
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
                idx1 = len(tag)
                if read2_seq.startswith(tag):
                    idx2 = len(tag)
                    nbAssignedPairsTwoTags += 1
                else:
                    nbAssignedPairsOneTag += 1
                break
            if read2_seq.startswith(tag):
                assigned = True
                ind = tagId
                idx2 = len(tag)
                nbAssignedPairsOneTag += 1
                break
        return assigned, ind, idx1, idx2, nbAssignedPairsTwoTags, \
            nbAssignedPairsOneTag
        
        
    def identifyIndividual_4(self, read1, read2):
        """
        Count if at least one read contains the tag next to the cut site (only fwd).
        """
        assigned = False
        ind = None
        for tagId in self.tags:
            if read1.seq.startswith(self.tags[tagId]) \
               or read2.seq.startswith(self.tags[tagId]):
                assigned = True
                ind = tagId
                break
            elif self.tags[tagId] in read1.seq[0:self.dist]:
                coord = read1.seq[0:self.dist].find(self.tags[tagId]) # will be > 0
                for remainingMotif in self.remainingMotifs:
                    if read1.seq[(coord+len(self.tags[tagId])):(coord+len(self.tags[tagId])+len(remainingMotif))] \
                       == remainingMotif:
                        assigned = True
                        ind = tagId
                        print("ind %s, tag %s, read1 %s (%i)" \
                              % (ind, self.tags[tagId], read1.description, coord))
                        print(read1.seq[0:(self.dist+5)])
                        break
            elif self.tags[tagId] in read2.seq[0:self.dist]:
                coord = read2.seq[0:self.dist].find(self.tags[tagId]) # will be > 0
                for remainingMotif in self.remainingMotifs:
                    if read2.seq[(coord+len(self.tags[tagId])):(coord+len(self.tags[tagId])+len(remainingMotif))] \
                       == remainingMotif:
                        assigned = True
                        ind = tagId
                        break
        return assigned, ind
        
        
    # def findCoordCutSite(self, read):
    #     fuzzy matching
    #     aligns = pairwise2.align.localms(read.seq,
    #                                      self.RestrictEnzyme.motifFwd,
    #                                      5.0, -4.0, -9.0, -0.5,
    #                                      # one_alignment_only=True,
    #                                      gap_char="-")
    #     readAlignedRegion = aligns[0][0][aligns[0][3]:aligns[0][4]]
    #     motifAlignedRegion = aligns[0][1][aligns[0][3]:aligns[0][4]]
    #     nbMatches = sum((1 if s == motifAlignedRegion[i] else 0) 
    #                     for i, s in enumerate(readAlignedRegion))
        
        
    def demultiplexPairedReads(self):
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
        meanQuals = []
        
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
            if self.method == "1":
                assigned, ind, idx1, idx2 = self.identifyIndividual_1(
                    read1_seq, read2_seq)
            elif self.method == "2":
                assigned, ind, idx1, idx2 = self.identifyIndividual_2(
                    read1_seq, read2_seq)
            elif self.method == "3":
                assigned, ind, idx1, idx2, t2, t1 = self.identifyIndividual_3(
                    read1_seq, read2_seq)
            elif self.method == "4":
                assigned, ind = self.identifyIndividual_4(read1, read2)
                
            if assigned:
                nbAssignedPairs += 1
                if self.method == "3":
                    nbAssignedPairsTwoTags += t2
                    nbAssignedPairsOneTag += t1
                if ind not in dOutFqHandles:
                    dOutFqHandles[ind] = [
                        gzip.open("%s_%s_R1.fastq.gz" %
                                  (self.outFqPrefix, ind), "w"),
                        gzip.open("%s_%s_R2.fastq.gz" %
                                  (self.outFqPrefix, ind), "w")]
                dOutFqHandles[ind][0].write("@%s\n%s\n+\n%s\n" %
                                            (read1_id,
                                             read1_seq[idx1:],
                                             read1_q[idx1:]))
                dOutFqHandles[ind][1].write("@%s\n%s\n+\n%s\n" %
                                            (read2_id,
                                             read2_seq[idx2:],
                                             read2_q[idx2:]))
            else: # not assigned
                if "unassigned" not in dOutFqHandles:
                    dOutFqHandles["unassigned"] = [
                        gzip.open("%s_unassigned_R1.fastq.gz" %
                                  self.outFqPrefix, "w"),
                        gzip.open("%s_unassigned_R2.fastq.gz" %
                                  self.outFqPrefix, "w")]
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
        
        if self.verbose > 0:
            msg = "total nb of read pairs: %i" % nbPairs
            msg += "\nnb of assigned read pairs: %i" % nbAssignedPairs
            if "3" in self.method:
                msg += "; 2t=%i 1t=%i" % (nbAssignedPairsTwoTags, \
                                          nbAssignedPairsOneTag)
            msg += "\nnb of unassigned read pairs: %i (%.2f%%" % (
                (nbPairs - nbAssignedPairs),
                100 * (nbPairs - nbAssignedPairs) / float(nbPairs))
            msg += ")"
            nbInds = len(dOutFqHandles) - 1
            msg += "\nnb of individuals with assigned reads: %i" % nbInds
            print(msg); sys.stdout.flush()
            
            
    def run(self):
        self.loadTags()
        if self.method == "4":
            self.prepareRemainingMotifs()
        self.demultiplexPairedReads()
        
        
if __name__ == "__main__":
    i = Demultiplex()
    
    i.setAttributesFromCmdLine()
    
    i.checkAttributes()
    
    if i.verbose > 0:
        startTime = time.time()
        msg = "START %s %s" % (os.path.basename(sys.argv[0]),
                               time.strftime("%Y-%m-%d %H:%M:%S"))
        msg += "\ncmd-line: %s" % ' '.join(sys.argv)
        msg += "\ncwd: %s" % os.getcwd()
        print(msg); sys.stdout.flush()
        
    i.run()
    
    if i.verbose > 0:
        msg = "END %s %s" % (os.path.basename(sys.argv[0]),
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
