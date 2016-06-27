#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: demultiplex samples in fastq files
# Copyright (C) 2014-2016 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut], Laurène Gay [ctb], Nicolas Rode [ctb]
# Versioning: https://github.com/timflutre/quantgen

# Inspired from:
# https://gist.github.com/seandavi/3015625 by Sean Devi
# http://bcb.io/2009/08/09/trimming-adaptors-from-short-read-sequences/ by Brad Chapman
# http://news.open-bio.org/news/2009/09/biopython-fast-fastq/ by Peter Cock

# Tests:
# $ python -m doctest demultiplex.py
# $ ./test_demultiplex.py -p ~/src/quantgen/demultiplex.py

# TODO:
# when --subst > 0, allow to decrease nb subst depending on the pattern (i.e. not for all)
# when --subst > 0, if several possible assignments, choose best (lower subst, higher length)
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
try:
    import regex
    regexIsImported = True
except ImportError:
    regexIsImported = False
    
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Restriction
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "1.13.0" # http://semver.org/


class Demultiplex(object):
    """
    Class performing demultiplexing.
    """
    
    def __init__(self):
        self.verbose = 1
        self.inDir = "."
        self.inFqFile1 = ""
        self.inFqFile2 = ""
        self.tagFile = ""
        self.outPrefix = ""
        self.method = ""
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
        self.nbSubstitutionsAllowed = 0 # requires module "regex"
        self.enforceSubst = "lenient"
        self.onlyComparePatterns = False
        
        
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
        msg += "      --idir\tpath to the input directory with the fastq files\n"
        msg += "\t\tdefault=. (current directory)\n"
        msg += "\t\toptional with --compp\n"
        msg += "      --ifq1\tpath to the first input fastq file\n"
        msg += "\t\tcan be compressed with gzip\n"
        msg += "\t\toptional with --compp\n"
        msg += "      --ifq2\tpath to the second input fastq file\n"
        msg += "\t\tcan be compressed with gzip\n"
        msg += "\t\terror raised if reads not in same order as --ifq1\n"
        msg += "\t\toptional with --compp\n"
        msg += "      --it\tpath to the tag file\n"
        msg += "\t\tonly A, T, G and C (no ambiguous nucleotide allowed)\n"
        msg += "\t\tcan be in 2 formats (automatically detected)\n"
        msg += "\t\t fasta: put sample names in the fasta headers\n"
        msg += "\t\t table: 2 columns, header line should be 'id\\ttag'\n"
        msg += "\t\talways compulsory\n"
        msg += "      --ofqp\tprefix for the output fastq files\n"
        msg += "\t\t2 per ind, 1 for unassigned, 1 for chimeras\n"
        msg += "\t\tsuffix will be \".fastq.gz\" (not \".fq\" because of FastQC)\n"
        msg += "\t\twill be compressed with gzip\n"
        msg += "\t\tthis prefix will also be used for a gzipped text file\n"
        msg += "\t\t counting assigned read pairs per individual\n"
        msg += "\t\toptional with --compp\n"
        msg += "      --met\tmethod to assign pairs of reads to individuals via tags\n"
        msg += "\t\tonly forward strand is considered\n"
        msg += "\t\talways compulsory\n"
        msg += "\t\t1: assign pair if both reads start with the tag\n"
        msg += "\t\t2: assign pair if at least one read starts with the tag\n"
        msg += "\t\t3: same as 2 but count if one or both reads start with the tag\n"
        msg += "\t\t4a: assign pair if first read starts with tag (ignore second read)\n"
        msg += "\t\t4b: assign pair if first read has tag in its first N bases\n"
        msg += "\t\t    (ignore second read, see --dist)\n"
        msg += "\t\t4c: assign pair if first read has tag and remaining cut site\n"
        msg += "\t\t    in its first N bases (ignore second read, see --dist and --re)\n"
        msg += "\t\t4d: assign pair if first and/or second read has tag and\n"
        msg += "\t\t    remaining cut site in its first N bases (see --dist and --re)\n"
        msg += "\t\t    PCR chimeras (R1 tag is different than R2 tag) are detected\n"
        msg += "\t\t    and saved in distinct files than the unassigned\n"
        msg += "      --subst\tnumber of substitutions allowed\n"
        msg += "\t\tdefault=1\n"
        msg += "\t\tif > 0, the 'regex' module is required\n"
        msg += "      --ensubst\tenforce the nb of substitutions allowed\n"
        msg += "\t\tdefault=lenient/strict\n"
        msg += "\t\t'lenient' starts from the value given via '--subst'\n"
        msg += "\t\tand decreases it until all tags are distinguishable\n"
        msg += "      --dist\tdistance from the read start to search for the tag (in bp)\n"
        msg += "\t\tdefault=20\n"
        msg += "\t\tset to '0' to disable for --met 4b/4c/4d\n"
        msg += "      --re\tname of the restriction enzyme (e.g. 'ApeKI')\n"
        msg += "      --chim\tsearch if full restriction site found in R1 and/or R2\n"
        msg += "\t\tdefault=1, see --re\n"
        msg += "\t\t0: don't search (some chimeras may still be detected if --met 4d)\n"
        msg += "\t\t1: if chimera, count as such, try to assign, and save in same files as others\n"
        msg += "\t\t2: if chimera, don't even try to assign and save in distinct files\n"
        msg += "      --nci\tdo not clip the tag when saving the assigned reads\n"
        msg += "      --compp\tonly compare patterns to be searched\n"
        msg += "\t\tuseful to choose how to set '--subst'\n"
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
        msg += "Copyright (C) 2014-2016 Institut National de la Recherche Agronomique (INRA).\n"
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
                                        "ofqp=", "met=", "subst=", "ensubst=",
                                        "re=", "chim=", "dist=", "nci",
                                        "compp"])
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
            elif o == "--subst":
                self.nbSubstitutionsAllowed = int(a)
            elif o == "--ensubst":
                self.enforceSubst = a
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
            elif o == "--compp":
                self.onlyComparePatterns = True
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if not self.onlyComparePatterns and self.inDir == "":
            msg = "ERROR: missing compulsory option --idir"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not self.onlyComparePatterns and not os.path.exists(self.inDir):
            msg = "ERROR: can't find dir %s" % self.inDir
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not self.onlyComparePatterns and self.inFqFile1 == "":
            msg = "ERROR: missing compulsory option --ifq1"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        self.inFqFile1 = "%s/%s" % (self.inDir, self.inFqFile1)
        if not self.onlyComparePatterns and not os.path.exists(self.inFqFile1):
            msg = "ERROR: can't find file %s" % self.inFqFile1
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not self.onlyComparePatterns and self.inFqFile2 == "":
            msg = "ERROR: missing compulsory option --ifq2"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        self.inFqFile2 = "%s/%s" % (self.inDir, self.inFqFile2)
        if not self.onlyComparePatterns and not os.path.exists(self.inFqFile2):
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
        if not self.onlyComparePatterns and self.outPrefix == "":
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
        if self.nbSubstitutionsAllowed > 0 and not regexIsImported:
            msg = "ERROR: --subst > 0 but module 'regex' can't be imported"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.enforceSubst not in ["lenient", "strict"]:
            msg = "ERROR: --ensubst %s is unknown" % self.enforceSubst
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
            
    def findTagFileFormatFromLine(self, line):
        """
        >>> i = Demultiplex()
        >>> line = "id\ttag"
        >>> i.findTagFileFormatFromLine(line)
        u'table'
        >>> line = ">ind001"
        >>> i.findTagFileFormatFromLine(line)
        u'fasta'
        """
        tagFileFormat = None
        
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
    
    
    def findTagFileFormat(self):
        """
        Look for a '>' in the tag file and decide to return a 'fasta' or 'table' file format.
        """
        tagFileFormat = "table"
        tagHandle = open(self.tagFile)
        line = tagHandle.readline()
        tagHandle.close()
        tagFileFormat = self.findTagFileFormatFromLine(line)
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
        """
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.restrictEnzyme = Restriction.__dict__["ApeKI"]
        >>> i.retrieveRestrictionEnzyme()
        >>> i.cutMotif
        u'G^CWG_C'
        >>> i.regexpCutMotif
        u'GC[AT]GC'
        """
        if self.verbose > 0:
            print("retrieve restriction enzyme from Biopython...")
            
        self.cutMotif = u"%s" % self.restrictEnzyme.elucidate()
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
                
            self.regexpCompilCutMotif = re.compile(self.regexpCutMotif,
                                                   flags=re.IGNORECASE)
            
            
    def prepareRemainingMotif(self):
        """
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.restrictEnzyme = Restriction.__dict__["ApeKI"]
        >>> i.retrieveRestrictionEnzyme()
        >>> i.prepareRemainingMotif()
        >>> i.regexpRemainMotif
        u'C[AT]GC'
        >>> i.lenRemainMotif
        4
        """
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
            
            
    def checkDist(self):
        """
        Check that length of tag (+ remaining motif) <= self.dist.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'AAA', u'ind2': u'TTT'}
        >>> i.dist = 4
        >>> i.checkDist()
        >>> i.method = u'4c'
        >>> i.restrictEnzyme = Restriction.__dict__["ApeKI"]
        >>> i.retrieveRestrictionEnzyme()
        >>> i.prepareRemainingMotif()
        >>> i.checkDist()
        Traceback (most recent call last):
        ValueError: --dist 4 is too short for tag ind1
        because the whole pattern is AAAC[AT]GC
        """
        if self.verbose > 0:
            msg = "check that '--dist' is compatible with pattern lengths..."
            print(msg); sys.stdout.flush()
            
        for tagId in self.tags:
            pattern = self.tags[tagId]
            if self.method in ["4c", "4d"]:
                pattern += self.regexpRemainMotif
            if self.dist > 0 and len(pattern) > self.dist:
                msg = "--dist %i is too short for tag %s" % (self.dist,
                                                                    tagId)
                msg += "\nbecause the whole pattern is %s" % pattern
                raise ValueError(msg)
            
            
    def compilePatterns(self):
        """
        Compile patterns, that is each tag (+ remaining right of cut site, depending on self.method) as regular expression for quicker searches.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'AAA'}
        >>> i.restrictEnzyme = Restriction.__dict__["ApeKI"]
        >>> i.retrieveRestrictionEnzyme()
        >>> i.prepareRemainingMotif()
        >>> i.method = u'1'
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'^AAA'
        >>> i.method = u'2'
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'^AAA'
        >>> i.method = u'3'
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'^AAA'
        >>> i.method = u'4a'
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'^AAA'
        >>> i.method = u'4b'
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'AAA'
        >>> i.method = u'4c'
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'AAAC[AT]GC'
        >>> i.method = u'4d'
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'AAAC[AT]GC'
        >>> i.nbSubstitutionsAllowed = 1
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'(AAAC[AT]GC){s<=1}'
        >>> i.dist = 0
        >>> i.compilePatterns()
        >>> i.patterns["ind1"].pattern
        u'^(AAAC[AT]GC){s<=1}'
        """
        if self.verbose > 0:
            msg = "compile patterns"
            msg += " (method %s" % self.method
            msg += ", distance %i bp" % self.dist
            msg += ", %i substitution" % self.nbSubstitutionsAllowed
            if self.nbSubstitutionsAllowed > 1:
                msg += "s"
            msg += " allowed)..."
            print(msg); sys.stdout.flush()
            
        for tagId in self.tags:
            
            ## build the pattern
            pattern = self.tags[tagId]
            if self.method in ["4c", "4d"]:
                pattern += self.regexpRemainMotif
            if self.nbSubstitutionsAllowed > 0:
                pattern = "(%s){s<=%i}" % (pattern,
                                           self.nbSubstitutionsAllowed)
            if self.method in ["1", "2", "3", "4a"] or self.dist == 0:
                pattern = "^" + pattern
                
            ## compile the pattern
            if self.nbSubstitutionsAllowed > 0:
                self.patterns[tagId] = regex.compile(pattern,
                                                     flags=regex.IGNORECASE)
            else:
                self.patterns[tagId] = re.compile(pattern,
                                                  flags=re.IGNORECASE)
                
                
    def makeSeqsToCompareTwoTags(self, tag2):
        """
        Return a list of DNA sequence(s) as string(s) to which "patseq" should be matched against.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.method = "1"
        >>> i.makeSeqsToCompareTwoTags("AAA")
        [u'AAA']
        >>> i.method = "4c"
        >>> i.regexpRemainMotif = "CAGC"
        >>> i.makeSeqsToCompareTwoTags("AAA")
        [u'AAACAGC']
        >>> i.regexpRemainMotif = "C[AT]GC"
        >>> i.lenRemainMotif = 4
        >>> i.makeSeqsToCompareTwoTags("AAA")
        [u'AAACAGC', u'AAACTGC']
        """
        seqs = []
        if self.method not in ["4c", "4d"]:
            seqs.append(tag2)
        else:
            if re.search("^[ATGCN]+$", self.regexpRemainMotif) is not None:
                seqs.append(tag2 + self.regexpRemainMotif)
            elif re.search("^[ATGCN\[\]]+$", self.regexpRemainMotif) is not None:
                idx1 = self.regexpRemainMotif.find("[")
                idx2 = self.regexpRemainMotif.find("]")
                nbSeqs = idx2 - (idx1 + 1)
                for idxSeq in range(nbSeqs):
                    seqs.append(tag2)
                    idxNt = 0
                    while True:
                        if idxNt >= len(self.regexpRemainMotif):
                            break
                        if idxNt < idx1 or idxNt > idx2:
                            seqs[idxSeq] += self.regexpRemainMotif[idxNt]
                        elif idxNt > idx1 and idxNt < idx2:
                            seqs[idxSeq] += self.regexpRemainMotif[idxNt+idxSeq]
                            idxNt = idx2
                        idxNt += 1
            else:
                msg = "self.regexpRemainMotif %s contains other symbols" % \
                      self.regexpRemainMotif
                msg += "\nthan only A, T, G, C, [, and ]"
                raise ValueError(msg)
        return seqs
    
    
    def comparePatternsTwoTags(self, tagId1, tagId2, seqs):
        """
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'TTAC', u'ind2': u'TTAGCTT'}
        >>> i.method = u'4c'
        >>> i.nbSubstitutionsAllowed = 2
        >>> i.restrictEnzyme = Restriction.__dict__["ApeKI"]
        >>> i.retrieveRestrictionEnzyme()
        >>> i.prepareRemainingMotif()
        >>> i.compilePatterns()
        >>> seqs = [u'TTAGCTTCAGC', u'TTAGCTTCTGC']
        >>> i.comparePatternsTwoTags(u'ind1', u'ind2', seqs)
        Traceback (most recent call last):
        ValueError: with 2 allowed substitutions, tag TTAC corresponding to ind1
        is indistinguishable from tag TTAGCTT corresponding to ind2
        pattern: (TTACC[AT]GC){s<=2}
        string: TTAGCTTCAGC
        start-end: 0-8
        >>> i.tags = {u'ind1': u'GTGGA', u'ind2': u'TTGCAGGA'}
        >>> i.compilePatterns()
        >>> seqs = [u'TTGCAGGACAGC', u'TTGCAGGACTGC']
        >>> i.comparePatternsTwoTags(u'ind1', u'ind2', seqs)
        Traceback (most recent call last):
        ValueError: with 2 allowed substitutions, tag GTGGA corresponding to ind1
        is indistinguishable from tag TTGCAGGA corresponding to ind2
        pattern: (GTGGAC[AT]GC){s<=2}
        string: TTGCAGGACAGC
        start-end: 3-12
        >>> i.dist = 0
        >>> i.compilePatterns()
        >>> i.comparePatternsTwoTags(u'ind1', u'ind2', seqs)
        """
        for seq in seqs:
            tmpRe = self.patterns[tagId1].search(seq)
            if tmpRe:
                msg = "with %i allowed substitution" % \
                      self.nbSubstitutionsAllowed
                if self.nbSubstitutionsAllowed > 1:
                    msg += "s"
                    msg += ", tag %s" % self.tags[tagId1]
                    msg += " corresponding to %s" % tagId1
                    msg += "\nis indistinguishable from tag %s" % self.tags[tagId2]
                    msg += " corresponding to %s" % tagId2
                    msg += "\npattern: %s" % self.patterns[tagId1].pattern
                    msg += "\nstring: %s" % seq
                    msg += "\nstart-end: %i-%i" % (tmpRe.start(), tmpRe.end())
                    raise ValueError(msg)
                
                
    def comparePatterns(self):
        """
        Raise an exception if two patterns of the same length are indistinguishable when allowing substitutions.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'TTAC', u'ind2': u'TTAGCTT'}
        >>> i.method = u'4c'
        >>> i.nbSubstitutionsAllowed = 2
        >>> i.restrictEnzyme = Restriction.__dict__["ApeKI"]
        >>> i.retrieveRestrictionEnzyme()
        >>> i.prepareRemainingMotif()
        >>> i.compilePatterns()
        >>> i.comparePatterns()
        Traceback (most recent call last):
        ValueError: with 2 allowed substitutions, tag TTAC corresponding to ind1
        is indistinguishable from tag TTAGCTT corresponding to ind2
        pattern: (TTACC[AT]GC){s<=2}
        string: TTAGCTTCAGC
        start-end: 0-8
        """
        if self.nbSubstitutionsAllowed > 0:
            if self.verbose > 0:
                print("check that searched patterns are distinguishable...")
                sys.stdout.flush()
                
            lTagIds = self.tags.keys()
            lTagIds.sort()
            
            for i in range(0, len(lTagIds)):
                patseq = self.patterns[lTagIds[i]].pattern.split(")")[0].replace("(","")
                
                for j in range(0, len(lTagIds)):
                    if j == i:
                        continue
                    
                    lenSeq = len(self.tags[lTagIds[j]])
                    if self.method in ["4c", "4d"]:
                        lenSeq += self.lenRemainMotif
                    if len(patseq) > lenSeq:
                        continue
                    
                    seqs = self.makeSeqsToCompareTwoTags(self.tags[lTagIds[j]])
                    
                    self.comparePatternsTwoTags(lTagIds[i], lTagIds[j], seqs)
                    
                    
    def flexCompileComparePatterns(self):
        """
        Run compilePatterns() and comparePatterns() by decreasing the number of substitutions allowed, until it doesn't raise any exception, if possible.
        """
        for nbSubsts in range(self.nbSubstitutionsAllowed, -1, -1):
            self.nbSubstitutionsAllowed = nbSubsts
            
            try:
                self.compilePatterns()
                self.comparePatterns()
            except Exception as e:
                print(e.args[0])
                if nbSubsts > 0:
                    continue
                else:
                    msg = "some tags are indistinguishable even when no substitution is allowed"
                    e.args += (msg,)
                    raise e
                
            if self.verbose > 0:
                msg = "nb of substitutions allowed: %i" % self.nbSubstitutionsAllowed
                print(msg); sys.stdout.flush()
                
            break
        
        
    def identifyChimeras(self, read1_seq, read2_seq):
        chimera = False
        if self.regexpCompilCutMotif.search(read1_seq) != None \
           or self.regexpCompilCutMotif.search(read2_seq) != None:
            chimera = True
        return chimera
        
        
    def identifyIndividual_1(self, read1_seq, read2_seq):
        """
        If both reads start with the tag, return "assign" as True, "ind" as the identifier of the individual to which reads are assigned, and "idx" as the position at which both reads should be saved.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'AAA', u'ind2': u'TTT'}
        >>> i.method = u'1'
        >>> i.compilePatterns()
        >>> i.clipIdx = False
        >>> i.identifyIndividual_1(u'AAANNNNN', u'AAANNNNN')
        (True, u'ind1', 0, 0)
        >>> i.identifyIndividual_1(u'AAANNNNN', u'NNNNNNNN')
        (False, None, 0, 0)
        >>> i.clipIdx = True
        >>> i.identifyIndividual_1(u'AAANNNNN', u'AAANNNNN')
        (True, u'ind1', 3, 3)
        """
        assigned = False
        ind = None
        idx = 0
        
        for tagId in self.tags:
            
            if self.patterns[tagId].search(read1_seq) \
               and self.patterns[tagId].search(read2_seq):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx = len(self.tags[tagId])
                break
            
        return assigned, ind, idx, idx
    
    
    def identifyIndividual_2(self, read1_seq, read2_seq):
        """
        If at least one read starts with the tag, return "assign" as True, "ind" as the identifier of the individual to which reads are assigned, and idx1 and/or idx2 as the position(s) at which the read(s) should be saved. Note that "idx1 != 0 and idx2 = 0", or "idx1 = 0 and idx2 != 0", if the tag is found in only one of the reads.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'AAA', u'ind2': u'TTT'}
        >>> i.method = u'2'
        >>> i.compilePatterns()
        >>> i.clipIdx = False
        >>> i.identifyIndividual_2(u'AAANNNNN', u'AAANNNNN')
        (True, u'ind1', 0, 0)
        >>> i.identifyIndividual_2(u'AAANNNNN', u'NNNNNNNN')
        (True, u'ind1', 0, 0)
        >>> i.clipIdx = True
        >>> i.identifyIndividual_2(u'AAANNNNN', u'AAANNNNN')
        (True, u'ind1', 3, 3)
        >>> i.identifyIndividual_2(u'AAANNNNN', u'NNNNNNNN')
        (True, u'ind1', 3, 0)
        """
        assigned = False
        ind = None
        idx1 = 0
        idx2 = 0
        
        for tagId in self.tags:
            
            if self.patterns[tagId].search(read1_seq):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = len(self.tags[tagId])
                    if self.patterns[tagId].search(read2_seq):
                        idx2 = len(self.tags[tagId])
                break
            
            if self.patterns[tagId].search(read2_seq):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx2 = len(self.tags[tagId])
                break
            
        return assigned, ind, idx1, idx2
    
    
    def identifyIndividual_3(self, read1_seq, read2_seq):
        """
        Same as identifyIndividual_2(), but count if one or both reads start with the tag.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'AAA', u'ind2': u'TTT'}
        >>> i.method = u'3'
        >>> i.compilePatterns()
        >>> i.clipIdx = False
        >>> i.identifyIndividual_3(u'AAANNNNN', u'AAANNNNN')
        (True, u'ind1', 0, 0, 1, 0)
        >>> i.identifyIndividual_3(u'AAANNNNN', u'NNNNNNNN')
        (True, u'ind1', 0, 0, 0, 1)
        >>> i.clipIdx = True
        >>> i.identifyIndividual_3(u'AAANNNNN', u'AAANNNNN')
        (True, u'ind1', 3, 3, 1, 0)
        >>> i.identifyIndividual_3(u'AAANNNNN', u'NNNNNNNN')
        (True, u'ind1', 3, 0, 0, 1)
        """
        assigned = False
        ind = None
        idx1 = 0
        idx2 = 0
        nbAssignedPairsTwoTags = 0
        nbAssignedPairsOneTag = 0
        
        for tagId in self.tags:
            
            if self.patterns[tagId].search(read1_seq):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = len(self.tags[tagId])
                if self.patterns[tagId].search(read2_seq):
                    if self.clipIdx:
                        idx2 = len(self.tags[tagId])
                    nbAssignedPairsTwoTags += 1
                else:
                    nbAssignedPairsOneTag += 1
                break
            
            if self.patterns[tagId].search(read2_seq):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx2 = len(self.tags[tagId])
                nbAssignedPairsOneTag += 1
                break
            
        return assigned, ind, idx1, idx2, nbAssignedPairsTwoTags, \
            nbAssignedPairsOneTag
    
    
    def identifyIndividual_4a(self, read1_seq):
        """
        If the first read starts with the tag, return "assign" as True, "ind" as the identifier of the individual to which reads are assigned, and idx1 as the position at which the first read should be saved. Note that the second read is ignored, thus idx2 = 0.
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'AAA', u'ind2': u'TTT'}
        >>> i.method = u'4a'
        >>> i.compilePatterns()
        >>> i.clipIdx = False
        >>> i.identifyIndividual_4a(u'AAANNNNN')
        (True, u'ind1', 0, 0)
        >>> i.identifyIndividual_4a(u'NNNNNNNN')
        (False, None, 0, 0)
        >>> i.clipIdx = True
        >>> i.identifyIndividual_4a(u'AAANNNNN')
        (True, u'ind1', 3, 0)
        """
        assigned = False
        ind = None
        idx1 = 0
        
        for tagId in self.tags:
            
            if self.patterns[tagId].search(read1_seq):
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = len(self.tags[tagId])
                break
            
        return assigned, ind, idx1, 0
    
    
    def identifyIndividual_4b(self, read1_seq):
        """
        If the first read has a tag in its first N bases (specified via self.dist), return "assign" as True, "ind" as the identifier of the individual to which reads are assigned, and idx1 as the position at which the first read should be saved. Note that the second read is ignored, thus idx2 = 0.
        """
        assigned = False
        ind = None
        idx1 = 0
        
        for tagId in self.tags:
            
            tmpRe = self.patterns[tagId].search(read1_seq[:self.dist])
            if tmpRe != None:
                assigned = True
                ind = tagId
                if self.clipIdx:
                    idx1 = tmpRe.end()
                break
            
        return assigned, ind, idx1, 0
        
        
    def identifyIndividual_4c(self, read1_seq):
        """
        Assign pair if first read has tag and remaining cut site in its first
        N bases (ignore second read, use self.dist and self.patterns).
        
        >>> i = Demultiplex()
        >>> i.verbose = 0
        >>> i.tags = {u'ind1': u'AAA', u'ind2': u'TTT'}
        >>> i.restrictEnzyme = Restriction.__dict__["ApeKI"]
        >>> i.retrieveRestrictionEnzyme()
        >>> i.prepareRemainingMotif()
        >>> i.method = u'4c'
        >>> i.compilePatterns()
        >>> i.clipIdx = False
        >>> i.identifyIndividual_4c(u'NNNNNNNNNN')
        (False, None, 0, 0)
        >>> i.identifyIndividual_4c(u'AAANNNNNNN')
        (False, None, 0, 0)
        >>> i.identifyIndividual_4c(u'AAACAGCNNN')
        (True, u'ind1', 0, 0)
        >>> i.identifyIndividual_4c(u'AAACTGCNNN')
        (True, u'ind1', 0, 0)
        >>> i.clipIdx = True
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
            
        return assigned, ind, idx1, 0
    
    
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
            
        return assigned, ind, idx1, idx2, AssignedPairsTwoTags, \
            AssignedPairsOneTag, chimera
    
    
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
            msg = "demultiplex paired-end reads (method=%s, subst=%i)..." % \
                  (self.method, self.nbSubstitutionsAllowed)
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
            
        self.checkDist()
        
        if self.enforceSubst == "lenient":
            self.flexCompileComparePatterns()
        elif self.enforceSubst == "strict":
            self.compilePatterns()
            self.comparePatterns()
            
        if not self.onlyComparePatterns:
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
