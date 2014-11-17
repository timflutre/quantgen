#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: test demultiplex.py
# Author: Timothée Flutre
# Copyright (C) 2014 Institut National de la Recherche Agronomique
# License: GPL-3+

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import sys
import os
import getopt
import time
import datetime
from subprocess import Popen, PIPE, check_output
import math
import gzip
import shutil
import itertools
# import numpy as np
# import scipy as sp

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data.IUPACData import ambiguous_dna_values

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
        
class TestDemultiplex(object):
    
    def __init__(self):
        self.verbose = 1
        self.pathToProg = ""
        self.clean = True
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' tests demultiplex.py.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "  -p, --p2p\tfull path to the program to be tested\n"
        msg += "  -n, --noclean\tkeep temporary directory with all files\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -p ~/src/demultiplex.py\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        """
        msg = "%s 1.0\n" % os.path.basename(sys.argv[0])
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
            opts, args = getopt.getopt( sys.argv[1:], "hVv:p:n",
                                        ["help", "version", "verbose=",
                                         "p2p=", "noclean"])
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
            
            
    def launchProg(self, ifq1, ifq2, it, met):
        args = [self.pathToProg,
                "--idir", "./",
                "--ifq1", ifq1,
                "--ifq2", ifq2,
                "--it", it,
                "--ofqp", "test",
                "--met", met,
                "-v", str(self.verbose - 1)]
        if self.verbose > 0:
            print(" ".join(args))
        msgs = check_output(args)
        return msgs
        
        
    def test_met1_prepare(self):
        ifq1 = "reads_R1.fastq.gz"
        ifq2 = "reads_R2.fastq.gz"
        ifq1Handle = gzip.open(ifq1, "w")
        ifq2Handle = gzip.open(ifq2, "w")
        
        # both reads have perfect tag of ind 2
        txt = "@INST1:1:FLOW1:2:2104:15343:197393 1:N:0\n"
        txt += "TTT"
        txt += "TCAACCTGGAGTTCCAC\n"
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197393 2:N:0\n"
        txt += "TTT"
        txt += "GTAGCTGAGATCGGAAG\n"
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        # only read 1 has perfect tag of ind 1
        txt = "@INST1:1:FLOW1:2:2104:15343:197393 1:N:0\n"
        txt += "AAA"
        txt += "TCAACCTGGAGTTCCAC\n"
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197393 2:N:0\n"
        txt += "ATA"
        txt += "GTAGCTGAGATCGGAAG\n"
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        ifq1Handle.close()
        ifq2Handle.close()
        
        for f in ["test_ind2_R1.fastq.gz", "test_ind2_R2.fastq.gz",
                  "test_unassigned_R1.fastq.gz", "test_unassigned_R2.fastq.gz"]:
            if os.path.isfile(f):
                os.remove(f)
                
        it = "tags.fa"
        with open(it, "w") as itHandle:
            txt = ">ind1\nAAA\n"
            txt += ">ind2\nTTT\n"
            txt += ">ind3\nGGG\n"
            itHandle.write(txt)
        
        return ifq1, ifq2, it
        
        
    def test_met1_comp(self, msgs):
        if not os.path.exists("test_ind2_R1.fastq.gz") or \
           not os.path.exists("test_ind2_R2.fastq.gz"):
            print("test_met1: fail (1)")
        else:
            with gzip.open("test_ind2_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind2_R2.fastq.gz") as inFqHandle2:
                reads1 = SeqIO.parse(inFqHandle1, "fastq",
                                     alphabet=IUPAC.ambiguous_dna)
                reads2 = SeqIO.parse(inFqHandle2, "fastq",
                                     alphabet=IUPAC.ambiguous_dna)
                for (read1, read2) in itertools.izip(reads1, reads2):
                    if read1.id != "INST1:1:FLOW1:2:2104:15343:197393" and \
                       read2.id != "INST1:1:FLOW1:2:2104:15343:197393":
                        print("test_met1: fail (2)")
            print("test_met1: pass")
            
            
    def test_met1(self):
        if self.verbose > 0:
            print("launch test 1 ...")
            sys.stdout.flush()
        ifq1, ifq2, it = self.test_met1_prepare()
        msgs = self.launchProg(ifq1, ifq2, it, "1")
        self.test_met1_comp(msgs)
        
        
    def run(self):
        cwd = os.getcwd()
        uniqId = os.getpid()
        testDir = "tmp_test_%s" % uniqId
        if os.path.exists(testDir):
            shutil.rmtree(testDir)
        os.mkdir(testDir)
        os.chdir(testDir)
        if self.verbose > 0:
            print("temp dir: %s" % os.getcwd()); sys.stdout.flush()
            
        self.test_met1()
        
        os.chdir(cwd)
        if self.clean:
            shutil.rmtree(testDir)
            
            
if __name__ == "__main__":
    i = TestDemultiplex()
    
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
