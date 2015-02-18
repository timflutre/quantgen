#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: test trimfilter.py
# Copyright (C) 2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Author: Timothée Flutre
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
from subprocess import Popen, PIPE, check_output
import math
import gzip
import tempfile
import shutil
import itertools

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data.IUPACData import ambiguous_dna_values

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "0.1.0" # http://semver.org/


class TestTrimfilter(object):
    
    def __init__(self):
        self.verbose = 0
        self.pathToProg = ""
        self.testsToRun = ["exact", "fuzzy"]
        self.clean = True
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' tests trimfilter.py.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (default=0/1/2/3)\n"
        msg += "  -p, --p2p\tfull path to the program to be tested\n"
        msg += "  -t, --test\tidentifiers of test(s) to run (default=exact-fuzzy)\n"
        msg += "  -n, --noclean\tkeep temporary directory with all files\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -p ~/src/trimfilter.py -n\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA).\n"
        msg += "License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:p:t:n",
                                        ["help", "version", "verbose=",
                                         "p2p=", "test=", "noclean"])
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
            if t not in ["exact", "fuzzy"]:
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
        
        
    def writeAdpFile(self, adpFile):
        """
        The adapter is 5'-AAA-3', for both R1 and R2 reads. Thus, in 
        case of read-through, we should look for TTT in the 3' end of both R1
        and R2. Moreover, for R2, we should also look for the barcode, which is
        5'-TAG-3' in R1, thus 5'-CTA-3' for R2.
        """
        with open(adpFile, "w") as handle:
            txt = ">adp2_rc/1\nTTT\n" # reverse-complement of the R2 adapter
            txt += ">adp1_rc/2\nCTATTT\n" # revcomp of R1 adp preceded by the tag
            handle.write(txt)
            
            
    def launchProg(self, ifq1, ifq2, adpFile, maxNbErrors):
        """
        Launch trimfilter.py with the specified arguments
        """
        args = [self.pathToProg,
                "--idir", "./",
                "--ifq1", ifq1,
                "--ifq2", ifq2,
                "--adp", adpFile,
                "--err", str(maxNbErrors),
                "--op", "test",
                "-v", str(self.verbose - 1)]
        if self.verbose > 0:
            print(" ".join(args))
        # msgs = check_output(args)
        msgs = Popen(args, stdout=PIPE, stderr=PIPE).communicate()
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
    
    
    def test_exact_prepare(self):
        ifq1 = "reads_R1.fq.gz"
        ifq2 = "reads_R2.fq.gz"
        ifq1Handle = gzip.open(ifq1, "w")
        ifq2Handle = gzip.open(ifq2, "w")
        
        # pair 1: both reads have exact adapter read-through
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 1:N:0:\n"
        txt += "GGGGGGGGGG" # insert
        txt += "TTT\n" # full adp read-through
        txt += "+\n"
        txt += "~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 2:N:0:\n"
        txt += "CCCCCCCCCC" # insert
        txt += "CTATTT\n" # full tag + adp read-through
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        # pair 2: both reads have a partial adapter read-through
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 1:N:0:\n"
        txt += "GGGGGGGGGG" # insert
        txt += "TT\n" # partial adp read-through
        txt += "+\n"
        txt += "~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 2:N:0:\n"
        txt += "CCCCCCCCCC" # insert
        txt += "CTATT\n" # partial tag + adp read-through
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        ifq1Handle.close()
        ifq2Handle.close()
        
        for f in ["test_paired_R1.fq.gz", "test_unpaired_R1.fq.gz",
                  "test_paired_R2.fq.gz", "test_unpaired_R2.fq.gz"]:
            if os.path.isfile(f):
                os.remove(f)
                
        adpFile = "adapters.fa"
        self.writeAdpFile(adpFile)
        
        return ifq1, ifq2, adpFile
        
        
    def test_exact_comp(self, msgs):
        testId = 1
        if not os.path.exists("test_paired_R1.fq.gz") \
           or not os.path.exists("test_unpaired_R1.fq.gz") \
           or not os.path.exists("test_paired_R2.fq.gz") \
           or not os.path.exists("test_unpaired_R2.fq.gz") \
           or not os.path.exists("test.log.gz"):
            print("test_exact: fail (%i)" % testId)
            return
        testId += 1
        with gzip.open("test.log.gz") as logHandle:
            lines = logHandle.readlines()
            if lines[1] != "INST1:1:FLOW1:2:2104:15343:197391 1:N:0:\tadp2_rc/1\t10\n":
                print("test_exact: fail (%i)" % testId)
                print(line[1])
                return
            testId += 1
            if lines[2] != "INST1:1:FLOW1:2:2104:15343:197391 2:N:0:\tadp1_rc/2\t10\n":
                print("test_exact: fail (%i)" % testId)
                print(line[2])
                return
            testId += 1
            if lines[3] != "INST1:1:FLOW1:2:2104:15343:197392 1:N:0:\tadp2_rc/1\t-1\n":
                print("test_exact: fail (%i)" % testId)
                print(line[3])
                return
            testId += 1
            if lines[4] != "INST1:1:FLOW1:2:2104:15343:197392 2:N:0:\tadp1_rc/2\t-1\n":
                print("test_exact: fail (%i)" % testId)
                print(line[4])
                return
            testId += 1
        with gzip.open("test_paired_R1.fq.gz") as inFqHandle1, \
             gzip.open("test_paired_R2.fq.gz") as inFqHandle2:
            reads1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
            reads2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
            if reads1[0].id != "INST1:1:FLOW1:2:2104:15343:197391" and \
               reads2[0].id != "INST1:1:FLOW1:2:2104:15343:197391":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[0].id))
                print(str(reads2[0].id))
                return
            testId += 1
            if str(reads1[0].seq) != "GGGGGGGGGG" or \
               str(reads2[0].seq) != "CCCCCCCCCC":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[0].seq))
                print(str(reads2[0].seq))
                return
            testId += 1
            if reads1[1].id != "INST1:1:FLOW1:2:2104:15343:197392" and \
               reads2[1].id != "INST1:1:FLOW1:2:2104:15343:197392":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[1].id))
                print(str(reads2[1].id))
                return
            testId += 1
            if str(reads1[1].seq) != "GGGGGGGGGGTT" or \
               str(reads2[1].seq) != "CCCCCCCCCCCTATT":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[1].seq))
                print(str(reads2[1].seq))
                return
            testId += 1
        print("test_exact: pass")
        
        
    def test_exact(self):
        if self.verbose > 0:
            print("launch test_exact ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, adpFile = self.test_exact_prepare()
        msgs = self.launchProg(ifq1, ifq2, adpFile, 0)
        self.test_exact_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================
    
    
    def test_fuzzy_comp(self, msgs):
        testId = 1
        if not os.path.exists("test_paired_R1.fq.gz") \
           or not os.path.exists("test_unpaired_R1.fq.gz") \
           or not os.path.exists("test_paired_R2.fq.gz") \
           or not os.path.exists("test_unpaired_R2.fq.gz") \
           or not os.path.exists("test.log.gz"):
            print("test_exact: fail (%i)" % testId)
            return
        testId += 1
        with gzip.open("test.log.gz") as logHandle:
            lines = logHandle.readlines()
            if lines[1] != "INST1:1:FLOW1:2:2104:15343:197391 1:N:0:\tadp2_rc/1\t10\n":
                print("test_exact: fail (%i)" % testId)
                print(line[1])
                return
            testId += 1
            if lines[2] != "INST1:1:FLOW1:2:2104:15343:197391 2:N:0:\tadp1_rc/2\t10\n":
                print("test_exact: fail (%i)" % testId)
                print(line[2])
                return
            testId += 1
            if lines[3] != "INST1:1:FLOW1:2:2104:15343:197392 1:N:0:\tadp2_rc/1\t10\n":
                print("test_exact: fail (%i)" % testId)
                print(line[3])
                return
            testId += 1
            if lines[4] != "INST1:1:FLOW1:2:2104:15343:197392 2:N:0:\tadp1_rc/2\t10\n":
                print("test_exact: fail (%i)" % testId)
                print(line[4])
                return
            testId += 1
        with gzip.open("test_paired_R1.fq.gz") as inFqHandle1, \
             gzip.open("test_paired_R2.fq.gz") as inFqHandle2:
            reads1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
            reads2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
            if reads1[0].id != "INST1:1:FLOW1:2:2104:15343:197391" and \
               reads2[0].id != "INST1:1:FLOW1:2:2104:15343:197391":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[0].id))
                print(str(reads2[0].id))
                return
            testId += 1
            if str(reads1[0].seq) != "GGGGGGGGGG" or \
               str(reads2[0].seq) != "CCCCCCCCCC":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[0].seq))
                print(str(reads2[0].seq))
                return
            testId += 1
            if reads1[1].id != "INST1:1:FLOW1:2:2104:15343:197392" and \
               reads2[1].id != "INST1:1:FLOW1:2:2104:15343:197392":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[1].id))
                print(str(reads2[1].id))
                return
            testId += 1
            if str(reads1[1].seq) != "GGGGGGGGGG" or \
               str(reads2[1].seq) != "CCCCCCCCCC":
                print("test_exact: fail (%i)" % testId)
                print(str(reads1[1].seq))
                print(str(reads2[1].seq))
                return
            testId += 1
        print("test_fuzzy: pass")
        
        
    def test_fuzzy(self):
        if self.verbose > 0:
            print("launch test_fuzzy ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, adpFile = self.test_exact_prepare()
        msgs = self.launchProg(ifq1, ifq2, adpFile, 1)
        self.test_fuzzy_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================
    
    
    def run(self):
        if "exact" in self.testsToRun:
            self.test_exact()
        if "fuzzy" in self.testsToRun:
            self.test_fuzzy()
            
            
if __name__ == "__main__":
    i = TestTrimfilter()
    
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
