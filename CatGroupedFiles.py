#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: concatenates files per group
# Copyright (C) 2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/quantgen

# TODO:
# read --input from stdin to allow usage with pipes

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

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "1.1.0" # http://semver.org/


def progressBar(progress):
    """
    Displays or updates a console progress bar.
    http://stackoverflow.com/a/15860757/597069
    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'. A value at 1 or bigger represents 100%.
    """
    barLength = 60
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "ERROR: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        # status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\r[{0}] {1:.2f}% {2}".format("="*block + "-"*(barLength-block),
                                     progress*100, status)
    sys.stdout.write(text)
    if progress == 1:
        sys.stdout.write("\n")
    sys.stdout.flush()
    
    
class CatGroupedFiles(object):
    
    def __init__(self):
        self.verbose = 1
        self.inFile = ""
        self.suffix = ""
        self.outDir = ""
        self.group2files = {}
        self.useSymLink = True
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' concatenates files per group.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "  -i, --input\tpath to the input file (2 columns sep. by a tab)\n"
        msg += "\t\tall files (col 2) with same group (col 1) will be\n"
        msg += "\t\t concatenated to a file named <col1>.<suffix>\n"
        msg += "\t\texample from sequencing applications:\n"
        msg += "\t\t indA<tab>run1/A.fastq.gz\n"
        msg += "\t\t indB<tab>run1/B.fastq.gz\n"
        msg += "\t\t indA<tab>run2/A.fastq.gz\n"
        msg += "  -s, --suffix\tsuffix for the output files (e.g. txt, fastq.gz, etc)\n"
        msg += "  -o, --outdir\toutput directory (default=\"\")\n"
        msg += "  -c, --copy\tcopy if group with single file (symlink otherwise)\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -i files_per_sample.txt -s txt\n" % os.path.basename(sys.argv[0])
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
        msg += "Copyright (C) 2015 Institut National de la Recherche Agronomique.\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "This is free software; see the source for copying conditions.  There is NO\n"
        msg += "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:i:s:o:c",
                                        ["help", "version", "verbose=",
                                         "input=", "suffix=", "output=",
                                         "copy="])
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
            elif o == "-i" or o == "--input":
                 self.inFile = a
            elif o == "-s" or o == "--suffix":
                self.suffix = a
            elif o == "-o" or o == "--output":
                 self.outDir = a
            elif o == "-c" or o == "--copy":
                self.useSymLink = False
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if self.inFile == "":
            msg = "ERROR: missing compulsory option --input"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.inFile):
            msg = "ERROR: can't find file %s" % self.inFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if len(self.suffix) == 0:
            msg = "ERROR: missing compulsory option --suffix"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.suffix[0] == ".":
            self.suffix = self.suffix[1:]
        if self.outDir != "" and not os.path.exists(self.outDir):
            msg = "ERROR: can't find directory %s" % self.outDir
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
            
    def loadInputFile(self):
        if self.verbose > 0:
            print("load input file ...")
            sys.stdout.flush()
        inH = open(self.inFile)
        lines = inH.readlines()
        nbFiles = 0
        for i in range(0,len(lines)):
            line = lines[i]
            tokens = line.rstrip().split("\t")
            if len(tokens) != 2:
                msg = "ERROR: line %i of %s should be as <col1><tab><col2>" \
                      % (i+1, self.inFile)
                sys.stderr.write("%s\n\n" % msg)
                sys.exit(1)
            if not os.path.exists(tokens[1]):
                msg = "ERROR: can't find file %s (line %i of %s)" \
                      % (tokens[1], i+1, self.inFile)
                sys.stderr.write("%s\n\n" % msg)
                sys.exit(1)
            # if not tokens[1].endswith(".gz"):
            #     msg = "ERROR: file %s should be gzipped (line %i of %s)" \
            #           % (tokens[1], i+1, self.inFile)
            #     sys.stderr.write("%s\n\n" % msg)
            #     sys.exit(1)
            if tokens[0] not in self.group2files:
                self.group2files[tokens[0]] = []
            self.group2files[tokens[0]].append(tokens[1])
            nbFiles += 1
        inH.close()
        if self.verbose > 0:
            print("%i groups ; %i files" % (len(self.group2files), nbFiles))
            
            
    def handleOneGroup(self, group, lFiles):
        """
        >>> i = CatGroupedFiles()
        >>> i.useSymLink = True; i.outDir = ""; i.suffix = "txt"
        >>> group = "indA"; lFiles=["run1/A.txt"]
        >>> i.handleOneGroup(group, lFiles)
        u'ln -s run1/A.txt indA.txt'
        >>> lFiles=["run1/A.txt", "run2/A.txt"]
        >>> i.handleOneGroup(group, lFiles)
        u'cat run1/A.txt run2/A.txt > indA.txt'
        """
        cmd = ""
        
        if len(lFiles) == 1:
            if self.useSymLink:
                cmd += "ln -s"
            else:
                cmd += "cp"
            cmd += " %s " % lFiles[0]
        else:
            cmd += "cat"
            for f in lFiles:
                cmd += " %s" % f
            cmd += " > "
            
        if self.outDir != "":
            cmd += "%s/" % self.outDir
        cmd += "%s.%s" % (group, self.suffix)
        
        if self.verbose > 1:
            print(cmd)
        return cmd
        
        
    def handleAllGroups(self):
        if self.verbose > 0:
            print("handle all groups ...")
            sys.stdout.flush()
            
        groups = self.group2files.keys()
        groups.sort()
        
        for i in range(0, len(groups)):
            if self.verbose == 1:
                progressBar((i+1) / float(len(groups)))
                
            group = groups[i]
            lFiles = self.group2files[group]
            
            cmd = self.handleOneGroup(group, lFiles)
            os.system(cmd)
            
            
    def run(self):
        self.loadInputFile()
        self.handleAllGroups()
        
        
if __name__ == "__main__":
    i = CatGroupedFiles()
    
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
