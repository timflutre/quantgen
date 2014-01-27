#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: does this and that
# choose between:
# Author: Timothee Flutre
# Not copyrighted -- provided to the public domain
# or:
# Copyright (C) 2011-2013 Timothee Flutre
# License: GPLv3+

import sys
import os
import getopt
import time
import datetime
import math
import gzip


class MyClass(object):
    
    def __init__(self):
        self.verbose = 1
        self.inFile = ""
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' does this and that.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "  -i, --input\tpath to the input file\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -i <input>\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Remarks:\n"
        msg += "  This is my typical template file for python."
        print msg; sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        """
        msg = "%s 1.0\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Written by Timothee Flutre.\n"
        msg += "\n"
# choose between:
        msg += "Not copyrighted -- provided to the public domain\n"
# or:
        msg += "Copyright (C) 2011-2013 Timothee Flutre.\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "This is free software; see the source for copying conditions.  There is NO\n"
        msg += "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
        print msg; sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:i:",
                                        ["help", "version", "verbose=",
                                         "input="])
        except getopt.GetoptError, err:
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
            
            
    def run(self):
        pass
    
    
if __name__ == "__main__":
    i = MyClass()
    
    i.setAttributesFromCmdLine()
    
    i.checkAttributes()
    
    if i.verbose > 0:
        startTime = time.time()
        msg = "START %s %s" % (os.path.basename(sys.argv[0]),
                               time.strftime("%Y-%m-%d %H:%M:%S"))
        msg += "\ncmd-line: %s" % ' '.join(sys.argv)
        msg += "\ncwd: %s" % os.getcwd()
        print msg; sys.stdout.flush()
        
    i.run()
    
    if i.verbose > 0:
        msg = "END %s %s" % (os.path.basename(sys.argv[0]),
                             time.strftime("%Y-%m-%d %H:%M:%S"))
        endTime = time.time()
        runLength = datetime.timedelta(seconds=
                                       math.floor(endTime - startTime))
        msg += " (%s)" % str(runLength)
        print msg; sys.stdout.flush()
        
