#!/usr/bin/env python

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
        self.input = ""
        
        
    def help(self):
        msg = "`%s' does this and that.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += " -h, --help\tdisplay the help and exit\n"
        msg += " -V, --version\toutput version information and exit\n"
        msg += " -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += " -i\tinput\n"
        msg += "\n"
        msg += "Examples:\n"
        print msg; sys.stdout.flush()
        
        
    def version(self):
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
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:i:",
                                        ["help", "version", "verbose="])
        except getopt.GetoptError, err:
            sys.stderr.write("%s\n" % str(err))
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
            elif o == "-i":
                 self.input = a
            else:
                assert False, "unhandled option"
                
                
    def checkAttributes(self):
        if self.input == "":
            msg = "ERROR: missing compulsory option -i"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.input):
            msg = "ERROR: can't find '%s'" % self.input
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
            
    def run(self):
        self.checkAttributes()
        
        if self.verbose > 0:
            startTime = time.time()
            msg = "START %s %s" % (os.path.basename(sys.argv[0]),
                                   time.strftime("%Y-%m-%d %H:%M:%S"))
            msg += "\ncmd-line: %s" % ' '.join(sys.argv)
            print msg; sys.stdout.flush()
            
        # ... specific code ...
        
        if self.verbose > 0:
            msg = "END %s %s" % (os.path.basename(sys.argv[0]),
                                 time.strftime("%Y-%m-%d %H:%M:%S"))
            endTime = time.time()
            runLength = datetime.timedelta(seconds=
                                           math.floor(endTime - startTime))
            msg += " (%s)" % str(runLength)
            print msg; sys.stdout.flush()
            
            
if __name__ == "__main__":
    i = MyClass()
    i.setAttributesFromCmdLine()
    i.run()
