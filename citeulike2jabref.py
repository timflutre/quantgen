#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: download data from CiteULike for usage with Jabref
# Copyright (C) 2014 Timothée Flutre
# Author: Timothée Flutre
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
from subprocess import Popen, PIPE
import math
import gzip
import getpass
import json
import glob
from pybtex.database.input import bibtex
# import bibtexparser

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
def user_input(msg):
    if sys.version_info[0] == 2:
        return raw_input(msg)
    elif sys.version_info[0] == 3:
        return input(msg)
    else:
        msg = "ERROR: Python's major version should be 2 or 3"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
        
class Citeulike2Jabref(object):
    
    def __init__(self):
        self.verbose = 1
        self.tasks = []
        self.identifier = "timflutre"
        self.email = "timflutre@gmail.com"
        self.cookieFile = "cookies_citeulike.txt"
        self.jsonFile = ""
        self.jsonRefs = None
        self.libDir = ""
        self.bibtexFile = ""
        self.bibtexRefs = None
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' downloads data from CiteULike for usage with Jabref.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "  -t, --task\ttask(s) to perform (can be '1+2' for instance)\n"
        msg += "\t\t1: save cookies (requires -i, -e)\n"
        msg += "\t\t2: download JSON file (requires -i, -e)\n"
        msg += "\t\t3: download Bibtex file (requires -i, -e)\n"
        msg += "\t\t4: download new files (requires -i, -e, -j, -p)\n"
        msg += "\t\t5: add 'file' field to Bibtex file (requires -j, -b)\n"
        msg += "  -i, --id\tyour CiteULike identifier (default=timflutre)\n"
        msg += "  -e, --email\tyour email (default=timflutre@gmail.com)\n"
        msg += "  -j, --json\tpath to the JSON file\n"
        msg += "  -p, --path\tpath to the directory with all the files\n"
        msg += "  -b, --bibtex\tpath to the Bibtex file\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -t 1+2+3\n" % os.path.basename(sys.argv[0])
        msg += "  %s -t 4+5 -j refs.json -p ~/work/biblio -b refs.bib\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timflutre@gmail.com>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        """
        msg = "%s 1.0\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Copyright (C) 2014 Timothée Flutre.\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "This is free software; see the source for copying conditions.  There is NO\n"
        msg += "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
        msg += "\n"
        msg += "Written by Timothée Flutre."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:t:i:e:j:b:p:",
                                        ["help", "version", "verbose=",
                                         "tasks=", "id=", "email=", "json=",
                                         "bibtex=", "path="])
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
            elif o == "-t" or o == "--tasks":
                self.tasks = a.split("+")
            elif o == "-i" or o == "--id":
                self.identifier = a
            elif o == "-e" or o == "--email":
                self.email = a
            elif o == "-j" or o == "--json":
                self.jsonFile = a
            elif o == "-b" or o == "--bibtex":
                self.bibtexFile = a
            elif o == "-p" or o == "--path":
                self.libDir = a
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if self.tasks == []:
            msg = "ERROR: missing compulsory option -t"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        for task in self.tasks:
            if task not in ["1","2","3","4","5"]:
                msg = "ERROR: -t %s is not recognized" % task
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        if "1" in self.tasks and (self.identifier == "" or self.email == ""):
            msg = "ERROR: missing compulsory option(s) -i or -e with -t 1"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if ("2" in self.tasks or "3" in self.tasks) and \
           (self.identifier == "" or self.email == ""):
            msg = "ERROR: missing compulsory option(s) -i or -e with -t 2 or -t 3"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "4" in self.tasks and self.jsonFile == "":
            msg = "ERROR: missing compulsory option -j with -t 4"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.jsonFile != "" and not os.path.exists(self.jsonFile):
            msg = "ERROR: can't find file %s" % self.jsonFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "4" in self.tasks and self.libDir == "":
            msg = "ERROR: missing compulsory option -p with -t 4"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.libDir != "" and not os.path.exists(self.libDir):
            msg = "ERROR: can't find directory %s" % self.libDir
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "5" in self.tasks and self.bibtexFile == "":
            msg = "ERROR: missing compulsory option -b with -t 5"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.bibtexFile != "" and not os.path.exists(self.bibtexFile):
            msg = "ERROR: can't find file %s" % self.bibtexFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
            
    def saveCookies(self):
        if self.verbose > 0:
            print("save cookies ...")
            sys.stdout.flush()
            
        pwd = getpass.getpass()
        
        cmd = "wget"
        cmd += " --header=\"User-Agent: %s" % self.identifier
        cmd += "/%s" % self.email
        cmd += " downloader/1.0\""
        cmd += " -O /dev/null"
        cmd += " --quiet"
        cmd += " --keep-session-cookies"
        cmd += " --save-cookies %s" % self.cookieFile
        cmd += " --post-data=\"username=%s" % self.identifier
        cmd += "&password=%s" % pwd
        cmd += "&perm=1\""
        cmd += " http://www.citeulike.org/login.do"
        if self.verbose > 1:
            print(cmd)
            
        os.system(cmd)
        
        if self.verbose > 0:
            print("cookies saved in file %s" % self.cookieFile)
            
            
    def downloadFile(self, fileFormat):
        if self.verbose > 0:
            print("download %s file ..." % fileFormat)
            sys.stdout.flush()
            
        cmd = "wget"
        cmd += " --header=\"User-Agent: %s" % self.identifier
        cmd += "/%s" % self.email
        cmd += " downloader/1.0\""
        cmd += " --load-cookies %s" % self.cookieFile
        cmd += " --quiet"
        cmd += " -O citeulike_%s.%s" % (self.identifier, fileFormat)
        cmd += " http://www.citeulike.org/%s" % fileFormat
        cmd += "/user/%s" % self.identifier
        
        if self.verbose > 1:
            print(cmd)
            
        os.system(cmd)
        
        if self.verbose > 0:
            print("library saved in file citeulike_%s.%s" \
                  % (self.identifier, fileFormat))
            
            
    def loadJsonFile(self):
        if self.jsonRefs == None:
            if self.verbose > 0:
                print("load JSON file ...")
                sys.stdout.flush()
            jsonHandle = open(self.jsonFile)
            self.jsonRefs = json.load(jsonHandle)
            jsonHandle.close()
            if self.verbose > 0:
                print("JSON file: %i entries" % len(self.jsonRefs))
                
                
    def downloadNewFiles(self):
        if self.verbose > 0:
            print("download the new files ...")
            sys.stdout.flush()
        totalNbFiles = 0
        for entry in self.jsonRefs:
            if "userfiles" in entry:
                totalNbFiles += len(entry["userfiles"])
                for f in entry["userfiles"]:
                    p = "%s/%s" % (self.libDir, f["name"])
                    root, ext = os.path.splitext(p)
                    if not os.path.exists(p):
                        print(entry["title"])
                        cmd = "wget"
                        cmd += " --header=\"User-Agent: %s" % self.identifier
                        cmd += "/%s" % self.email
                        cmd += " downloader/1.0\""
                        cmd += " --load-cookies %s" % self.cookieFile
                        cmd += " --quiet"
                        cmd += " -O %s" % p
                        cmd += " http://www.citeulike.org/%s" % f["path"]
                        os.system(cmd)
        if self.verbose > 0:
            print("total nb of files: %i" % totalNbFiles)
            sys.stdout.flush()
            
            
    def rmvOldFiles(self):
        if self.verbose > 0:
            print("remove the old files ...")
            sys.stdout.flush()
        lOldFiles = glob.glob("%s/*" % self.libDir)
        sNewFiles = set()
        for entry in self.jsonRefs:
            if "userfiles" in entry:
                for f in entry["userfiles"]:
                    p = "%s/%s" % (self.libDir, f["name"])
                    sNewFiles.add(p)
        for f in lOldFiles:
            if f not in sNewFiles:
                wantRmvFile = user_input("Do you want to remove the file %s? [y/n] " % f)
                if wantRmvFile.lower() == "y" or wantRmvFile.lower() == "yes":
                    os.remove(f)
                    print("=> done!")
                    sys.stdout.flush()
        for f in sNewFiles:
            if f not in lOldFiles:
                print(f)
                
                
    def loadBibtexFile(self):
        if self.bibtexRefs == None:
            if self.verbose > 0:
                print("load Bibtex file ...")
                sys.stdout.flush()
            self.bibtexRefs = bibtex.Parser().parse_file(self.bibtexFile)
            # bibtexHandle = open(self.bibtexFile)
            # bibtexRefs = bibtexparser.load(bibtexHandle)
            # bibtexHandle.close()
            if self.verbose > 0:
                print("Bibtex file: %i entries" % len(self.bibtexRefs.entries))


    def addFileFieldToBibtex(self):
        # file = {:~/work/refs/main.pdf:PDF;:~/work/refs/supp.pdf:PDF},
        if self.verbose > 0:
            print("add 'file' field to Bibtex entries ...")
            sys.stdout.flush()
        for entry in self.jsonRefs:
            if "citation_keys" not in entry \
               or len(entry["citation_keys"]) == 0:
                msg = "ERROR: '%s' has no citation_keys" % entry["title"]
                sys.stderr.write("%s\n\n" % msg)
                sys.exit(1)
            ck_idx = 0
            while ck_idx < len(entry["citation_keys"]):
                if entry["citation_keys"][ck_idx] in self.bibtexRefs.entries:
                    break
                ck_idx += 1
            if ck_idx >= len(entry["citation_keys"]):
                msg = "ERROR: can't find entry '%s' in Bibtex file" % \
                      entry["citation_keys"][ck_idx]
                msg += "\n%s" % entry["title"]
                sys.stderr.write("%s\n\n" % msg)
                sys.exit(1)
                continue
            self.bibtexRefs.entries[entry["citation_keys"][ck_idx]].fields["file"] = []
            
            
    def writeBibtexFile(self):
        if self.verbose > 0:
            print("write new Bibtex file ...")
            sys.stdout.flush()
        newBibtexFile = "%s_new" % self.bibtexFile
        # newBibtexHandle = open(self.newBibtexFile, "w")
        # bibtexparser.dump(bibtexRefs, newBibtexHandle)
        # newBibtexHandle.close()
        # TODO: mv newfile oldfile
        
        
    def run(self):
        if "1" in self.tasks:
            self.saveCookies()
        if "2" in self.tasks:
            self.downloadFile("json")
        if "3" in self.tasks:
            self.downloadFile("bibtex")
        if "4" in self.tasks:
            self.loadJsonFile()
            self.downloadNewFiles()
            self.rmvOldFiles()
        if "5" in self.tasks:
            self.loadJsonFile()
            self.loadBibtexFile()
            self.addFileFieldToBibtex()
            self.writeBibtexFile()
            
            
if __name__ == "__main__":
    i = Citeulike2Jabref()
    
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
