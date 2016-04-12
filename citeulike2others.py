#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: download data from CiteULike for usage with Jabref or Zotero
# Copyright (C) 2014-2016 Timothée Flutre
# Persons: Timothée Flutre [cre,aut]
# License: GPL-3+
# Versioning: https://github.com/timflutre/quantgen

# http://wiki.citeulike.org/index.php/Importing_and_Exporting

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
from pybtex.database.input import bibtex as bib_i
from pybtex.database.output import bibtex as bib_o
# import bibtexparser
import RISparser # from https://github.com/sacabuche/RISparser

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n" % msg)
        sys.exit(1)
        
progVersion = "1.3.3" # http://semver.org/


def user_input(msg):
    if sys.version_info[0] == 2:
        return raw_input(msg)
    elif sys.version_info[0] == 3:
        return input(msg)
    else:
        msg = "ERROR: Python's major version should be 2 or 3"
        sys.stderr.write("%s\n" % msg)
        sys.exit(1)
        
        
class Citeulike2Others(object):
    
    def __init__(self):
        self.verbose = 1
        self.tasks = []
        self.identifier = "timflutre"
        self.email = "timflutre@gmail.com"
        self.cookieFile = ""
        self.outFormat = "json"
        self.jsonFile = ""
        self.jsonRefs = None
        self.libDir = ""
        self.bibtexFile = ""
        self.bibtexRefs = None
        self.risFile = ""
        self.risRefs = None
        self.otherTool = "zotero"
        self.tag = None
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' downloads data from CiteULike for usage with Jabref or Zotero.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "  -t, --task\ttask(s) to perform (can be '1+2' for instance)\n"
        msg += "\t\t1: save cookies (requires -i, -e)\n"
        msg += "\t\t2: download file with references (requires -i, -e, -f; can use -a)\n"
        msg += "\t\t3: download new attached files (requires -i, -e, -j, -p)\n"
        msg += "\t\t4: adapt the reference file for another tool (requires -j, -p, -b/-r, -o)\n"
        msg += "\t\tadd link(s) to file(s)\n"
        msg += "\t\te.g. rename 'comment' into 'annote' in Bibtex for Zotero\n"
        msg += "\t\tor turn note formatting into HTML in RIS for Zotero\n"
        msg += "  -i, --id\tyour CiteULike identifier (default=timflutre)\n"
        msg += "  -e, --email\tyour email (default=timflutre@gmail.com)\n"
        msg += "  -f, --fmt\tformat of the file to be downloaded\n"
        msg += "\t\tdefault=json/bibtex/ris\n"
        msg += "  -j, --json\tpath to the JSON file\n"
        msg += "  -p, --path\tpath to the directory with all the files\n"
        msg += "  -b, --bibtex\tpath to the Bibtex file\n"
        msg += "  -r, --ris\tpath to the RIS file\n"
        msg += "  -o, --other\tother tool (default=zotero/jabref)\n"
        msg += "  -a, --tag\tspecific CiteULike tag to retrieve (whole library otherwise)\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -t 1\n" % os.path.basename(sys.argv[0])
        msg += "  %s -t 2 -f json\n" % os.path.basename(sys.argv[0])
        msg += "  %s -t 3 -j refs.json -p ~/work/biblio\n" % os.path.basename(sys.argv[0])
        msg += "  %s -t 2 -f bibtex\n" % os.path.basename(sys.argv[0])
        msg += "  %s -t 4 -j refs.json -p ~/work/biblio -b refs.bib\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timflutre@gmail.com>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        
        The person roles complies with R's guidelines (The R Journal Vol. 4/1, June 2012).
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2014-2016 Timothée Flutre.\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt(sys.argv[1:],
                                       "hVv:t:i:e:f:j:p:b:r:o:a:",
                                       ["help", "version", "verbose=",
                                        "tasks=", "id=", "email=", "format=",
                                        "json=", "path=", "bibtex=", "ris=",
                                        "out=", "tag="])
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
            elif o == "-f" or o == "--format":
                self.outFormat = a
            elif o == "-j" or o == "--json":
                self.jsonFile = a
            elif o == "-p" or o == "--path":
                self.libDir = a
            elif o == "-b" or o == "--bibtex":
                self.bibtexFile = a
            elif o == "-r" or o == "--ris":
                self.risFile = a
            elif o == "-o" or o == "--other":
                self.otherTool = a
            elif o == "-a" or o == "--tag":
                self.tag = a
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
            if task not in ["1","2","3","4","5","6"]:
                msg = "ERROR: -t %s is not recognized" % task
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
        if "1" in self.tasks and (self.identifier == "" or self.email == ""):
            msg = "ERROR: missing compulsory option(s) -i or -e with -t 1"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "2" in self.tasks and (self.identifier == "" or self.email == "" \
                                  or self.outFormat == ""):
            msg = "ERROR: missing compulsory option(s) -i or -e or -f with -t 2"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "2" in self.tasks and self.outFormat not in ["json", "bibtex", "ris"]:
            msg = "ERROR: unknown -f %s with -t 2" % self.outFormat
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "3" in self.tasks and self.jsonFile == "":
            msg = "ERROR: missing compulsory option -j with -t 3"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.jsonFile != "" and not os.path.exists(self.jsonFile):
            msg = "ERROR: can't find file %s" % self.jsonFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "3" in self.tasks and self.libDir == "":
            msg = "ERROR: missing compulsory option -p with -t 3"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.libDir != "" and not os.path.exists(self.libDir):
            msg = "ERROR: can't find directory %s" % self.libDir
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "4" in self.tasks and self.bibtexFile == "" and self.risFile == "":
            msg = "ERROR: missing compulsory option -b or -r with -t 4"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.bibtexFile != "" and not os.path.exists(self.bibtexFile):
            msg = "ERROR: can't find file %s" % self.bibtexFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.risFile != "" and not os.path.exists(self.risFile):
            msg = "ERROR: can't find file %s" % self.risFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "4" in self.tasks and self.libDir == "":
            msg = "ERROR: missing compulsory option -p with -t 3"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "4" in self.tasks and self.otherTool not in ["jabref", "zotero"]:
            msg = "ERROR: missing compulsory option -o with -t 4"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
        if "1" in self.tasks or "2" in self.tasks or "3" in self.tasks or \
           "4" in self.tasks:
            self.cookieFile = "cookies_citeulike_%s.txt" % self.identifier
            
            
    def saveCookies(self):
        if self.verbose > 0:
            print("save cookies ...")
            sys.stdout.flush()
            
        if os.path.exists(self.cookieFile):
            os.remove(self.cookieFile)
            
        pwd = getpass.getpass()
        
        cmd = "wget"
        cmd += " --header=\"User-Agent: %s" % self.identifier
        cmd += "/%s" % self.email
        cmd += " downloader/%s\"" % progVersion
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
            print("cookies saved in file '%s'" % self.cookieFile)
            print("CHECK IT! (to see if it's not empty)")
            
            
    def downloadFile(self):
        if self.verbose > 0:
            print("download %s file ..." % self.outFormat)
            sys.stdout.flush()
            
        if not os.path.exists(self.cookieFile):
            msg = "ERROR: missing cookie file '%s'" % self.cookieFile
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
        if self.verbose > 0:
            print("cookie file: %s" % self.cookieFile)
            
        outFile = "citeulike_%s" % self.identifier
        if self.tag:
            outFile += "_%s" % self.tag
        if self.outFormat == "bibtex":
            fileExt = "bib"
        else:
            fileExt = self.outFormat
        outFile += ".%s" % fileExt
        
        cmd = "wget"
        cmd += " --header=\"User-Agent: %s" % self.identifier
        cmd += "/%s" % self.email
        cmd += " downloader/%s\"" % progVersion
        cmd += " --load-cookies %s" % self.cookieFile
        cmd += " --quiet"
        cmd += " -O %s" % outFile
        cmd += " http://www.citeulike.org/%s" % self.outFormat
        cmd += "/user/%s" % self.identifier
        if self.tag:
            cmd += "/tag/%s" % self.tag
            
        if self.verbose > 1:
            print(cmd)
            
        os.system(cmd)
        
        if self.verbose > 0:
            print("library saved in file '%s'" % outFile)
            
            
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
            
        if not os.path.exists(self.cookieFile):
            msg = "ERROR: missing cookie file '%s'" % self.cookieFile
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
        if self.verbose > 0:
            print("cookie file: %s" % self.cookieFile)
            
        totalNbFiles = 0
        for entry in self.jsonRefs:
            if "userfiles" in entry:
                totalNbFiles += len(entry["userfiles"])
                for f in entry["userfiles"]:
                    p = "%s/%s" % (self.libDir, f["name"])
                    root, ext = os.path.splitext(p)
                    if not os.path.exists(p):
                        print("download '%s'" % entry["title"])
                        cmd = "wget"
                        cmd += " --header=\"User-Agent: %s" % self.identifier
                        cmd += "/%s" % self.email
                        cmd += " downloader/%s\"" % progVersion
                        cmd += " --load-cookies %s" % self.cookieFile
                        cmd += " --quiet"
                        cmd += " -O %s" % p
                        cmd += " http://www.citeulike.org/%s" % f["path"]
                        os.system(cmd)
                        
        if self.verbose > 0:
            print("total nb of files in JSON: %i" % totalNbFiles)
            sys.stdout.flush()
            
            
    def rmvOldFiles(self):
        if self.verbose > 0:
            print("remove the old files ...")
            sys.stdout.flush()
            
        lOldFiles = glob.glob("%s/*" % self.libDir)
        if self.verbose > 0:
            print("total nb of files in libDir: %i" % len(lOldFiles))
            
        sNewFiles = set()
        for entry in self.jsonRefs:
            if "userfiles" in entry:
                for f in entry["userfiles"]:
                    p = "%s/%s" % (self.libDir, f["name"])
                    sNewFiles.add(p)
        if self.verbose > 0:
            print("total nb of files in JSON: %i" % len(sNewFiles))
        if len(sNewFiles) == 0:
            msg = "ERROR: no new files, possibly due to incomplete JSON file?"
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
            
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
                
            self.bibtexRefs = bib_i.Parser().parse_file(self.bibtexFile)
            
            # for bibtexparser:
            # bibtexHandle = open(self.bibtexFile)
            # bibtexRefs = bibtexparser.load(bibtexHandle)
            # bibtexHandle.close()
            
            if self.verbose > 0:
                print("Bibtex file: %i entries" % len(self.bibtexRefs.entries))
                
                
    def loadRisFile(self):
        if self.risRefs == None:
            if self.verbose > 0:
                print("load RIS file ...")
                sys.stdout.flush()
                
            tmp = None
            with open(self.risFile, "r") as risHandle:
                tmp = list(RISparser.readris(risHandle))
                
            self.risRefs = {}
            for i in range(len(tmp)):
                entry = tmp[i]
                if "id" not in entry:
                    msg = "ERROR: no 'id' field in the entry #%i of %s" \
                          % (i, self.risFile)
                    sys.stderr.write("%s\n" % msg)
                    sys.exit(1)
                self.risRefs[entry["id"]] = entry
                
            if self.verbose > 0:
                print("RIS file: %i entries" % len(self.risRefs))
                
                
    def getSharedCitationKey(self, entryJson, nonJsonEntries):
        """
        Among possibly several citation keys in the JSON entry, return the one
        which is alsos in the non-JSON file (Bibtex or RIS).
        """
        ck = None
        
        if "citation_keys" not in entryJson \
           or len(entryJson["citation_keys"]) == 0:
            msg = "ERROR: '%s' has no citation_keys" % entryJson["title"]
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
            
        ck_idx = 0
        while ck_idx < len(entryJson["citation_keys"]):
            if entryJson["citation_keys"][ck_idx] in nonJsonEntries:
                break
            ck_idx += 1
        if ck_idx >= len(entryJson["citation_keys"]):
            msg = "ERROR: can't find entry in non-JSON file for '%s'" % \
                  entryJson["title"]
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
            
        ck = entryJson["citation_keys"][ck_idx]
        return ck
        
        
    def getFileExtension(self, fName):
        return os.path.splitext(fName)[1][1:].strip().lower()
        
        
    def getMainPdf(self, entryJson):
        fMain = None
        fMainExt = None
        for idx, f in enumerate(entryJson["userfiles"]):
            fExt = self.getFileExtension(f["name"])
            if fExt == "pdf":
                fMain = f
                fMainExt = fExt
                break
        return fMain, fMainExt
        
        
    def formatFileForBibtex(self, f, fExt):
        if self.otherTool == "jabref":
            return "%s:%s/%s:%s" % (f["name"],
                                    self.libDir,
                                    f["name"],
                                    fExt.upper())
        elif self.otherTool == "zotero":
            return "%s:%s/%s:application/%s" % (f["name"],
                                                self.libDir,
                                                f["name"],
                                                fExt)
            
            
    def formatFilesForBibtex(self, entryJson):
        txt = None
        
        tmp = []
        fMain = None
        if len(entryJson["userfiles"]) > 1:
            fMain, fMainExt = self.getMainPdf(entryJson)
            if fMain != None:
                tmp.append(self.formatFileForBibtex(fMain, fMainExt))
                
        for f in entryJson["userfiles"]:
            if fMain != None and f["name"] == fMain["name"]:
                continue
            fExt = self.getFileExtension(f["name"])
            # if fExt not in ["pdf", "png", "jpg", "jpeg", "gif"]:
            #     msg = "ERROR: unknown file extension '%s'" % fExt
            #     msg += "\n%s" % entryJson["title"]
            #     sys.stderr.write("%s\n" % msg)
            #     sys.exit(1)
            tmp.append(self.formatFileForBibtex(f, fExt))
            
        if len(tmp) == 1:
            txt = tmp[0]
        else:
            txt = ";".join(tmp)
        if self.otherTool == "zotero":
            txt = "{%s}" % txt
            
        return txt
        
        
    def addFileFieldToBibtex(self):
        """
        To each Bibtex entry, add a field 'file' if the Json entry has any.
        """
        if self.verbose > 0:
            print("add 'file' field to Bibtex entries ...")
            sys.stdout.flush()
            
        for entryJson in self.jsonRefs:
            if "userfiles" not in entryJson or \
               len(entryJson["userfiles"]) == 0:
                continue
            ck = self.getSharedCitationKey(entryJson, self.bibtexRefs.entries)
            self.bibtexRefs.entries[ck].fields["file"] \
                = self.formatFilesForBibtex(entryJson)
            
            
    def formatFilesForRis(self, entryJson):
        out = []
        
        fMain = None
        if len(entryJson["userfiles"]) > 1:
            fMain, fMainExt = self.getMainPdf(entryJson)
            if fMain != None:
                out.append("%s/%s" % (self.libDir, fMain["name"]))
                
        for f in entryJson["userfiles"]:
            if fMain != None and f["name"] == fMain["name"]:
                continue
            out.append("%s/%s" % (self.libDir, f["name"]))
            
        return out
    
    
    def addFileFieldToRis(self):
        """
        To each RIS entry, add a field 'L1' if the Json entry has any.
        """
        if self.verbose > 0:
            print("add 'L1' field to RIS entries ...")
            sys.stdout.flush()
            
        for entryJson in self.jsonRefs:
            if "userfiles" not in entryJson or \
               len(entryJson["userfiles"]) == 0:
                continue
            ck = self.getSharedCitationKey(entryJson, self.risRefs)
            self.risRefs[ck]["L1"] \
                = self.formatFilesForRis(entryJson)
            
            
    def editCommentField(self):
        """
        Zotero assumes notes to be in field 'annote'.
        """
        if self.verbose > 0:
            print("edit 'comment' field of Bibtex entries ...")
            sys.stdout.flush()
            
        for ck in self.bibtexRefs.entries:
            entryBibtex = self.bibtexRefs.entries[ck]
            if "comment" in entryBibtex.fields.keys():
                self.bibtexRefs.entries[ck].fields["annote"] \
                    = entryBibtex.fields["comment"].replace("* ", "\n* ")
                del self.bibtexRefs.entries[ck].fields["comment"]
                
                
    def editKeywordsField(self):
        """
        Zotero keeps "\" before "_".
        """
        if self.verbose > 0:
            print("edit 'keywords' field of Bibtex entries ...")
            sys.stdout.flush()
            
        for ck in self.bibtexRefs.entries:
            entryBibtex = self.bibtexRefs.entries[ck]
            if "keywords" in entryBibtex.fields.keys():
                self.bibtexRefs.entries[ck].fields["keywords"] \
                    = entryBibtex.fields["keywords"].replace("\\", "")
                
                
    def writeBibtexFile(self):
        if self.verbose > 0:
            print("write new Bibtex file ...")
            sys.stdout.flush()
            
        newBibtexFile = "%s_for-%s.bib" % \
                        (os.path.splitext(self.bibtexFile)[0],
                         self.otherTool)
        if os.path.exists(newBibtexFile):
            os.remove(newBibtexFile)
            
        bibtexRefsWriter = bib_o.Writer().write_file(self.bibtexRefs,
                                                     newBibtexFile)
        # for bibtexparser:
        # newBibtexHandle = open(self.newBibtexFile, "w")
        # bibtexparser.dump(bibtexRefs, newBibtexHandle)
        # newBibtexHandle.close()
        
        if self.otherTool == "zotero":
            cmd = "sed -i 's/file = \\\"{/file = {/g' %s" % newBibtexFile
            os.system(cmd)
            
        if self.verbose > 0:
            print("library saved in file '%s'" % newBibtexFile)
            
            
    def writeRisFile(self):
        if self.verbose > 0:
            print("write new RIS file ...")
            sys.stdout.flush()
            
        newRisFile = "%s_for-%s.ris" % \
                        (os.path.splitext(self.risFile)[0],
                         self.otherTool)
        if os.path.exists(newRisFile):
            os.remove(newRisFile)
            
        newRisHandle = open(newRisFile, "w")
        for entryRis in self.risRefs:
            continue # TODO
        newRisHandle.close()
        
        if self.verbose > 0:
            print("library saved in file '%s'" % newRisFile)
            
            
    def run(self):
        if "1" in self.tasks:
            self.saveCookies()
            
        if "2" in self.tasks:
            self.downloadFile()
            
        if "3" in self.tasks:
            self.loadJsonFile()
            self.downloadNewFiles()
            self.rmvOldFiles()
            
        if "4" in self.tasks:
            self.loadJsonFile()
            if self.bibtexFile != "":
                self.loadBibtexFile()
                self.addFileFieldToBibtex()
                if self.otherTool == "zotero":
                    self.editKeywordsField()
                self.writeBibtexFile()
            elif self.risFile != "":
                self.loadRisFile()
                self.addFileFieldToRis()
                # if self.otherTool == "zotero":
                #     self.editNoteField()
                self.writeRisFile()
                
                
if __name__ == "__main__":
    i = Citeulike2Others()
    
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
