#!/usr/bin/env python

# Author: Tim Flutre
# Aim: retrieve data from the Gene Expression Omnibus at the NCBI

import urllib2
import lxml.html
import sys
import os
import getopt
import shutil
import time

BASEURL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"


class Sample( object ):
    def __init__(self, a="", d="", u=""):
        self.acc = a
        self.desc = d
        self.url_data = u
    def getUrl(self):
        if self.acc == "":
            sys.stderr.write("error: empty URL\n")
            sys.exit(1)
        return "%sacc=%s" % (BASEURL, self.acc)
    def __repr__(self):
        txt = ""
        if self.acc != "":
            txt += "%s" % self.acc
        if self.desc != "":
            txt += " %s" % self.desc
        if self.url_data != "":
            txt += " %s" % self.url_data
        return txt


class GetGeoData(object):
    
    def __init__(self):
        self.study = ""
        self.dSamples = {}
        self.verbose = 0


    def help(self):
        print
        print "usage: %s.py [options]" % self.__class__.__name__
        print "options:"
        print "     -h: this help"
        print "     -s: study (e.g. GSE17080)"
        print "     -v: verbosity level (default=0/1)"
        print


    def setAttributesFromCmdLine(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:],"hs:v:")
        except getopt.GetoptError, err:
            sys.stderr.write("%s\n" % str(err))
            self.help()
            sys.exit(1)
        for o,a in opts:
            if o == "-h":
                self.help()
                sys.exit(0)
            elif o == "-s":
                self.study = a
            elif o == "-v":
                self.verbose = int(a)
        if self.study == "":
            msg = "error: missing study (-s)"
            sys.stderr.write("%s\n" % msg)
            self.help()
            sys.exit(1)

        
    def getSourceCodeFromWebPage(self, url):
        """
        Return the source code of a web page.
        @param url: URL to retrieve
        @type url: string
        @return: source code corresponding to the given URL
        @rtype: string
        """
        try:
            result = urllib2.urlopen( url )
        except urllib2.HTTPError:
            if self.verbose > 1:
                sys.stderr.write("HTTPError\n")
            sys.exit(1)
        except urllib2.URLError:
            if self.verbose > 1:
                sys.stderr.write("URLError\n")
            sys.exit(1)
        return result.read()


    def getSampleList(self, outFile):
        """
        Retrieve the list of all samples in the study.
        @note: load this list into the attribute 'self.dSamples'
        """
        if self.verbose > 0:
            print "get sample list"; sys.stdout.flush()
        url = "%sacc=%s" % (BASEURL, self.study)
        html = self.getSourceCodeFromWebPage(url)
        root = lxml.html.fromstring(html)
        #tds = root.cssselect("td")
        i = 0
        previ = 0
        prevSample = ""
        for elem in root.iter():
            i += 1
            if elem.tag == "a" and "acc.cgi?acc=GSM" in elem.attrib["href"]:
                #print elem.attrib["href"].split("=")[1]
                iSample = Sample( elem.attrib["href"].split("=")[1] )
                self.dSamples[ iSample.acc ] = iSample
                previ = i
                prevSample = iSample.acc
            elif i == previ + 1 and elem.tag == "td":
                #print elem.text
                self.dSamples[prevSample].desc = elem.text
#        print len(self.dSamples.keys())
#         lSamples = self.dSamples.keys()
#         lSamples.sort()
#         for sample in lSamples:
#             print self.dSamples[sample]
        if self.verbose > 0:
            print "%i samples listed" % len(self.dSamples.keys())
            sys.stdout.flush()


    def saveSamples(self, outFile):
        if self.verbose > 0:
            print "save samples in file '%s'" % outFile
            sys.stdout.flush()
        outH = open(outFile, "w")
        outH.write("acc\tdesc\turl\n")
        lAccessions = self.dSamples.keys()
        lAccessions.sort()
        for acc in lAccessions:
            iSample = self.dSamples[acc]
            txt = "%s" % iSample.acc
            txt += "\t%s" % iSample.desc
            txt += "\t%s" % iSample.url_data
            outH.write("%s\n" % txt)
        outH.close()
        if self.verbose > 0:
            print "done"
            sys.stdout.flush()


    def loadSamples(self, inFile):
        if self.verbose > 0:
            print "load samples from file '%s'" % inFile
            sys.stdout.flush()
        inH = open(inFile)
        line = inH.readline()
        while True:
            line = inH.readline()
            if line == "":
                break
            tok = line.replace("\n","").split("\t")
            iSample = Sample(tok[0], tok[1])
            if len(tok) > 2:
                iSample.url_data = tok[2]
            self.dSamples[ tok[0] ] = iSample
        inH.close()
        if self.verbose > 0:
            print "%i samples loaded" % len(self.dSamples.keys())
            sys.stdout.flush()


    def getProbeIntensitiesUrl(self):
        if self.verbose > 0:
            print "get URLs for probe intensities"
            sys.stdout.flush()
        lAccessions = self.dSamples.keys()
        lAccessions.sort()
        nbSamples = len(lAccessions)
        for i in range(0,nbSamples):
            acc = lAccessions[i]
            iSample = self.dSamples[acc]
            if iSample.acc != acc:
                sys.stderr.write("error\n")
                sys.exit(1)
            if self.verbose > 0:
                print "sample #%i/%i: %s" % (i+1, nbSamples, iSample.acc)
                sys.stdout.flush()
            html = self.getSourceCodeFromWebPage(iSample.getUrl())
            root = lxml.html.fromstring(html)
            inputs = root.cssselect("input")
            for j in inputs:
                txt = lxml.html.tostring(j)
                if "fulltable" in txt:
                    partial_url = txt.split("acc.cgi?")[1].split("', '")[0]
                    iSample.url_data = "%s%s" % (BASEURL, partial_url)
        if self.verbose > 0:
            print "%i URLs retrieved" % (i+1)
            sys.stdout.flush()


    def saveProbeIntensities(self):
        newDir = "%s_all_samples" % self.study
        if os.path.exists(newDir):
            shutil.rmtree(newDir)
        os.mkdir(newDir)
        os.chdir(newDir)
        lAccessions = self.dSamples.keys()
        lAccessions.sort()
        nbSamples = len(lAccessions)
        for i in range(0,nbSamples):
            acc = lAccessions[i]
            iSample = self.dSamples[acc]
            outFile = "%s_sample_%s_%s.txt" % (self.study, iSample.acc, iSample.desc)
            if self.verbose > 0:
                print "save data into file '%s'" % outFile
                sys.stdout.flush()
            outH = open(outFile, "w")
            outH.write("probe\tvalue\n")
            html = self.getSourceCodeFromWebPage(iSample.url_data)
            root = lxml.html.fromstring(html)
            pres = root.cssselect("pre")
            txt = lxml.html.tostring(pres[0])
            for t in txt.split("\n"):
                tok = t.split("\t")
                if len(tok) == 2:
                    try:
                        float(tok[1])
                    except ValueError:
                        continue
                    outH.write("%s\n" % t)
            outH.close()
        os.chdir("..")


    def run(self):
        if self.verbose > 0:
            print "begin %s.py (%s)" % (self.__class__.__name__, time.strftime("%Y-%m-%d %H:%M:%S"))
            sys.stdout.flush()

        outFile = "%s_acc_desc_url.txt" % self.study
        if not os.path.exists(outFile):
            self.getSampleList(outFile)
            self.getProbeIntensitiesUrl()
            self.saveSamples(outFile)
        else:
            self.dSample = self.loadSamples(outFile)

        self.saveProbeIntensities()

        if self.verbose > 0:
            print "end %s.py (%s)" % (self.__class__.__name__, time.strftime("%Y-%m-%d %H:%M:%S"))
            sys.stdout.flush()


if __name__ == "__main__":
    i = GetGeoData()
    i.setAttributesFromCmdLine()
    i.run()
