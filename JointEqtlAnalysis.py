#!/usr/bin/env python

# Author: Timothee Flutre
# License: GPL-3
# Aim: implement the Bayesian meta-analysis model of Wen and Stephens
# and perform permutations in order to obtain feature-level P-values
# and thereby estimate the FDR of the joint analysis
# help2man -o JointEqtlAnalysis.man ./JointEqtlAnalysis.py
# groff -mandoc JointEqtlAnalysis.man > JointEqtlAnalysis.ps


import sys
import os
import getopt
import time
import datetime
import math
import glob
import scipy
import gzip
import shutil
import multiprocessing


class JointEqtlAnalysis(object):
    
    def __init__(self):
        self.verbose = 1
        self.step = 0
        self.nbPermutations = 10000
        self.seed = 1859
        self.pathToPhenoDir = ""
        self.pathToGenoFile = ""
        self.pathToFtrCoords = ""
        self.pathToFtrsSnpsLinks = ""
        self.pathToFtrsToKeep = ""
        self.pathToGridFile = ""
        
        self.originalAnalysisDir = "original_analysis"
        self.originalAbfsFile = "abfs_meta.txt"
        self.permDir = "permutations"
        self.lPathToPhenoFiles = []
        self.nbIndividuals = -1
        self.outFile = "ftr-level_permutation_pvalues.txt"
        
        
    def help(self):
        msg = "`%s' implements the Bayesian meta-analysis model of Wen and Stephens\n" % os.path.basename(sys.argv[0])
        msg += "and performs permutations in order to obtain feature-level P-values\n"
        msg += "and thereby estimate the FDR of the joint analysis.\n"
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += " -h, --help\tdisplay the help and exit\n"
        msg += " -V, --version\toutput version information and exit\n"
        msg += " -v, --verbose\tverbosity level (default=1)\n"
        msg += " -s, --step\tstep of the pipeline (1/2/3/4/5)\n"
        msg += "\t\t1: give commands to launch full analysis on the original data\n"
        msg += "\t\t2: write all permuted phenotypes (one file per subgroup)\n"
        msg += "\t\t3: calculate the summary statistics for each permutation and subgroup\n"
        msg += "\t\t4: calculate the ABF_meta for each permutation\n"
        msg += "\t\t5: calculate the gene-level joint-analysis permutation P-values\n"
        msg += " -g, --geno\tfile with genotypes in IMPUTE format\n"
        msg += "\t\ta header line with sample names is required\n"
        msg += "\t\tsamples in columns should be in same order as in phenotype file\n"
        msg += " -p, --pheno\trelative path to directory with one phenotypes file per subgroup\n"
        msg += "\t\trow 1 for sample names, column 1 for feature names\n"
        msg += "     --fcoord\tBED file with the features coordinates\n"
        msg += "\t\tfeatures should be in same order than in phenotypes files\n"
        msg += " -l, --links\tfile with links between genes and SNPs\n"
        msg += "\t\tcustom format: feature<space/tab>SNP|coord\n"
        msg += "\t\tfeatures should be in same order than in phenotypes files\n"
        msg += "\t\tuseful to focus on genetic variants in cis (use windowBed)\n"
        msg += " -f, --ftr\tfile with a list of features to analyze\n"
        msg += "\t\tone feature name per line\n"
        msg += "     --perm\tnumber of phenotype permutations at each feature\n"
        msg += "\t\tdefault=10000\n"
        msg += "     --seed\tseed for the random number generator\n"
        msg += "\t\tdefault=1859\n"
        msg += "     --grid\tfile with the grid of values for phi2 and omega2 (ES model)\n"
        msg += "\t\tsee GetGridPhiOmega() in package Rquantgen\n"
        msg += "\n"
        msg += "Example:\n"
        print msg; sys.stdout.flush()
        
        
    def version(self):
        msg = "%s 0.1\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "This is free software; see the source for copying conditions.  There is NO\n"
        msg += "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
        print msg; sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], "hVv:s:g:p:l:f:",
                                       ["help", "version", "verbose=",
                                        "step=", "geno=", "pheno=", "fcoord=",
                                        "links=", "ftr=", "perm=", "seed=",
                                        "grid="])
        except getopt.GetoptError, err:
            sys.stderr.write("%s\n" % str(err))
            self.help()
            sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                self.help()
                sys.exit(0)
            elif o in ("-V", "--version"):
                self.version()
                sys.exit(0)
            elif o in ("-v", "--verbose"):
                self.verbose = int(a)
            elif o in ("-s", "--step"):
                self.step = int(a)
            elif o in ("-g", "--geno"):
                if os.path.dirname(a) == "":
                    self.pathToGenoFile = "%s/%s" % (os.getcwd(), a)
                else:
                    self.pathToGenoFile = a
            elif o in ("-p", "--pheno"):
                if os.path.dirname(a) == "":
                    self.pathToPhenoDir = "%s/%s" % (os.getcwd(), a)
                else:
                    self.pathToPhenoDir = a
            elif o in ("--fcoord"):
                if os.path.dirname(a) == "":
                    self.pathToFtrCoords = "%s/%s" % (os.getcwd(), a)
                else:
                    self.pathToFtrCoords = a
            elif o in ("-l", "--links"):
                if os.path.dirname(a) == "":
                    self.pathToFtrsSnpsLinks = "%s/%s" % (os.getcwd(), a)
                else:
                    self.pathToFtrsSnpsLinks = a
            elif o in ("-f", "--ftr"):
                if os.path.dirname(a) == "":
                    self.pathToFtrsToKeep = "%s/%s" % (os.getcwd(), a)
                else:
                    self.pathToFtrsToKeep = a
            elif o in ("--perm"):
                self.nbPermutations = int(a)
            elif o in ("--seed"):
                seed = int(a)
            elif o in ("--grid"):
                if os.path.dirname(a) == "":
                    self.pathToGridFile = "%s/%s" % (os.getcwd(), a)
                else:
                    self.pathToGridFile = a
            else:
                assert False, "unhandled option"
                
                
    def checkAttributes(self):
        if self.step == 0:
            msg = "ERROR: step is missing (-s)"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.pathToPhenoDir == "":
            msg = "ERROR: phenotypes directory is missing (-p)"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.pathToPhenoDir):
            msg = "ERROR: directory '%s' doesn't exist" % self.pathToPhenoDir
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.pathToGenoFile == "":
            msg = "ERROR: genotypes file is missing (-g)"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.pathToGenoFile):
            msg = "ERROR: file '%s' doesn't exist" % self.pathToGenoFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.pathToFtrCoords == "":
            msg = "ERROR: features coordinates file is missing (--fcoord)"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.pathToFtrCoords):
            msg = "ERROR: file '%s' doesn't exist" % self.pathToFtrCoords
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.pathToFtrsSnpsLinks == "":
            msg = "ERROR: links feature-SNP file is missing (-l)"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.pathToFtrsSnpsLinks):
            msg = "ERROR: file '%s' doesn't exist" % self.pathToFtrsSnpsLinks
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if self.pathToGridFile == "":
            msg = "ERROR: grid file is missing (-g)"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.pathToGridFile):
            msg = "ERROR: file '%s' doesn't exist" % self.pathToGridFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        scipy.random.seed(self.seed)
        
        
    def parsePhenoDir(self):
        """Set the absolute path to the phenotypes files,
        as well as the number of individuals."""
        self.lPathToPhenoFiles = glob.glob("%s/%s/*" % (os.getcwd(),
                                                        self.pathToPhenoDir))
        for inFile in self.lPathToPhenoFiles:
            inH = open(inFile)
            line = inH.readline()
            tokens = line.split()
            if self.nbIndividuals == -1:
                self.nbIndividuals = len(tokens)
            elif len(tokens) != self.nbIndividuals:
                msg = "ERROR: directory '%s' doesn't exist" % self.pathToPhenoDir
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
            inH.close()
        if self.verbose > 1:
            print "nb of individuals: %i" % self.nbIndividuals
            sys.stdout.flush()
            
            
    def giveCommandsForFullAnalysisOnOriginalData(self):
        scriptFile = "commands_for_step_1.txt"
        scriptH = open(scriptFile, "w")
        
        cmd = "qsub -cwd -j y -V -l h_vmem=1G"
        
        txt = ""
        txt += "rm -rf %s\n" % self.originalAnalysisDir
        txt += "mkdir %s\n" % self.originalAnalysisDir
        txt += "cd %s\n" % self.originalAnalysisDir
        pathToSumstatsDir = "all_summary_stats"
        txt += "rm -rf %s\n" % pathToSumstatsDir
        txt += "mkdir %s\n" % pathToSumstatsDir
        scriptH.write("%s\n" % txt)
        
        # commands to compute summary stats in each subgroup
        txt = ""
        subgroup = 0
        for pathToPhenoFile in self.lPathToPhenoFiles:
            subgroup += 1
            phenoFile = os.path.basename(pathToPhenoFile)
            txt += "echo \""
            txt += "get_summary_stats"
            txt += " -g %s" % self.pathToGenoFile
            txt += " -p %s" % pathToPhenoFile
            txt += " --fcoord %s" % self.pathToFtrCoords
            txt += " -l %s" % self.pathToFtrsSnpsLinks
            if self.pathToFtrsToKeep != "":
                txt += " -f %s" % self.pathToFtrsToKeep
            txt += " -o %s/%s_subgroup%i_sumstats" % (pathToSumstatsDir, 
                                                      phenoFile, subgroup)
            txt += "\" | %s -N get_summary_stats_subgroup%i\n" % (cmd,
                                                                  subgroup)
        scriptH.write("%s\n" % txt)
        
        # commands to compute the ABF "meta" (joint analysis)
        txt = ""
        txt += "echo \""
        txt += "get_abf_meta"
        txt += " -i %s" % pathToSumstatsDir
        txt += " -g %s" % self.pathToGridFile
        txt += " -o %s" % self.originalAbfsFile
        txt += "\" | %s -N get_abf_meta\n" % cmd
        scriptH.write("%s\n" % txt)
        
        txt = ""
        txt += "cd ..\n"
        scriptH.write("%s\n" % txt)
        
        scriptH.close()
        
        print "you need to submit yourself the jobs from the file '%s'" % scriptFile
        print "in the mean time, you can launch steps 2, 3 and 4"
        
        
    def prepareDataForPermutations(self):
        if self.verbose > 0:
            print "prepare data for permutations ..."
            sys.stdout.flush()
            
        if os.path.exists(self.permDir):
            msg = "ERROR: directory '%s' already exists" % self.permDir
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
        os.mkdir(self.permDir)
        os.chdir(self.permDir)
        
        for permId in xrange(1, self.nbPermutations+1):
            newDir = "permutation_%i" % permId
            os.mkdir(newDir)
            os.chdir(newDir)
            permFile = "permutation_%i.txt" % permId
            aPerm = scipy.random.permutation(self.nbIndividuals)
            scipy.savetxt(fname=permFile, X=aPerm, fmt="%i")
            os.chdir("..")
            
        os.chdir("..")
        
        if self.verbose > 0:
            print "see directory '%s/'" % self.permDir
            sys.stdout.flush()
            
            
    def calcSumStatsForEachPermAndEachSubgroup(self):
        """Write a shell script to be launched as a job array with
        one job per permutation."""
        scriptFile = "script_for_array_job_step_3.sh"
        scriptH = open(scriptFile, "w")
        txt = "#!/bin/sh\n"
        txt += "#$ -t 1-%i\n" % self.nbPermutations
        txt += "uname -a\n"
        txt += "cd permutation_${SGE_TASK_ID}\n"
        for pathToPhenoFile in self.lPathToPhenoFiles:
            phenoFile = os.path.basename(pathToPhenoFile)
            txt += "get_summary_stats"
            txt += " -g %s" % self.pathToGenoFile
            txt += " -p %s" % pathToPhenoFile
            txt += " --fcoord %s" % self.pathToFtrCoords
            txt += " -l %s" % self.pathToFtrsSnpsLinks
            if self.pathToFtrsToKeep != "":
                txt += " -f %s" % self.pathToFtrsToKeep
            txt += " -o %s_perm${SGE_TASK_ID}_sumstats" % phenoFile
            txt += " --permf permutation_${SGE_TASK_ID}.txt"
            txt += "\n"
        txt += "cd .."
        scriptH.write("%s\n" % txt)
        
        cmd = "cd permutations; qsub -cwd -j y -V -N get_summary_stats_perm -l h_vmem=1G ../%s" % scriptFile
        scriptH.write("#%s\n" % cmd)
        scriptH.close()
        
        print "you need to submit the array job to SGE yourself:"
        print cmd
        print "once all jobs are finished, launch step 4"
        
        
    def calcAbfMetaForEachPermutation(self):
        """Write a shell script to be launched as a job array with
        one job per permutation."""
        scriptFile = "script_for_array_job_step_4.sh"
        scriptH = open(scriptFile, "w")
        txt = "#!/bin/sh\n"
        txt += "#$ -t 1-%i\n" % self.nbPermutations
        txt += "uname -a\n"
        txt += "cd permutation_${SGE_TASK_ID}\n"
        txt += "rm -rf all_summary_stats_perm${SGE_TASK_ID}\n"
        txt += "mkdir all_summary_stats_perm${SGE_TASK_ID}\n"
        txt += "cd all_summary_stats_perm${SGE_TASK_ID}\n"
        for pathToPhenoFile in self.lPathToPhenoFiles:
            phenoFile = os.path.basename(pathToPhenoFile)
            txt += "ln -s ../%s_perm${SGE_TASK_ID}_sumstats .\n" % phenoFile
        txt += "cd ..\n"
        txt += "get_abf_meta -i all_summary_stats_perm${SGE_TASK_ID}"
        txt += " -g %s" % self.pathToGridFile
        txt += " -o abfs_meta_perm${SGE_TASK_ID}.txt\n"
        txt += "cd .."
        scriptH.write("%s\n" % txt)
        
        cmd = "cd permutations; qsub -cwd -j y -V -N get_abf_meta_perm -l h_vmem=1G ../%s" % scriptFile
        scriptH.write("#%s\n" % cmd)
        scriptH.close()
        
        print "you need to submit the array job to SGE yourself:"
        print cmd
        print "once all jobs are finished, launch step 5"
        
        
    def getMaxAbfForEachNonPermutedFeature(self):
        dFtr2MaxAbf = {}
        
        abfsFile = "%s/%s" % (self.originalAnalysisDir,
                              self.originalAbfsFile)
        abfsH = open(abfsFile)
        line = abfsH.readline() # skip the header
        
        while True:
            line = abfsH.readline()
            if line =="":
                break
            tokens = line.rstrip().split()
            ftrName = tokens[0]
            abfMeta = float(tokens[5])
            if ftrName not in dFtr2MaxAbf:
                dFtr2MaxAbf[ftrName] = 0
            if abfMeta > dFtr2MaxAbf[ftrName]:
                dFtr2MaxAbf[ftrName] = abfMeta
                
        abfsH.close()
        
        return dFtr2MaxAbf
    
    
    def calcPermutationPvalueForOneFeature(self, lPathToAbfH, dFtr2MaxAbf):
        currentFtr = ""
        permPval = None
        
        # fill array with ABFs (SNPs in rows, permutations in columns)
        aAbfs = None
        nbSnps = 0
        endOfFile = False
        while True:
            lLines = []
            lPositions = []
            for pathToAbfH in lPathToAbfH:
                pos = pathToAbfH.tell()
                lPositions.append(pos)
                line = pathToAbfH.readline()
                if line == "":
                    break
                lLines.append(line)
            if lLines == []:
                endOfFile = True
                break
            if len(lLines) != self.nbPermutations:
                msg = "ERROR: ABF output files have different numbers"
                msg += " of lines for all permutations"
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
                
            tmp = scipy.zeros(self.nbPermutations)
            nextFtr = False
            for i in xrange(0, len(lLines)):
                tokens = lLines[i].rstrip().split()
                if currentFtr == "":
                    currentFtr = tokens[0]
                else:
                    if i == 0 and tokens[0] != currentFtr:
                        for i in xrange(0, len(lPathToAbfH)):
                            lPathToAbfH[i].seek(lPositions[i])
                        nextFtr = True
                        break
                    if i != 0 and tokens[0] != currentFtr:
                        msg = "ERROR: ABF output files have different numbers of SNPs"
                        msg += " for feature '%s'" % currentFtr
                        sys.stderr.write("%s\n" % msg)
                        sys.exit(1)
                tmp[i] = float(tokens[6])
            if nextFtr:
                break
            if aAbfs == None:
                aAbfs = tmp
            else:
                aAbfs = scipy.concatenate((aAbfs, tmp), axis=0)
            nbSnps += 1
            
        # calculate the permutation P-value
        if nbSnps != 0 or not endOfFile:
            aAbfs = scipy.reshape(aAbfs, (nbSnps, self.nbPermutations))
            aMaxPerPerm = scipy.amax(aAbfs, axis=0)
            permPval = (1 + sum(aMaxPerPerm >= dFtr2MaxAbf[currentFtr])) \
                / float(1 + self.nbPermutations)
            
        return currentFtr, permPval
    
    
    def calcPermutationPvaluesForAllFeatures(self):
        if self.verbose > 0:
            print "calculate feature-level permutation P-values ..."
            sys.stdout.flush()
            
        outH = open(self.outFile, "w")
        txt = "ftr permPval"
        outH.write("%s\n" % txt)
        
        dFtr2MaxAbf = self.getMaxAbfForEachNonPermutedFeature()
        if self.verbose > 0:
            print "nb of features: %i" % len(dFtr2MaxAbf)
            sys.stdout.flush()
            
        lPathToAbfFiles = glob.glob("%s/permutations/permutation_*/abfs_meta_perm*.txt" % os.getcwd())
        if len(lPathToAbfFiles) != self.nbPermutations:
            msg = "ERROR: %i ABF file(s)," % len(lPathToAbfFiles)
            msg +=" check if step 3 went well for all permutations"
            sys.stderr.write("%s\n" % msg)
            sys.exit(1)
        lPathToAbfH = []
        for pathToAbfFile in lPathToAbfFiles:
            pathToAbfH = open(pathToAbfFile)
            line = pathToAbfH.readline() # skip the header
            lPathToAbfH.append(pathToAbfH)
            
        while True:
            currentFtr, permPval = self.calcPermutationPvalueForOneFeature(lPathToAbfH,
                                                                           dFtr2MaxAbf)
            if currentFtr == "" and permPval == None:
                break
            txt = "%s %f" % (currentFtr, permPval)
            outH.write("%s\n" % txt)
            
        for pathToAbfH in lPathToAbfH:
            pathToAbfH.close()
            
        if self.verbose > 0:
            print "results written in file '%s'" % self.outFile
            sys.stdout.flush()
            
            
    def run(self):
        self.checkAttributes()
        
        if self.verbose > 0:
            msg = "START %s" % time.strftime("%Y-%m-%d %H:%M:%S")
            startTime = time.time()
            print msg; sys.stdout.flush()
            
        self.parsePhenoDir()
        
        if self.step == 1:
            self.giveCommandsForFullAnalysisOnOriginalData()
            
        if self.step == 2:
            self.prepareDataForPermutations()
            
        if self.step == 3:
            self.calcSumStatsForEachPermAndEachSubgroup()
            
        if self.step == 4:
            self.calcAbfMetaForEachPermutation()
            
        if self.step == 5:
            self.calcPermutationPvaluesForAllFeatures()
            
        if self.verbose > 0:
            msg = "END %s" % time.strftime("%Y-%m-%d %H:%M:%S")
            endTime = time.time()
            runLength = datetime.timedelta(seconds=
                                           math.floor(endTime - startTime))
            msg += " (%s)" % str(runLength)
            print msg; sys.stdout.flush()
            
            
if __name__ == "__main__":
    i = JointEqtlAnalysis()
    i.setAttributesFromCmdLine()
    i.run()
