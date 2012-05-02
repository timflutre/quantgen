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
        self.step = ""
        self.nbPermutations = 10000
        self.seed = 1859
        self.pathToPhenoDir = ""
        self.pathToGenoFile = ""
        self.pathToFtrCoords = ""
        self.pathToFtrsSnpsLinks = ""
        self.pathToGridFile = ""
        
        self.permFile = "permutations.txt.gz"
        self.permDir = "permutations"
        self.lPathToPhenoFiles = []
        self.nbIndividuals = -1
        
        
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
        msg += " -s, --step\tstep of the pipeline (1/2/3/4)\n"
        msg += "\t\t1: write all permuted phenotypes (one file per subgroup)\n"
        msg += "\t\t2: calculate the summary statistics for each permutation and subgroup\n"
        msg += "\t\t3: calculate the ABF_meta for each permutation\n"
        msg += "\t\t4: calculate the gene-level permutation P-values for the joint analysis\n"
        msg += " -g, --geno\tfile with genotypes in IMPUTE format"
        msg += "\t\ta header line with sample names is required\n"
        msg += "\t\tsamples in columns should be in same order as in phenotype file\n"
        msg += " -p, --pheno\trelative path to directory with one phenotypes file per subgroup\n"
        msg += "\t\trow 1 for sample names, column 1 for feature names, delimiter=<space/tab>\n"
        msg += "     --fcoord\tBED file with the features coordinates\n"
        msg += "\t\tfeatures should be in same order than in phenotypes files\n"
        msg += " -l, --links\tfile with links between genes and SNPs\n"
        msg += "\t\tcustom format: feature<space/tab>SNP|coord\n"
        msg += "\t\tfeatures should be in same order than in phenotypes files\n"
        msg += "\t\tuseful to focus on genetic variants in cis (use windowBed)\n"
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
            opts, args = getopt.getopt(sys.argv[1:], "hVv:s:g:p:l:",
                                       ["help", "version", "verbose=",
                                        "step=", "geno=", "pheno=", "fcoord=",
                                        "links=", "perm=", "seed=", "grid="])
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
                self.step = a
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
        if self.step == "":
            msg = "ERROR: step is missing (-s)"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if "1" in self.step or "2" in self.step:
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
        if "2" in self.step:
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
        if "3" in self.step:
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
        self.lPathToPhenoFiles = glob.glob("%s/%s/*" % (os.getcwd(), self.pathToPhenoDir))
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
        if self.verbose > 0:
            print "nb of individuals: %i" % self.nbIndividuals
            sys.stdout.flush()
            
            
    def writePhenotypesForOnePermutation(self, permId, aPerm, lPathToPhenoH):
        """Write the permuted phenotypes for each subgroup for a given permutation."""
        newDir = "permutation_%s" % str(permId).zfill(len(str(self.nbPermutations)))
        os.mkdir(newDir)
        os.chdir(newDir)
        
        for i in xrange(0, len(lPathToPhenoH)):
            pathToPhenoFile = self.lPathToPhenoFiles[i]
            pathToPhenoH = lPathToPhenoH[i]
            pathToPhenoH.seek(0)
            phenoFile = "%s_perm%i" % (os.path.basename(pathToPhenoFile), permId)
            phenoH = open(phenoFile, "w")
            
            header = pathToPhenoH.readline()
            phenoH.write(header)
            
            while True:
                line = pathToPhenoH.readline()
                if line == "":
                    break
                tokens = line.rstrip().split()
                txt = tokens[0]
                for j in xrange(0, self.nbIndividuals):
                    txt += " %s" % tokens[aPerm[j]]
                phenoH.write("%s\n" % txt)
                
            phenoH.close()
            
        os.chdir("..")
        
        
    def writePhenotypesForAllPermutations(self):
        """Write the permuted phenotypes for each subgroup and each permutation."""
        if self.verbose > 0:
            print "write permuted phenotypes ..."
            sys.stdout.flush()
            
        self.parsePhenoDir()
        lPathToPhenoH = []
        for pathToPhenoFile in self.lPathToPhenoFiles:
            pathToPhenoH = open(pathToPhenoFile)
            lPathToPhenoH.append(pathToPhenoH)
            
        if os.path.exists(self.permDir):
            shutil.rmtree(self.permDir)
        os.mkdir(self.permDir)
        os.chdir(self.permDir)
        
        permH = gzip.open(self.permFile, "w")
        
        for permId in xrange(1, self.nbPermutations+1):
            if self.verbose > 1:
                print "%s / %i" % (str(permId).zfill(len(str(self.nbPermutations))),
                                   self.nbPermutations)
                
            aPerm = scipy.random.permutation(self.nbIndividuals)
            scipy.savetxt(fname=permH, X=aPerm, fmt="%i", newline=" ")
            permH.write("\n")
            
            self.writePhenotypesForOnePermutation(permId, aPerm, lPathToPhenoH)
            
        permH.close()
        
        for pathToPhenoH in lPathToPhenoH:
            pathToPhenoH.close()
            
        os.chdir("..")
        
        if self.verbose > 0:
            print "permuted phenotypes written in '%s/'" % self.permDir
            sys.stdout.flush()
            
            
    def calcSumStatsForEachPermAndEachSubgroup(self):
        """Write a shell script to be launched as a job array with
        one job per permutation."""
        self.parsePhenoDir()
        scriptFile = "script_for_array_job_step_1.sh"
        scriptH = open(scriptFile, "w")
        txt = "#!/bin/sh\n"
        txt += "#$ -t 1-%i\n" % self.nbPermutations
        txt += "uname -a\n"
        txt += "cd permutations/permutation_$SGE_TASK_ID\n"
        for pathToPhenoFile in self.lPathToPhenoFiles:
            phenoFile = os.path.basename(pathToPhenoFile)
            txt += "get_summary_stats"
            txt += " -g %s" % self.pathToGenoFile
            txt += " -p %s_perm$SGE_TASK_ID" % phenoFile
            txt += " --fcoord %s" % self.pathToFtrCoords
            txt += " -l %s" % self.pathToFtrsSnpsLinks
            txt += " -o %s_perm$SGE_TASK_ID_sumstats" % phenoFile
            txt += "\n"
        txt += "cd .."
        scriptH.write("%s\n" % txt)
        
        cmd = "qsub -cwd -j y -V -N get_summary_stats -l h_vmem=1G %s" % scriptFile
        scriptH.write("#%s\n" % cmd)
        scriptH.close()
        
        print "you need to submit the array job to SGE yourself, like this:"
        print cmd
        print "once all jobs are finished (monitor them with 'qstat'), launch step 3"
        
        
    def calcAbfMetaForEachPermutation(self):
        """Write a shell script to be launched as a job array with
        one job per permutation."""
        self.parsePhenoDir()
        scriptFile = "script_for_array_job_step_2.sh"
        scriptH = open(scriptFile, "w")
        txt = "#!/bin/sh\n"
        txt += "#$ -t 1-%i\n" % self.nbPermutations
        txt += "uname -a\n"
        txt += "cd permutations/permutation_$SGE_TASK_ID\n"
        txt += "rm -rf all_summary_stats_$SGE_TASK_ID\n"
        txt += "mkdir all_summary_stats_$SGE_TASK_ID\n"
        txt += "cd all_summary_stats_$SGE_TASK_ID\n"
        for pathToPhenoFile in self.lPathToPhenoFiles:
            phenoFile = os.path.basename(pathToPhenoFile)
            txt += "ln -s ../%s_perm$SGE_TASK_ID_sumstats .\n" % phenoFile
        txt += "cd .."
        txt += "get_abf_meta -i all_summary_stats_$SGE_TASK_ID\n"
        txt += " -g %s\n" % self.pathToGridFile
        txt += " -o abfs_meta_perm$SGE_TASK_ID.txt.gz\n"
        txt += "cd .."
        scriptH.write("%s\n" % txt)
        
        cmd = "qsub -cwd -j y -V -N get_abf_meta -l h_vmem=1G %s" % scriptFile
        scriptH.write("#%s\n" % cmd)
        scriptH.close()
        
        print "you need to submit the array job to SGE yourself, like this:"
        print cmd
        print "once all jobs are finished (monitor them with 'qstat'), launch step 3"
        
        
    def run(self):
        self.checkAttributes()
        
        if self.verbose > 0:
            msg = "START %s" % time.strftime("%Y-%m-%d %H:%M:%S")
            startTime = time.time()
            print msg; sys.stdout.flush()
            
        if "1" in self.step:
            self.writePhenotypesForAllPermutations()
            
        if "2" in self.step:
            self.calcSumStatsForEachPermAndEachSubgroup()
            
        if "3" in self.step:
            self.calcAbfMetaForEachPermutation()
            
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
