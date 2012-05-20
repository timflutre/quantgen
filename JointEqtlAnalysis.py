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


class JointEqtlAnalysis(object):
    
    def __init__(self):
        self.verbose = 1
        self.step = 0
        self.nbPermutationsSep = 0
        self.nbPermutationsJoint = 10000
        self.seed = 1859
        self.pathToPhenoDir = ""
        self.pathToGenoFile = ""
        self.pathToFtrCoords = ""
        self.pathToFtrsSnpsLinks = ""
        self.pathToFtrsToKeep = ""
        self.pathToGridFile = ""
        self.permDir = "permutations"
        self.nbThreads = 1
        self.whichAbf = "abf.meta"
        
        self.originalAnalysisDir = "original_analysis"
        self.originalAbfsFile = "abfs_meta.txt"
        self.lPathToPhenoFiles = []
        self.nbIndividuals = -1
        self.outFile = ""
        
        
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
        msg += "\t\t1: give commands to launch full analysis on the original data\n"
        msg += "\t\t2: write the permutations\n"
        msg += "\t\t3: give commands to launch full analysis on the permuted data\n"
        msg += "\t\t4: calculate the gene-level joint-analysis permutation P-values\n"
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
        msg += "     --perm1\tnumber of phenotype permutations at each feature\n"
        msg += "\t\tfor the separate analysis analysis of step 1, default=0\n"
        msg += "     --perm2\tnumber of phenotype permutations at each feature\n"
        msg += "\t\tfor the joint analysis of steps 2-3-4, default=10000\n"
        msg += "     --permdir\tname of the directory for the permutations\n"
        msg += "\t\tdefault=permutations\n"
        msg += "     --seed\tseed for the random number generator\n"
        msg += "\t\tdefault=1859\n"
        msg += "     --grid\tfile with the grid of values for phi2 and omega2 (ES model)\n"
        msg += "\t\tsee GetGridPhiOmega() in package Rquantgen\n"
        msg += " -t, --thread\tnumber of threads (default=1)\n"
        msg += "\t\tused for SNPs in cis of the same feature (get_summary_stats)\n"
        msg += "     --abf\twhich ABF to use as the test statistic for the permutations\n"
        msg += "\t\tdefault=abf.meta/abf.fix/abf.meta.avg.all/abf.meta.avg.subset\n"
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
            opts, args = getopt.getopt(sys.argv[1:], "hVv:s:g:p:l:f:t:",
                                       ["help", "version", "verbose=",
                                        "step=", "geno=", "pheno=", "fcoord=",
                                        "links=", "ftr=", "perm1=", "perm2=", 
                                        "seed=", "grid=", "permdir=",
                                        "thread=", "abf="])
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
            elif o in ("--perm1"):
                self.nbPermutationsSep = int(a)
            elif o in ("--perm2"):
                self.nbPermutationsJoint = int(a)
            elif o in ("--seed"):
                seed = int(a)
            elif o in ("--grid"):
                if os.path.dirname(a) == "":
                    self.pathToGridFile = "%s/%s" % (os.getcwd(), a)
                else:
                    self.pathToGridFile = a
            elif o in ("--permdir"):
                self.permDir = a
            elif o in ("-t", "--thread"):
                self.nbThreads = int(a) if int(a) > 1 else 1
            elif o in ("--abf"):
                self.whichAbf = a
            else:
                assert False, "unhandled option"
        self.outFile = "ftr-level-pvalues_%s_%s.txt" % (self.permDir,
                                                        self.whichAbf.replace(".","-"))
        
        
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
        if self.whichAbf not in ("abf.meta", "abf.fix", "abf.meta.avg.all",
                                 "abf.meta.avg.subset"):
            msg = "ERROR: can't recognize option --abf %s" % self.whichAbf
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
        
        cmd = "qsub -cwd -j y -V -l h_vmem=2G"
        
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
            if self.nbThreads > 1:
                txt += " -t %i" % self.nbThreads
            if self.nbPermutationsSep > 0:
                txt += " --perm %i" % self.nbPermutationsSep
            txt += "\" | %s" % cmd
            if self.nbThreads > 1:
                txt += " -pe simple_pe %i" % self.nbThreads
            txt += " -N get_summary_stats_subgroup%i\n" % subgroup
            txt += "\n"
        scriptH.write("%s\n" % txt)
        
        # commands to compute the ABF "meta" (joint analysis)
        txt = ""
        txt += "echo \""
        txt += "get_abf_meta"
        txt += " -i %s" % pathToSumstatsDir
        txt += " -g %s" % self.pathToGridFile
        txt += " -o %s" % self.originalAbfsFile
        if self.whichAbf in ("abf.meta.avg.all", "abf.meta.avg.subset"):
            txt += " -c"
        txt += "\" | %s -N get_abf_meta\n" % cmd
        scriptH.write("%s\n" % txt)
        
        txt = ""
        txt += "cd ..\n"
        scriptH.write("%s\n" % txt)
        
        scriptH.close()
        
        print "you need to submit yourself the jobs from the file '%s'" % scriptFile
        print "in the mean time, you can launch steps 2 and then 3"
        
        
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
        
        for permId in xrange(1, self.nbPermutationsJoint+1):
            permFile = "permutation_%i.txt" % permId
            aPerm = scipy.random.permutation(self.nbIndividuals)
            scipy.savetxt(fname=permFile, X=aPerm, fmt="%i")
            
        os.chdir("..")
        
        if self.verbose > 0:
            print "see directory '%s/'" % self.permDir
            sys.stdout.flush()
            
            
    def giveCommandsForFullAnalysisOnPermutedData(self):
        scriptFile = "jobarray_%s.bash" % self.permDir
        scriptH = open(scriptFile, "w")
        
        cmd = "qsub -cwd -j y -V -l h_vmem=2G -t 1-%i" % self.nbPermutationsJoint
        cmd += " -N jobarray_%s" % self.permDir
        if self.nbThreads > 1:
            cmd += " -pe simple_pe %i" % self.nbThreads
            
        txt = "#!/usr/bin/env bash\n"
        txt += "set -o errexit -o pipefail\n"
        txt += "mytmpdir=${TMPDIR}/tmp_${USER}_$$\n"
        txt += "trap \"cd; rm -rf $mytmpdir; exit\" INT TERM EXIT\n"
        txt += "uname -a\n"
        txt += "mycwd=$(pwd); echo \"mycwd: \"$mycwd\n"
        txt += "echo \"mytmpdir: \"$mytmpdir\n"
        txt += "rm -rf $mytmpdir; mkdir $mytmpdir; cd $mytmpdir\n"
        scriptH.write("%s\n" % txt)
        
        # commands to compute summary stats in each subgroup
        txt = ""
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
            txt += " --permf ${mycwd}/permutation_${SGE_TASK_ID}.txt"
            if self.nbThreads > 1:
                txt += " -t %i" % self.nbThreads
            txt += "\n"
        scriptH.write("%s\n" % txt)
        
        # commands to compute the ABF "meta" (joint analysis)
        txt = ""
        txt += "rm -rf all_sumstats_perm${SGE_TASK_ID};"
        txt += " mkdir all_sumstats_perm${SGE_TASK_ID};"
        txt += " cd all_sumstats_perm${SGE_TASK_ID}\n"
        for pathToPhenoFile in self.lPathToPhenoFiles:
            phenoFile = os.path.basename(pathToPhenoFile)
            txt += "ln -s ../%s_perm${SGE_TASK_ID}_sumstats .\n" % phenoFile
        txt += "cd ..\n"
        txt += "get_abf_meta"
        txt += " -i all_sumstats_perm${SGE_TASK_ID}"
        txt += " -g %s" % self.pathToGridFile
        txt += " -o abfs_perm${SGE_TASK_ID}.txt"
        if self.whichAbf in ("abf.meta.avg.all", "abf.meta.avg.subset"):
            txt += " -c"
        txt += "\n"
        txt += "rm -rf *sumstats*\n"
        scriptH.write("%s\n" % txt)
        
        # commands to extract results and remove temporary files
        txt = ""
        txt += "awk 'NR>1 {if(! ($1 in a)){a[$1] = $6}"
        txt += " else{if($6 > a[$1]){a[$1]=$6}}}"
        txt += " END{for(i in a) print i,a[i]}' abfs_perm${SGE_TASK_ID}.txt"
        txt += " | sort -k1,1"
        txt += " | gzip > best_abf-meta_per_ftr_perm${SGE_TASK_ID}.txt.gz\n"
        txt += "awk 'NR>1 {if(! ($1 in a)){a[$1] = $7}"
        txt += " else{if($7 > a[$1]){a[$1]=$7}}}"
        txt += " END{for(i in a) print i,a[i]}' abfs_perm${SGE_TASK_ID}.txt"
        txt += " | sort -k1,1"
        txt += " | gzip > best_abf-fix_per_ftr_perm${SGE_TASK_ID}.txt.gz\n"
        if self.whichAbf in ("abf.meta.avg.all", "abf.meta.avg.subset"):
            txt += "awk 'NR>1 {if(! ($1 in a)){a[$1] = $9}"
            txt += " else{if($9 > a[$1]){a[$1]=$9}}}"
            txt += " END{for(i in a) print i,a[i]}' abfs_perm${SGE_TASK_ID}.txt"
            txt += " | sort -k1,1"
            txt += " | gzip > best_abf-meta-avg-all_per_ftr_perm${SGE_TASK_ID}.txt.gz\n"
            txt += "awk 'NR>1 {if(! ($1 in a)){a[$1] = $10}"
            txt += " else{if($10 > a[$1]){a[$1]=$10}}}"
            txt += " END{for(i in a) print i,a[i]}' abfs_perm${SGE_TASK_ID}.txt"
            txt += " | sort -k1,1"
            txt += " | gzip > best_abf-meta-avg-subset_per_ftr_perm${SGE_TASK_ID}.txt.gz\n"
        txt += "cp best_abf-*_per_ftr_perm${SGE_TASK_ID}.txt.gz ${mycwd}/\n"
        txt += "rm -f abfs_perm${SGE_TASK_ID}.txt\n"
        txt += "cd $mycwd\n"
        scriptH.write("%s\n" % txt)
        
        cmd2 = "cd %s;" % self.permDir
        cmd2 += " echo \"time bash ../%s\" | %s" % (scriptFile, cmd)
        cmd2 += "; cd .."
        scriptH.write("#%s\n" % cmd2)
        
        scriptH.close()
        
        print "you need to submit yourself the job array:"
        print cmd2
        
        
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
            abfFix = float(tokens[6])
            if ftrName not in dFtr2MaxAbf:
                dFtr2MaxAbf[ftrName] = [0, 0]
            if abfMeta > dFtr2MaxAbf[ftrName][0]:
                dFtr2MaxAbf[ftrName][0] = abfMeta
            if abfFix > dFtr2MaxAbf[ftrName][1]:
                dFtr2MaxAbf[ftrName][1] = abfFix
                
        abfsH.close()
        
        return dFtr2MaxAbf
    
    
    def calcPermutationPvaluesForAllFeatures(self):
        if self.verbose > 0:
            print "calculate feature-level permutation P-values ..."
            sys.stdout.flush()
            
        outH = open(self.outFile, "w")
        txt = "ftr permPvalAbfMeta permPvalAbfFix"
        outH.write("%s\n" % txt)
        
        dFtr2MaxAbf = self.getMaxAbfForEachNonPermutedFeature()
        if self.verbose > 0:
            print "nb of features: %i" % len(dFtr2MaxAbf)
            sys.stdout.flush()
            
        dFtr2PermPval = {}
        for i in xrange(0, self.nbPermutationsJoint):
            pathToAbfFile = "%s/best_abf-meta_per_ftr_perm%i.txt.gz" % (
                self.permDir, i+1)
            pathToAbfH = gzip.open(pathToAbfFile)
            while True:
                line = pathToAbfH.readline()
                if line == "":
                    break
                tokens = line.rstrip().split()
                if len(tokens) != 2:
                    msg = "ERROR: check step 3 for permutation %i" % (i+1)
                    sys.stderr.write("%s\n" % msg)
                    sys.exit(1)
                if not dFtr2PermPval.has_key(tokens[0]):
                    dFtr2PermPval[tokens[0]] = [1, 1]
                if float(tokens[1]) >= dFtr2MaxAbf[tokens[0]][0]:
                    dFtr2PermPval[tokens[0]][0] += 1
                if float(tokens[1]) >= dFtr2MaxAbf[tokens[0]][1]:
                    dFtr2PermPval[tokens[0]][1] += 1
            pathToAbfH.close()
            
        lFeatures = dFtr2PermPval.keys()
        lFeatures.sort()
        for ftr in lFeatures:
            txt = "%s %f %f" % (ftr,
                                dFtr2PermPval[ftr][0] / float(1 + self.nbPermutationsJoint),
                                dFtr2PermPval[ftr][1] / float(1 + self.nbPermutationsJoint))
            outH.write("%s\n" % txt)
            
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
            self.giveCommandsForFullAnalysisOnPermutedData()
            
        if self.step == 4:
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
