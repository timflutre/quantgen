/** \file get_ftr-level_perm-pval.cpp
 *
 *  `get_ftr-level_perm-pval' computes feature-level permutation P-values.
 *  Copyright (C) 2011,2012 Timothee Flutre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  g++ -Wall -fopenmp -O3 get_ftr-level_perm-pval.cpp -lgsl -lgslcblas -o get_ftr-level_perm-pval
 *  help2man -o .man ./get_ftr-level_perm-pval
 *  groff -mandoc get_ftr-level_perm-pval.man > get_ftr-level_perm-pval.ps
*/

#include <cmath>
#include <ctime>
#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
using namespace std;

#include "utils.cpp"

/** \brief Display the help on stdout.
*/
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " computes feature-level permutation " << endl
       << "P-values for the joint analysis." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "  -g, --geno\tfile with genotypes in IMPUTE format" << endl
       << "\t\ta header line with sample names is required" << endl
       << "\t\tsamples in columns should be in same order as in phenotype file" << endl
       << "  -p, --pheno\trelative path to directory with one phenotypes file per subgroup" << endl
       << "\t\trow 1 for sample names, column 1 for feature names" << endl
       << "      --fcoord\tBED file with the features coordinates" << endl
       << "\t\tfeatures should be in same order than in phenotypes files" << endl
       << "  -l, --links\tfile with links between genes and SNPs" << endl
       << "\t\tcustom format: feature<space/tab>SNP|coord" << endl
       << "\t\tfeatures should be in same order than in phenotypes files" << endl
       << "\t\tuseful to focus on genetic variants in cis (use windowBed)" << endl
       << "      --grid\tfile with the grid of values for phi2 and omega2 (ES model)" << endl
       << "\t\tsee GetGridPhiOmega() in package Rquantgen" << endl
       << "  -f, --ftr\tfile with a list of features to analyze" << endl
       << "\t\tone feature name per line" << endl
       << "\t\tcompulsory to force parallelization on a cluster" << endl
       << "      --truth\tfile with the ABFs on the original dataset" << endl
       << "\t\toutput from get_abf_meta" << endl
       << "      --abf\twhich ABF to use as the test statistic for the permutations" << endl
       << "\t\tdefault=abf.meta/abf.fix/abf.meta.avg.all/abf.meta.avg.subset" << endl
       << "  -o, --out\tfile that will contain the permutation P-values" << endl
       << "      --perm\tnumber of phenotype permutations at each feature" << endl
       << "\t\tdefault=10000" << endl
       << "      --trick\tstop after the tenth permutation for which the test statistic" << endl
       << "\t\tis larger than or equal to the true value, and sample from a" << endl
       << "\t\tuniform between 11/(nbPerm+2) and 11/(nbPerm+1)" << endl
       << "      --seed\tseed for the random number generator" << endl
       << "\t\tdefault=1859" << endl
       << "  -t, --thread\tnumber of threads (default=1)" << endl
       << "\t\tused for SNPs in cis of the same feature (get_summary_stats)" << endl
       << endl
    ;
}

/** \brief Display version and license information on stdout.
 */
void version (char ** argv)
{
  cout << argv[0] << " 0.1" << endl
       << endl
       << "Copyright (C) 2011,2012 T. Flutre." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by T. Flutre." << endl;
}

/** \brief Parse the command-line arguments and check the values of the 
 *  compulsory ones.
 */
void
parse_args (
  int argc,
  char ** argv,
  string & genoFile,
  string & phenoDir,
  string & ftrCoordsFile,
  string & linksFile,
  string & gridFile,
  string & ftrFile,
  string & truthFile,
  string & whichAbf,
  string & outFile,
  size_t & nbPermutations,
  bool & trickPerm,
  int & seed,
  int & nbThreads,
  int & verbose)
{
  int c = 0;
  while (1)
  {
    static struct option long_options[] =
      {
	{"help", no_argument, 0, 'h'},
	{"version", no_argument, 0, 'V'},
	{"verbose", required_argument, 0, 'v'},
	{"geno", required_argument, 0, 'g'},
	{"pheno", required_argument, 0, 'p'},
	{"fcoord", required_argument, 0, 0},
	{"links", required_argument, 0, 'l'},
	{"grid", required_argument, 0, 0},
	{"ftr", required_argument, 0, 'f'},
	{"truth", required_argument, 0, 0},
	{"abf", required_argument, 0, 0},
	{"out", required_argument, 0, 'o'},
	{"perm", required_argument, 0, 0},
	{"trick", required_argument, 0, 0},
	{"seed", required_argument, 0, 0},
	{"thread", required_argument, 0, 't'},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:l:f:o:t:",
		     long_options, &option_index);
    if (c == -1)
      break;
    switch (c)
    {
    case 0:
      if (long_options[option_index].flag != 0)
	break;
      if (strcmp(long_options[option_index].name, "fcoord") == 0)
      {
	ftrCoordsFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "grid") == 0)
      {
	gridFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "truth") == 0)
      {
	truthFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "abf") == 0)
      {
	whichAbf = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "perm") == 0)
      {
	nbPermutations = atol(optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "trick") == 0)
      {
	trickPerm = true;
	break;
      }
      if (strcmp(long_options[option_index].name, "seed") == 0)
      {
	seed = atoi(optarg);
	break;
      }
    case 'h':
      help (argv);
      exit (0);
    case 'V':
      version (argv);
      exit (0);
    case 'v':
      verbose = atoi(optarg);
      break;
    case 'g':
      genoFile = optarg;
      break;
    case 'p':
      phenoDir = optarg;
      break;
    case 'l':
      linksFile = optarg;
      break;
    case 'f':
      ftrFile = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 't':
      nbThreads = atoi(optarg);
      break;
    case '?':
      printf ("\n"); help (argv);
      abort ();
    default:
      printf ("\n"); help (argv);
      abort ();
    }
  }
  if (genoFile.empty())
  {
    fprintf (stderr, "ERROR: missing file with genotypes (-g).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (genoFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", genoFile.c_str());
    help (argv);
    exit (1);
  }
  if (phenoDir.empty())
  {
    fprintf (stderr, "ERROR: missing directory with phenotypes (-p).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (phenoDir))
  {
    fprintf (stderr, "ERROR: can't find directory '%s'.\n\n", phenoDir.c_str());
    help (argv);
    exit (1);
  }
  if (ftrCoordsFile.empty())
  {
    fprintf (stderr, "ERROR: missing feature coordinates file (--fcoord).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (ftrCoordsFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", ftrCoordsFile.c_str());
    help (argv);
    exit (1);
  }
  if (linksFile.empty())
  {
    fprintf (stderr, "ERROR: missing links file (-l).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (linksFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", linksFile.c_str());
    help (argv);
    exit (1);
  }
  if (ftrFile.empty())
  {
    fprintf (stderr, "ERROR: missing feature file (-f).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (ftrFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", ftrFile.c_str());
    help (argv);
    exit (1);
  }
  if (truthFile.empty())
  {
    fprintf (stderr, "ERROR: missing truth file (--truth).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (truthFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", truthFile.c_str());
    help (argv);
    exit (1);
  }
  if (outFile.empty())
  {
    fprintf (stderr, "ERROR: missing output file (-o).\n\n");
    help (argv);
    exit (1);
  }
}

map<string, double>
getTrueAbfForEachFtr (
  const string & truthFile,
  const string & whichAbf,
  const vector<string> & vFtrsToKeep,
  const int & verbose)
{
  size_t lineId = 1;
  string line;
  vector<string> tokens;
  map<string, double> mFtr2Abf;
  ifstream truthStream;
  
  truthStream.open(truthFile.c_str());
  if (! truthStream.is_open())
  {
    cerr << "ERROR: can't open file " << truthFile << endl;
    exit (1);
  }
  if (verbose > 0)
    cout <<"retrieve true ABFs from file " << truthFile << "..." << endl;
  
  // check first line is proper header
  getline (truthStream, line);
  if (line.empty())
  {
    cerr << "ERROR: file " << truthFile << " is empty" << endl;
    exit (1);
  }
  split (line, ' ', tokens);
  if ((whichAbf.compare("abf.meta") == 0
       || whichAbf.compare("abf.fix") == 0)
      && tokens.size() < 7)
  {
    cerr << "ERROR: file " << truthFile << " should have at least"
	 << " 7 columns in his header" << endl;
    exit (1);
  }
  if ((whichAbf.compare("abf.meta.avg.all") == 0
       || whichAbf.compare("abf.meta.avg.subset") == 0)
      && tokens.size() < 10)
  {
    cerr << "ERROR: file " << truthFile << " should have at least"
	 << " 10 columns in his header" << endl;
    exit (1);
  }
  if (tokens[0].compare("ftr") != 0
      || tokens[1].compare("snp") != 0
      || tokens[2].compare("coord") != 0
      || tokens[3].compare("nb.subgroups") != 0
      || tokens[4].compare("nb.samples") != 0
      || tokens[5].compare("l10abf.meta") != 0
      || tokens[6].compare("l10abf.fix") != 0
      || tokens[7].compare("l10abf.maxh") != 0)
  {
    cerr << "ERROR: header line of file " << truthFile
	 << " should start like this:" << endl
	 << "ftr snp coord nb.subgroups nb.samples l10abf.meta l10abf.fix l10abf.maxh" << endl;
    exit (1);
  }
  if ((whichAbf.compare("abf.meta.avg.all") == 0
       || whichAbf.compare("abf.meta.avg.subset") == 0)
      && (tokens[8].compare("l10abf.meta.avg.all.c") != 0
	  || tokens[9].compare("l10abf.meta.avg.subset.c") != 0))
  {
    cerr << "ERROR: header line of file " << truthFile
	 << " should start like this:" << endl
	 << "ftr snp coord nb.subgroups nb.samples l10abf.meta l10abf.fix l10abf.maxh l10abf.meta.avg.all.c l10abf.meta.avg.subset.c" << endl;
    exit (1);
  }
  ++lineId;
  
  while (truthStream.good())
  {
    getline (truthStream, line);
    if (line.empty())
      break;
    split (line, ' ', tokens);
    if (find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0])
	!= vFtrsToKeep.end())
    {
      if (whichAbf.compare("abf.meta") == 0)
	mFtr2Abf.insert ( make_pair (tokens[0], atof(tokens[5].c_str())));
      if (whichAbf.compare("abf.fix") == 0)
	mFtr2Abf.insert ( make_pair (tokens[0], atof(tokens[6].c_str())));
      if (whichAbf.compare("abf.meta.avg.all") == 0)
	mFtr2Abf.insert ( make_pair (tokens[0], atof(tokens[8].c_str())));
      if (whichAbf.compare("abf.meta.avg.subset") == 0)
	mFtr2Abf.insert ( make_pair (tokens[0], atof(tokens[9].c_str())));
    }
    ++lineId;
  }
  
  truthStream.close();
  
  if (mFtr2Abf.size() == 0)
  {
    cerr << "ERROR: can't retrieve any true ABFs" << endl;
    exit (1);
  }
  
  if (verbose > 0)
    cout << "ABFs retrieved: " << mFtr2Abf.size() << endl;
  
  return mFtr2Abf;
}

void 
computeJointAnalysisPermPvaluesForOneFtr (
  const string & genoFile,
  const string & phenoDir,
  const string & ftrCoordsFile,
  const string & linksFile,
  const string & gridFile,
  const string & ftrToKeep,
  const double & trueAbf,
  const string & whichAbf,
  double & jointPermPval,
  const size_t & nbPermutations,
  const bool & trickPerm,
  const int & seed,
  const int & nbThreads,
  const int & verbose)
{
  size_t permId;
  
  jointPermPval = 0.0;
  
  for (permId=0; permId<nbPermutations; ++permId)
  {
    // get permutations and write it into file
    
    // calc sumstats in each subgroup for this permutation
    
    // calc ABFs from these sumstats
    
    // compare test stat to true value
    
  }
}

void computeJointAnalysisPermPvaluesFtrPerFtr (
  const string & genoFile,
  const string & phenoDir,
  const string & ftrCoordsFile,
  const string & linksFile,
  const string & gridFile,
  const vector<string> & vFtrsToKeep,
  const map<string, double> & mFtr2Abf,
  const string & whichAbf,
  const string & outFile,
  const size_t & nbPermutations,
  const bool & trickPerm,
  const int & seed,
  const int & nbThreads,
  const int & verbose)
{
  vector<string>::const_iterator it;
  double jointPermPval;
  ofstream outStream;
  
  outStream.open(outFile.c_str());
  if (! outStream.is_open())
  {
    cerr << "ERROR: can't open file " << outFile << endl;
    exit (1);
  }
  outStream << "ftr " << whichAbf << " jointPermPval" << endl;
  
  for (it = vFtrsToKeep.begin(); it != vFtrsToKeep.end(); ++it)
  {
    computeJointAnalysisPermPvaluesForOneFtr (genoFile,
					      phenoDir,
					      ftrCoordsFile,
					      linksFile,
					      gridFile,
					      *it,
					      mFtr2Abf.find(*it)->second,
					      whichAbf,
					      jointPermPval,
					      nbPermutations,
					      trickPerm,
					      seed,
					      nbThreads,
					      verbose-1);
    outStream << *it << " " << mFtr2Abf.find(*it)->second
	      << " " << jointPermPval << endl;
  }
  
  outStream.close();
}

int main (int argc, char ** argv)
{
  string genoFile, phenoDir, ftrCoordsFile, linksFile, gridFile, ftrFile,
    truthFile, whichAbf = "abf.meta", outFile;
  size_t nbPermutations = 10000;
  bool trickPerm = false;
  int verbose = 1, seed = 1859, nbThreads = 1;
  
  parse_args (argc, argv, genoFile, phenoDir, ftrCoordsFile, linksFile,
	      gridFile, ftrFile, truthFile, whichAbf, outFile, nbPermutations,
	      trickPerm, seed, nbThreads, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl;
  }
  
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrFile, verbose);
  map<string, double> mFtr2Abf = getTrueAbfForEachFtr (truthFile,
						       whichAbf,
						       vFtrsToKeep,
						       verbose);
  
  computeJointAnalysisPermPvaluesFtrPerFtr (genoFile,
					    phenoDir,
					    ftrCoordsFile,
					    linksFile,
					    gridFile,
					    vFtrsToKeep,
					    mFtr2Abf,
					    whichAbf,
					    outFile,
					    nbPermutations,
					    trickPerm,
					    seed,
					    nbThreads,
					    verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime)
	 << ": elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << ")" << endl;
  }
  
  return EXIT_SUCCESS;
}
