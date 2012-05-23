/** \file get_ftr-level_perm-pval.cpp
 *
 *  `get_ftr-level_perm-pval' computes feature-level permutation P-values.
 *  Copyright (C) 2012 Timothee Flutre
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
 *  g++ -Wall -O3 get_ftr-level_perm-pval.cpp -lgsl -lgslcblas -o get_ftr-level_perm-pval
 *  help2man -o .man ./get_ftr-level_perm-pval
 *  groff -mandoc get_ftr-level_perm-pval.man > get_ftr-level_perm-pval.ps
*/

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>
#include <sys/stat.h>  // for mkdir
#include <sys/types.h> // for mkdir
#include <unistd.h>    // for basename

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

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
       << "  -e, --extr\textract from --geno only the SNPs in cis of the features in --pheno" << endl
       << "\t\tto speed up, tmp files are removed" << endl
       << "      --seed\tseed for the random number generator" << endl
       << "\t\tby default, the RNG is initialized via microseconds from epoch" << endl
       << "  -t, --thread\tnumber of threads (default=1)" << endl
       << "\t\tused for SNPs in cis of the same feature (get_summary_stats)" << endl
       << endl;
}

/** \brief Display version and license information on stdout.
 */
void version (char ** argv)
{
  cout << argv[0] << " 0.1" << endl
       << endl
       << "Copyright (C) 2012 T. Flutre." << endl
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
  bool & extractSnps,
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
	{"trick", no_argument, 0, 0},
	{"extr", no_argument, 0, 'e'},
	{"seed", required_argument, 0, 0},
	{"thread", required_argument, 0, 't'},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:l:f:o:et:",
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
    case 'e':
      extractSnps = true;
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

vector<string>
getListAllInputFeatures (
  const vector<string> & vPhenoFiles,
  const int & verbose)
{
  vector<string> vFtrsToKeep, tokens;
  string line;
  ifstream phenoStream;
  
  if (verbose > 0)
  {
    cout << "retrieve list of features ..." << endl;
    fflush (stdout);
  }
  
  for (vector<string>::const_iterator it = vPhenoFiles.begin();
       it != vPhenoFiles.end(); ++it)
  {
    phenoStream.open(it->c_str());
    if (! phenoStream.is_open())
    {
      cerr << "ERROR: can't open file " << *it << endl;
      exit (1);
    }
    getline (phenoStream, line); // skip header
    while (true)
    {
      getline (phenoStream, line);
      if (line.empty())
	break;
      split (line, ' ', tokens);
      if (find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0])
	  == vFtrsToKeep.end())
	vFtrsToKeep.push_back (tokens[0]);
    }
    phenoStream.close();
  }
  
  if (verbose > 0)
  {
    cout << "nb of features: " << vFtrsToKeep.size() << endl;
    fflush (stdout);
  }
  
  return vFtrsToKeep;
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
  map<string, double> mFtr2TrueL10Abf;
  ifstream truthStream;
  
  truthStream.open(truthFile.c_str());
  if (! truthStream.is_open())
  {
    cerr << "ERROR: can't open file " << truthFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout <<"retrieve true ABFs from file " << truthFile << " ..." << endl;
    fflush (stdout);
  }
  
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
	mFtr2TrueL10Abf.insert (make_pair (tokens[0], atof(tokens[5].c_str())));
      if (whichAbf.compare("abf.fix") == 0)
	mFtr2TrueL10Abf.insert (make_pair (tokens[0], atof(tokens[6].c_str())));
      if (whichAbf.compare("abf.meta.avg.all") == 0)
	mFtr2TrueL10Abf.insert (make_pair (tokens[0], atof(tokens[8].c_str())));
      if (whichAbf.compare("abf.meta.avg.subset") == 0)
	mFtr2TrueL10Abf.insert (make_pair (tokens[0], atof(tokens[9].c_str())));
    }
    ++lineId;
  }
  
  truthStream.close();
  
  if (mFtr2TrueL10Abf.size() == 0)
  {
    cerr << "ERROR: can't retrieve any true ABFs" << endl;
    exit (1);
  }
  
  if (verbose > 0)
  {
    cout << "ABFs retrieved: " << mFtr2TrueL10Abf.size() << endl;
    fflush (stdout);
  }
  
  return mFtr2TrueL10Abf;
}

size_t
getNbSamples (
  const string & genoFile,
  const int & verbose)
{
  size_t nbSamples = 0;
  
  ifstream genoStream;
  genoStream.open(genoFile.c_str());
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  
  string line;
  getline (genoStream, line);
  if (line.empty())
  {
    cerr << "ERROR: file " << genoFile << " is empty" << endl;
    exit (1);
  }
  
  nbSamples = count (line.begin(), line.end(), ' ') - 4;
  if (verbose > 0)
  {
    cout << "nb of samples: " << nbSamples << endl;
    fflush (stdout);
  }
  
  return nbSamples;
}

void
extractLinksForOneFtr (
  const string & linksFile,
  const string & ftrToKeep,
  stringstream & ssTmpLinksFile,
  vector<string> & vCisSnps)
{
  string line, snpCoord;
  vector<string> tokens, tokens2;
  ifstream linksStream;
  ofstream tmpLinksStream;
  
  linksStream.open(linksFile.c_str());
  if (! linksStream.is_open())
  {
    cerr << "ERROR: can't open file " << linksFile << endl;
    exit (1);
  }
  
  ssTmpLinksFile.str("");
  ssTmpLinksFile << "tmp_" << ftrToKeep << "_links.txt";
  tmpLinksStream.open(ssTmpLinksFile.str().c_str());
  if (! linksStream.is_open())
  {
    cerr << "ERROR: can't open file " << linksFile << endl;
    exit (1);
  }
  
  while (true)
  {
    getline (linksStream, line);
    if (line.empty())
      break;
    if (line.find ('\t') != string::npos)
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    if (ftrToKeep.compare(tokens[0]) == 0)
    {
      tmpLinksStream << line << endl;
      split (tokens[1], '|', tokens2);
      vCisSnps.push_back (tokens2[0]);
    }
  }
  
  tmpLinksStream.close();
  linksStream.close();
}

void
extractGenoForOneFtr (
  const string & genoFile,
  const string & ftrToKeep,
  stringstream & ssTmpGenoFile,
  const vector<string> & vCisSnps)
{
  string line;
  vector<string> tokens;
  ifstream genoStream;
  ofstream tmpGenoStream;
  
  genoStream.open(genoFile.c_str());
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  
  ssTmpGenoFile.str("");
  ssTmpGenoFile << "tmp_" << ftrToKeep << "_geno.impute";
  tmpGenoStream.open(ssTmpGenoFile.str().c_str());
  if (! tmpGenoStream.is_open())
  {
    cerr << "ERROR: can't open file " << ssTmpGenoFile.str() << endl;
    exit (1);
  }
  
  // copy header
  getline (genoStream, line);
  tmpGenoStream << line << endl;
  
  while (true)
  {
    getline (genoStream, line);
    if (line.empty())
      break;
    split (line, ' ', tokens);
    if (find(vCisSnps.begin(), vCisSnps.end(), tokens[1])
	!= vCisSnps.end())
      tmpGenoStream << line << endl;
  }
  
  tmpGenoStream.close();
  genoStream.close();
}

void
extractCisSnpsForOneFtr (
  const string & genoFile,
  const string & linksFile,
  const string & ftrToKeep,
  stringstream & ssTmpGenoFile,
  stringstream & ssTmpLinksFile,
  const int & verbose)
{
  vector<string> vCisSnps;
  
  extractLinksForOneFtr (linksFile, ftrToKeep, ssTmpLinksFile, vCisSnps);
  if (verbose > 0)
  {
    cout << "SNPs in cis: " << vCisSnps.size() << endl;
    fflush (stdout);
  }
  
  extractGenoForOneFtr (genoFile, ftrToKeep, ssTmpGenoFile, vCisSnps);
}

void
extractLinks (
  const string & linksFile,
  const vector<string> & vFtrsToKeep,
  const string & studyId,
  string & tmpLinksFile,
  vector<string> & vCisSnps)
{
  ifstream linksStream;
  linksStream.open(linksFile.c_str());
  if (! linksStream.is_open())
  {
    cerr << "ERROR: can't open file " << linksFile << endl;
    exit (1);
  }
  
  stringstream ssTmpLinksFile;
  ssTmpLinksFile << "tmp_" << studyId << "_links.txt";
  tmpLinksFile = ssTmpLinksFile.str();
  ofstream tmpLinksStream;
  tmpLinksStream.open(tmpLinksFile.c_str());
  if (! tmpLinksStream.is_open())
  {
    cerr << "ERROR: can't open file " << tmpLinksFile << endl;
    exit (1);
  }
  
  string line, snpCoord;
  vector<string> tokens, tokens2;
  while (true)
  {
    getline (linksStream, line);
    if (line.empty())
      break;
    if (line.find ('\t') != string::npos)
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    if (find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0])
	!= vFtrsToKeep.end())
    {
      tmpLinksStream << line << endl;
      split (tokens[1], '|', tokens2);
      vCisSnps.push_back (tokens2[0]);
    }
  }
  
  tmpLinksStream.close();
  linksStream.close();
}

void
extractGeno (
  const string & genoFile,
  const vector<string> & vFtrsToKeep,
  const string & studyId,
  string & tmpGenoFile,
  const vector<string> & vCisSnps)
{
  ifstream genoStream;
  genoStream.open(genoFile.c_str());
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }

  stringstream ssTmpGenoFile;
  ssTmpGenoFile << "tmp_" << studyId << "_geno.impute";
  tmpGenoFile = ssTmpGenoFile.str();
  ofstream tmpGenoStream;
  tmpGenoStream.open(tmpGenoFile.c_str());
  if (! tmpGenoStream.is_open())
  {
    cerr << "ERROR: can't open file " << tmpGenoFile << endl;
    exit (1);
  }
  
  // copy header
  string line;
  getline (genoStream, line);
  tmpGenoStream << line << endl;
  
  vector<string> tokens;
  while (true)
  {
    getline (genoStream, line);
    if (line.empty())
      break;
    split (line, ' ', tokens);
    if (find(vCisSnps.begin(), vCisSnps.end(), tokens[1])
	!= vCisSnps.end())
      tmpGenoStream << line << endl;
  }
  
  tmpGenoStream.close();
  genoStream.close();
}

void
extractCisSnps (
  const string & genoFile,
  const string & linksFile,
  const vector<string> & vFtrsToKeep,
  const string & studyId,
  string & tmpGenoFile,
  string & tmpLinksFile,
  const int & verbose)
{
  if (verbose > 0)
  {
    cout << "extract SNPs ..." << endl;
    fflush (stdout);
  }
  
  vector<string> vCisSnps;
  extractLinks (linksFile, vFtrsToKeep, studyId, tmpLinksFile, vCisSnps);
  extractGeno (genoFile, vFtrsToKeep, studyId, tmpGenoFile, vCisSnps);
  
  if (verbose > 0)
  {
    cout << "nb of SNPs: " << vCisSnps.size() << endl;
    fflush (stdout);
  }
}

/** \note TODO: refactor and split in several smaller functions.
 */
void 
computeJointAnalysisPermPvaluesForOneFtr (
  const string & genoFile,
  const vector<string> & vPhenoFiles,
  const string & ftrCoordsFile,
  const string & linksFile,
  const string & gridFile,
  const string & ftrToKeep,
  const double & trueL10Abf,
  const string & whichAbf,
  double & jointPermPval,
  const size_t & nbPermutations,
  const bool & trickPerm,
  gsl_rng * rng,
  const int & nbThreads,
  const int & verbose)
{
  size_t permId = 0;
  gsl_permutation * perm = NULL;
  FILE * pt_permFile = NULL;
  stringstream ssFtrFile, ssTmpGenoFile, ssTmpLinksFile, ssPermFile,
    ssSumStatsDir, ssCmd, ssAbfFile;
  ofstream ftrStream;
  double maxPermL10AbfOverSnps = -1.0;
  
  if (verbose > 0)
  {
    cout << "feature " << ftrToKeep << endl;
    fflush (stdout);
  }
  
  jointPermPval = 0;
  
  ssFtrFile << "tmp_" << ftrToKeep << "_ftr.txt";
  ftrStream.open(ssFtrFile.str().c_str());
  if (! ftrStream.is_open())
  {
    cerr << "ERROR: can't open file " << ssFtrFile.str() << endl;
    exit (1);
  }
  ftrStream << ftrToKeep << endl;
  ftrStream.close();
  
  // to speed-up get_summary_stats
  extractCisSnpsForOneFtr (genoFile, linksFile, ftrToKeep,
			   ssTmpGenoFile, ssTmpLinksFile,
    verbose);
  
  perm = gsl_permutation_calloc (getNbSamples (genoFile, verbose));
  if (perm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the permutation" << endl;
    exit (1);
  }
  
  for (permId = 1; permId <= nbPermutations; ++permId)
  {
    // get permutations and write it into file
    ssPermFile.str("");
    ssPermFile << "tmp_" << ftrToKeep << "_perm" << permId << "_idx.txt";
    pt_permFile = fopen (ssPermFile.str().c_str(), "w");
    if (pt_permFile == NULL)
    {
      cerr << "ERROR: can't open file " << ssPermFile.str() << endl;
      exit (1);
    }
    gsl_ran_shuffle (rng, perm->data, perm->size, sizeof(size_t));
    gsl_permutation_fprintf (pt_permFile, perm, "%i\n");
    fclose (pt_permFile);
    
    // calc sumstats in each subgroup for this permutation
    ssSumStatsDir.str("");
    ssSumStatsDir << "tmp_" << ftrToKeep << "_perm" << permId << "_sumstats2";
    if (removeDir (ssSumStatsDir.str()) != 0)
    {
      cerr << "ERROR: can't remove directory " << ssSumStatsDir.str() << endl;
      exit (1);
    }
    if (mkdir (ssSumStatsDir.str().c_str(), 0755) != 0)
    {
      cerr << "ERROR: can't create directory " << ssSumStatsDir.str() << endl;
      fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
      exit (1);
    }
    for (vector<string>::const_iterator it = vPhenoFiles.begin();
	 it != vPhenoFiles.end(); ++it)
    {
      ssCmd.str("");
      ssCmd << "get_summary_stats"
	    << " -g " << ssTmpGenoFile.str()
	    << " -p " << *it
	    << " -o " << ssSumStatsDir.str() << "/"
	    << basename(it->c_str()) << "_sumstats"
	    << " --fcoord " << ftrCoordsFile
	    << " -l " << ssTmpLinksFile.str()
	    << " -f " << ssFtrFile.str()
	    << " --permf " << ssPermFile.str()
	    << " -v " << verbose - 1;
      if (verbose > 0)
	cout << "launch: " << ssCmd.str() << endl;
      if (system (ssCmd.str().c_str()) != 0)
      {
	cerr << "ERROR: get_summary_stats failed" << endl;
	cerr << ssCmd.str() << endl;
	exit (1);
      }
    }
    
    // calc ABFs from these sumstats
    ssAbfFile.str("");
    ssAbfFile << "tmp_" << ftrToKeep << "_perm" << permId << "_abfs.txt";
    ssCmd.str("");
    ssCmd << "get_abf_meta"
	  << " -i " << ssSumStatsDir.str()
	  << " -g " << gridFile
	  << " -o " << ssAbfFile.str()
	  << " -v " << verbose - 1;
    if (whichAbf.find("avg") != string::npos)
      ssCmd << " -c";
    if (verbose > 0)
      cout << "launch: " << ssCmd.str() << endl;
    if (system (ssCmd.str().c_str()) != 0)
    {
      cerr << "ERROR: get_abf_meta failed" << endl;
      cerr << ssCmd.str() << endl;
      exit (1);
    }
    
    // retrieve the max log10(ABF) across SNPs
    ifstream abfStream;
    abfStream.open (ssAbfFile.str().c_str());
    if (! abfStream.is_open())
    {
      cerr << "ERROR: can't open file " << ssAbfFile.str() << endl;
      exit (1);
    }
    string line;
    getline (abfStream, line); // skip header
    if (line.empty())
    {
      cerr << "ERROR: file " << ssAbfFile.str() << " is empty" << endl;
      exit (1);
    }
    vector<string> tokens = split (line, ' ');
    while (abfStream.good())
    {
      getline (abfStream, line);
      if (line.empty())
	break;
      split (line, ' ', tokens);
      if (whichAbf.compare("abf.meta") == 0
	  && atof(tokens[5].c_str()) > maxPermL10AbfOverSnps)
	maxPermL10AbfOverSnps = atof(tokens[5].c_str());
      else if (whichAbf.compare("abf.fix") == 0
	  && atof(tokens[6].c_str()) > maxPermL10AbfOverSnps)
	maxPermL10AbfOverSnps = atof(tokens[6].c_str());
      else if (whichAbf.compare("abf.meta.avg.all") == 0
	  && atof(tokens[8].c_str()) > maxPermL10AbfOverSnps)
	maxPermL10AbfOverSnps = atof(tokens[8].c_str());
      else if (whichAbf.compare("abf.meta.avg.subset") == 0
	  && atof(tokens[9].c_str()) > maxPermL10AbfOverSnps)
	maxPermL10AbfOverSnps = atof(tokens[9].c_str());
    }
    abfStream.close();
    if (maxPermL10AbfOverSnps >= trueL10Abf)
      ++jointPermPval;
    
    // clean
    if (remove (ssPermFile.str().c_str()) != 0)
    {
      cerr << "ERROR: can't remove file " << ssPermFile.str() << endl;
      exit (1);
    }
    if (removeDir (ssSumStatsDir.str()) != 0)
    {
      cerr << "ERROR: can't remove directory " << ssSumStatsDir.str() << endl;
      exit (1);
    }
    if (remove (ssAbfFile.str().c_str()) != 0)
    {
      cerr << "ERROR: can't remove file " << ssAbfFile.str() << endl;
      exit (1);
    }
    
    if (trickPerm && jointPermPval == 11)
      break;
  }
  
  if (! trickPerm)
    jointPermPval /= (nbPermutations + 1);
  else
    jointPermPval = gsl_ran_flat (rng, 11 / ((double) (permId + 2)),
				  11 / ((double) (permId + 1)));
  
  // clean
  if (remove (ssTmpGenoFile.str().c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << ssTmpGenoFile.str() << endl;
    exit (1);
  }
  if (remove (ssTmpLinksFile.str().c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << ssTmpLinksFile.str() << endl;
    exit (1);
  }
  if (remove (ssFtrFile.str().c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << ssFtrFile.str() << endl;
    exit (1);
  }
  
  gsl_permutation_free (perm);
}

void computeJointAnalysisPermPvaluesFtrPerFtr (
  const string & genoFile,
  const vector<string> & vPhenoFiles,
  const string & ftrCoordsFile,
  const string & linksFile,
  const string & gridFile,
  const vector<string> & vFtrsToKeep,
  const map<string, double> & mFtr2TrueL10Abf,
  const string & whichAbf,
  const string & outFile,
  const size_t & nbPermutations,
  const bool & trickPerm,
  const int & seed,
  const int & nbThreads,
  const int & verbose)
{
  vector<string>::const_iterator it;
  double trueL10Abf, jointPermPval;
  ofstream outStream;
  gsl_rng * rng = NULL;
  
  if (verbose > 0)
  {
    cout << "compute joint-analysis permutation P-values"
	 << " (" << vFtrsToKeep.size()
	 << " features) ..." << endl;
    fflush (stdout);
  }
  
  outStream.open(outFile.c_str());
  if (! outStream.is_open())
  {
    cerr << "ERROR: can't open file " << outFile << endl;
    exit (1);
  }
  outStream << "ftr jointPermPval" << endl;
  
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  if (rng == NULL)
  {
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit (1);
  }
  if (seed < 0)
    gsl_rng_set (rng, getSeed());
  else
    gsl_rng_set (rng, seed);
  
  for (it = vFtrsToKeep.begin(); it != vFtrsToKeep.end(); ++it)
  {
    computeJointAnalysisPermPvaluesForOneFtr (genoFile,
					      vPhenoFiles,
					      ftrCoordsFile,
					      linksFile,
					      gridFile,
					      *it,
					      trueL10Abf,
					      whichAbf,
					      jointPermPval,
					      nbPermutations,
					      trickPerm,
					      rng,
					      nbThreads,
					      verbose-1);
    outStream << *it << " " << jointPermPval << endl;
  }
  
  outStream.close();
  gsl_rng_free (rng);
  if (verbose > 0)
    cout << "results written in " << outFile << endl;
}

void 
makeNewPermutationAndWriteToFile (
  const string & studyId,
  const size_t & permId,
  string & permFile,
  gsl_rng * rng,
  gsl_permutation * perm)
{
  stringstream ssPermFile;
  ssPermFile << "tmp_" << studyId << "_perm" << permId << "_idx.txt";
  permFile = ssPermFile.str();
  FILE * pt_permFile = fopen (ssPermFile.str().c_str(), "w");
  if (pt_permFile == NULL)
  {
    cerr << "ERROR: can't open file " << ssPermFile.str() << endl;
    exit (1);
  }
  gsl_ran_shuffle (rng, perm->data, perm->size, sizeof(size_t));
  gsl_permutation_fprintf (pt_permFile, perm, "%i\n");
  fclose (pt_permFile);
}

void writeFeaturesToAnalyze (
  const string & studyId,
  const size_t & permId,
  const map<string, double> & mFtr2MaxPermAbf,
  string & tmpFtrFile)
{
  stringstream ssTmpFtrFile;
  ssTmpFtrFile << "tmp_" << studyId << "_perm" << permId << "_ftrs.txt";
  tmpFtrFile = ssTmpFtrFile.str();
  
  ofstream tmpFtrStream;
  tmpFtrStream.open(tmpFtrFile.c_str());
  if (! tmpFtrStream.is_open())
  {
    cerr << "ERROR: can't open file " << tmpFtrFile << endl;
    exit (1);
  }
  
  for (map<string, double>::const_iterator it = mFtr2MaxPermAbf.begin();
       it != mFtr2MaxPermAbf.end(); ++it)
    tmpFtrStream << it->first << endl;
  
  tmpFtrStream.close();
}

void getSumStatsInEachSubgroup (
  const string & studyId,
  const size_t & permId,
  const string & permFile,
  const string & genoFile,
  const vector<string> & vPhenoFiles,
  const string & ftrCoordsFile,
  const string & linksFile,
  const string & ftrFile,
  string & sumStatsDir,
  string & stdoutFile,
  const int & verbose)
{
  stringstream ssSumStatsDir;
  ssSumStatsDir << "tmp_" << studyId << "_perm" << permId << "_sumstats";
  sumStatsDir = ssSumStatsDir.str();
  if (removeDir (sumStatsDir) != 0)
  {
    cerr << "ERROR: can't remove directory " << sumStatsDir << endl;
    exit (1);
  }
  if (mkdir (sumStatsDir.c_str(), 0755) != 0)
  {
    cerr << "ERROR: can't create directory " << sumStatsDir << endl;
    fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
    exit (1);
  }
  
  stringstream ssStdoutFile;
  ssStdoutFile << "tmp_" << studyId << "_perm" << permId << "_stdout.txt";
  stdoutFile = ssStdoutFile.str();
  if (doesFileExist (stdoutFile))
    remove (stdoutFile.c_str());
  
  for (vector<string>::const_iterator it = vPhenoFiles.begin();
       it != vPhenoFiles.end(); ++it)
  {
    stringstream ssCmd, ssStdoutFile;
    ssCmd << "get_summary_stats"
	  << " -g " << genoFile
	  << " -p " << *it
	  << " -o " << sumStatsDir << "/"
	  << basename(it->c_str()) << "_sumstats"
	  << " --fcoord " << ftrCoordsFile
	  << " -l " << linksFile;
    if (! ftrFile.empty())
      ssCmd << " -f " << ftrFile;
    ssCmd << " --permf " << permFile
	  << " >> " << stdoutFile
	  << " 2>> " << stdoutFile;
    if (verbose > 0)
      cout << "launch: " << ssCmd.str() << endl;
    if (system (ssCmd.str().c_str()) != 0)
    {
      cerr << "ERROR: get_summary_stats failed" << endl;
      cerr << ssCmd.str() << endl;
      exit (1);
    }
  }
}

void getAbfs (
  const string & studyId,
  const size_t & permId,
  string & abfsFile,
  const string & sumStatsDir,
  const string & gridFile,
  const string & whichAbf,
  const string & stdoutFile,
  const int & verbose)
{
  stringstream ssAbfFile, ssCmd;
  ssAbfFile << "tmp_" << studyId << "_perm" << permId << "_abfs.txt";
  abfsFile = ssAbfFile.str();
  ssCmd << "get_abf_meta"
	<< " -i " << sumStatsDir
	<< " -g " << gridFile
	<< " -o " << abfsFile;
  if (whichAbf.find("avg") != string::npos)
    ssCmd << " -c";
  ssCmd << " >> " << stdoutFile
	<< " 2>> " << stdoutFile;
  if (verbose > 0)
    cout << "launch: " << ssCmd.str() << endl;
  if (system (ssCmd.str().c_str()) != 0)
  {
    cerr << "ERROR: get_abf_meta failed" << endl;
    cerr << ssCmd.str() << endl;
    exit (1);
  }
}

void getMaxAbfAcrossSnps (
  const string & abfsFile,
  const string & whichAbf,
  map<string, double> & mFtr2MaxPermAbf)
{
  ifstream abfStream;
  abfStream.open (abfsFile.c_str());
  if (! abfStream.is_open())
  {
    cerr << "ERROR: can't open file " << abfsFile << endl;
    exit (1);
  }
  
  string line, ftr;
  getline (abfStream, line); // skip header
  if (line.empty())
  {
    cerr << "ERROR: file " << abfsFile << " is empty" << endl;
    exit (1);
  }
  
  vector<string> tokens = split (line, ' ');
  while (abfStream.good())
  {
    getline (abfStream, line);
    if (line.empty())
      break;
    split (line, ' ', tokens);
    ftr = tokens[0];
    if (whichAbf.compare("abf.meta") == 0
	&& atof(tokens[5].c_str()) > mFtr2MaxPermAbf[ftr])
      mFtr2MaxPermAbf[ftr] = atof(tokens[5].c_str());
    else if (whichAbf.compare("abf.fix") == 0
	     && atof(tokens[6].c_str()) > mFtr2MaxPermAbf[ftr])
      mFtr2MaxPermAbf[ftr] = atof(tokens[6].c_str());
    else if (whichAbf.compare("abf.meta.avg.all") == 0
	     && atof(tokens[8].c_str()) > mFtr2MaxPermAbf[ftr])
      mFtr2MaxPermAbf[ftr] = atof(tokens[8].c_str());
    else if (whichAbf.compare("abf.meta.avg.subset") == 0
	     && atof(tokens[9].c_str()) > mFtr2MaxPermAbf[ftr])
      mFtr2MaxPermAbf[ftr] = atof(tokens[9].c_str());
  }
  
  abfStream.close();
}

void
cleanPermutationTmp (
  const string & permFile,
  const string & tmpFtrFile,
  const string & sumStatsDir,
  const string & abfsFile,
  const string & stdoutFile)
{
  if (remove (permFile.c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << permFile << endl;
    exit (1);
  }
  if (remove (tmpFtrFile.c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << tmpFtrFile << endl;
    exit (1);
  }
  if (removeDir (sumStatsDir) != 0)
  {
    cerr << "ERROR: can't remove directory " << sumStatsDir << endl;
    exit (1);
  }
  if (remove (abfsFile.c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << abfsFile << endl;
    exit (1);
  }
  if (remove (stdoutFile.c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << stdoutFile << endl;
    exit (1);
  }
}

/** \brief http://stackoverflow.com/a/263958/597069
 */
void updateFtrsToAnalyzeAndCalcPermPval (
  const map<string, double> & mFtr2TrueL10Abf,
  const size_t & permId,
  map<string, double> & mFtr2MaxPermAbf,
  map<string, vector<double> > & mFtr2PermPval,
  const bool & trick)
{
  string ftr;
  double maxTrueL10Abf;
  map<string, double>::iterator it = mFtr2MaxPermAbf.begin();
  while (it != mFtr2MaxPermAbf.end())
  {
    ftr = it->first;
    mFtr2PermPval[ftr][0] = permId;
    maxTrueL10Abf = mFtr2TrueL10Abf.find(ftr)->second;
    if (it->second >= maxTrueL10Abf)
      ++(mFtr2PermPval[ftr][1]);
    if (trick && mFtr2PermPval[ftr][1] == 10)
      mFtr2MaxPermAbf.erase (it++); // note the post increment
    else
      ++it;
  }
}

void computeJointAnalysisPermPvaluesFtrPerFtr2 (
  const string & genoFile,
  const vector<string> & vPhenoFiles,
  const string & ftrCoordsFile,
  const string & linksFile,
  const string & gridFile,
  const string & ftrFile,
  const vector<string> & vFtrsToKeep,
  const map<string, double> & mFtr2TrueL10Abf,
  const string & whichAbf,
  const string & outFile,
  const size_t & nbPermutations,
  const bool & trickPerm,
  const bool & extractSnps,
  const int & seed,
  const int & nbThreads,
  const int & verbose)
{
  if (verbose > 0)
  {
    cout << "compute joint-analysis permutation P-values"
	 << " (" << vFtrsToKeep.size()
	 << " features) ..." << endl;
    fflush (stdout);
  }
  
  ofstream outStream;
  outStream.open(outFile.c_str());
  if (! outStream.is_open())
  {
    cerr << "ERROR: can't open file " << outFile << endl;
    exit (1);
  }
  outStream << "ftr jointPermPval" << endl;
  
  gsl_rng_env_setup();
  gsl_rng * rng = gsl_rng_alloc (gsl_rng_default);
  if (rng == NULL)
  {
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit (1);
  }
  if (seed < 0)
    gsl_rng_set (rng, getSeed());
  else
    gsl_rng_set (rng, seed);
  
  gsl_permutation * perm = gsl_permutation_calloc (getNbSamples (genoFile,
								 verbose));
  if (perm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the permutation" << endl;
    exit (1);
  }
  
  map<string, double> mFtr2MaxPermAbf;
  map<string, vector<double> > mFtr2PermPval;
  for (vector<string>::const_iterator it = vFtrsToKeep.begin();
       it != vFtrsToKeep.end(); ++it)
  {
    mFtr2MaxPermAbf.insert (make_pair (*it, -1.0));
    
    // 0: nb performed permutations, 1: nb with perm ABF >= true ABF
    vector<double> vStats (2, 0);
    mFtr2PermPval.insert (make_pair (*it, vStats));
  }
  
  stringstream ssStudyId;
  ssStudyId << getSeed() << "-" << getpid();
  string tmpLinksFile, tmpGenoFile;
  if (extractSnps)
    extractCisSnps (genoFile, linksFile, vFtrsToKeep, ssStudyId.str(),
		    tmpGenoFile, tmpLinksFile, verbose);
  
  string permFile, tmpFtrFile, stdoutFile, sumStatsDir, abfsFile;
  for (size_t permId = 1; permId <= nbPermutations; ++permId)
  {
    makeNewPermutationAndWriteToFile (ssStudyId.str(), permId, permFile, rng,
				      perm);
    writeFeaturesToAnalyze (ssStudyId.str(), permId, mFtr2MaxPermAbf,
			    tmpFtrFile);
    getSumStatsInEachSubgroup (ssStudyId.str(), permId, permFile,
			       (extractSnps ? tmpGenoFile : genoFile),
			       vPhenoFiles, ftrCoordsFile,
			       (extractSnps ? tmpLinksFile : linksFile),
			       tmpFtrFile, sumStatsDir, stdoutFile, verbose-1);
    getAbfs (ssStudyId.str(), permId, abfsFile, sumStatsDir, gridFile,
	     whichAbf, stdoutFile, verbose-1);
    getMaxAbfAcrossSnps (abfsFile, whichAbf, mFtr2MaxPermAbf);
    cleanPermutationTmp (permFile, tmpFtrFile, sumStatsDir, abfsFile, stdoutFile);
    updateFtrsToAnalyzeAndCalcPermPval (mFtr2TrueL10Abf, permId,
					mFtr2MaxPermAbf, mFtr2PermPval,
					trickPerm);
    if (mFtr2MaxPermAbf.size() == 0)
      break;
  }
  
  for (map<string, vector<double> >::iterator it = mFtr2PermPval.begin();
       it != mFtr2PermPval.end(); ++it)
  {
    outStream << it->first << " ";
    if (! trickPerm || it->second[0] == nbPermutations)
      outStream << ((it->second[1] + 1) /
		    ((double) it->second[0] + 1));
    else
      outStream << gsl_ran_flat (rng,
				 ((it->second[1] + 1) /
				  ((double) (it->second[0] + 2))),
				 ((it->second[1] + 1) /
				  ((double) (it->second[0] + 1))));
    outStream << endl;
  }
  outStream.close();
  if (verbose > 0)
    cout << "results written in " << outFile << endl;
  
  gsl_permutation_free (perm);
  gsl_rng_free (rng);
  if (remove (tmpGenoFile.c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << tmpGenoFile << endl;
    exit (1);
  }
  if (remove (tmpLinksFile.c_str()) != 0)
  {
    cerr << "ERROR: can't remove file " << tmpLinksFile << endl;
    exit (1);
  }
}

int main (int argc, char ** argv)
{
  string genoFile, phenoDir, ftrCoordsFile, linksFile, gridFile, ftrFile,
    truthFile, whichAbf = "abf.meta", outFile;
  size_t nbPermutations = 10000;
  bool trickPerm = false, extractSnps = false;
  int verbose = 1, seed = -1, nbThreads = 1;
  
  parse_args (argc, argv, genoFile, phenoDir, ftrCoordsFile, linksFile,
	      gridFile, ftrFile, truthFile, whichAbf, outFile, nbPermutations,
	      trickPerm, extractSnps, seed, nbThreads, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl;
  }
  
  vector<string> vPhenoFiles = scanInputDirectory (phenoDir, verbose);
  vector<string> vFtrsToKeep;
  if (! ftrFile.empty())
    vFtrsToKeep = loadOneColumnFile (ftrFile, verbose);
  else
    vFtrsToKeep = getListAllInputFeatures (vPhenoFiles, verbose);
  map<string, double> mFtr2TrueL10Abf = getTrueAbfForEachFtr (truthFile,
							      whichAbf,
							      vFtrsToKeep,
							      verbose);
  
  computeJointAnalysisPermPvaluesFtrPerFtr2 (genoFile,
					    vPhenoFiles,
					    ftrCoordsFile,
					    linksFile,
					    gridFile,
					     ftrFile,
					    vFtrsToKeep,
					    mFtr2TrueL10Abf,
					    whichAbf,
					    outFile,
					    nbPermutations,
					    trickPerm,
					    extractSnps,
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
