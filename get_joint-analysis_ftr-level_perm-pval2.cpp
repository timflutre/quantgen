/** \file get_joint-analysis_ftr-level_perm-pval2.cpp
 *
 *  `get_joint-analysis_ftr-level_perm-pval2' computes feature-level permutation P-values.
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
 *  g++ -Wall -O3 -fopenmp utils.cpp get_joint-analysis_ftr-level_perm-pval2.cpp -lgsl -lgslcblas -o get_joint-analysis_ftr-level_perm-pval2
 *  help2man -o get_joint-analysis_ftr-level_perm-pval2.man ./get_joint-analysis_ftr-level_perm-pval2
 *  groff -mandoc get_joint-analysis_ftr-level_perm-pval2.man > get_joint-analysis_ftr-level_perm-pval2.ps
*/

#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include "utils.h"
#include "get_summary_stats2.cpp"
#include "get_abf_meta.cpp"

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
       << "  -g, --geno\tfile with absolute paths to genotype files" << endl
       << "\t\tcan be a single line (eg. for multiple tissues)" << endl
       << "\t\tshould be in IMPUTE format (delimiter: space or tab)" << endl
       << "\t\ta header line with sample names is required" << endl
       << "  -p, --pheno\tfile with absolute paths to phenotype files" << endl
       << "\t\tcan be a single line (single subgroup)" << endl
       << "\t\trow 1 for sample names, column 1 for feature names" << endl
       << "\t\tsubgroups can have different features" << endl
       << "\t\tall features should be in the --fcoord file" << endl
       << "      --fcoord\tfile with the features coordinates" << endl
       << "\t\tshould be in the BED format (delimiter: tab)" << endl
       << "      --grid\tfile with the grid of values for phi2 and omega2 (ES model)" << endl
       << "\t\tsee GetGridPhiOmega() in package Rquantgen" << endl
       << "      --abf\twhich ABF to use as the test statistic for the permutations" << endl
       << "\t\tdefault=abf.const/abf.fix/abf.avg.all/abf.avg.subset" << endl
       << "  -o, --out\tprefix for the output files" << endl
       << "  -P, --perm\tnumber of permutations" << endl
       << "\t\tthe phenotype labels are permuted" << endl
       << "\t\tdefault=10000" << endl
       << "      --seed\tseed for the two random number generators" << endl
       << "\t\tone for the permutations, another for the trick" << endl
       << "\t\tby default, both are initialized via microseconds from epoch" << endl
       << "  -t, --trick\tapply trick to speed-up permutations" << endl
       << "\t\tstop after the tenth permutation for which the test statistic" << endl
       << "\t\tis better than or equal to the true value, and sample from" << endl
       << "\t\ta uniform between 11/(nbPermsSoFar+2) and 11/(nbPermsSoFar+1)" << endl
       << "\t\tif '1', the permutations really stops" << endl
       << "\t\tif '2', all permutations are done but the test statistics are not computed" << endl
       << "\t\tallow to compare different test statistics on the same permutations" << endl
       << "  -f, --ftr\tfile with a list of features to analyze" << endl
       << "\t\tone feature name per line" << endl
       << "  -s, --snp\tfile with a list of SNPs to analyze" << endl
       << "\t\tone SNP name per line" << endl
       << "  -q, --qnorm\tquantile-normalize the phenotypes" << endl
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
parseArgs (
  int argc,
  char ** argv,
  string & genoPathsFile,
  string & phenoPathsFile,
  string & ftrCoordsFile,
  string & gridFile,
  string & whichBf,
  string & outPrefix,
  size_t & nbPerms,
  size_t & seed,
  int & trick,
  string & ftrsToKeepFile,
  string & snpsToKeepFile,
  bool & needQnorm,
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
	{"grid", required_argument, 0, 0},
	{"abf", required_argument, 0, 0},
	{"out", required_argument, 0, 'o'},
	{"perm", required_argument, 0, 'P'},
	{"seed", required_argument, 0, 0},
	{"trick", required_argument, 0, 't'},
	{"ftr", required_argument, 0, 'f'},
	{"snp", required_argument, 0, 's'},
	{"qnorm", no_argument, 0, 'q'},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:o:P:t:f:s:q",
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
      if (strcmp(long_options[option_index].name, "abf") == 0)
      {
	whichBf = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "seed") == 0)
      {
	seed = atol (optarg);
	break;
      }
    case 'h':
      help (argv);
      exit (0);
    case 'V':
      version (argv);
      exit (0);
    case 'v':
      verbose = atoi (optarg);
      break;
    case 'g':
      genoPathsFile = optarg;
      break;
    case 'p':
      phenoPathsFile = optarg;
      break;
    case 'o':
      outPrefix = optarg;
      break;
    case 'P':
      nbPerms = atol (optarg);
      break;
    case 't':
      trick = atoi (optarg);
      break;
    case 'f':
      ftrsToKeepFile = optarg;
      break;
    case 's':
      snpsToKeepFile = optarg;
      break;
    case 'q':
      needQnorm = true;
      break;
    case '?':
      printf ("\n"); help (argv);
      abort ();
    default:
      printf ("\n"); help (argv);
      abort ();
    }
  }
  if (genoPathsFile.empty())
  {
    fprintf (stderr, "ERROR: missing compulsory option -g\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (genoPathsFile))
  {
    fprintf (stderr, "ERROR: can't file '%s'\n\n", genoPathsFile.c_str());
    help (argv);
    exit (1);
  }
  if (phenoPathsFile.empty())
  {
    fprintf (stderr, "ERROR: missing compulsory option -p\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (phenoPathsFile))
  {
    fprintf (stderr, "ERROR: can't find '%s'\n\n", phenoPathsFile.c_str());
    help (argv);
    exit (1);
  }
  if (ftrCoordsFile.empty())
  {
    fprintf (stderr, "ERROR: missing compulsory options --fcoord\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (ftrCoordsFile))
  {
    fprintf (stderr, "ERROR: can't find '%s'\n\n", ftrCoordsFile.c_str());
    help (argv);
    exit (1);
  }
  if (gridFile.empty())
  {
    fprintf (stderr, "ERROR: missing compulsory options --grid\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (gridFile))
  {
    fprintf (stderr, "ERROR: can't find '%s'\n\n", gridFile.c_str());
    help (argv);
    exit (1);
  }
  if (outPrefix.empty())
  {
    fprintf (stderr, "ERROR: missing compulsory option -o\n\n");
    help (argv);
    exit (1);
  }
  if (whichBf.compare("abf.const") == 0 && whichBf.compare("abf.fix") == 0
      && whichBf.compare("abf.avg.subset") == 0
      && whichBf.compare("abf.avg.all") == 0)
  {
    fprintf (stderr, "ERROR: --abf should be abf.const, abf.fix, abf.avg.subset or abf.avg.all\n\n");
    help (argv);
    exit (1);
  }
  if (trick != 0 && trick != 1 && trick != 2)
  {
    fprintf (stderr, "ERROR: unrecognized --trick %d\n\n", trick);
    help (argv);
    exit (1);
  }
  if (! ftrsToKeepFile.empty() && ! doesFileExist (ftrsToKeepFile))
  {
    fprintf (stderr, "ERROR: can't find '%s'\n\n", ftrsToKeepFile.c_str());
    help (argv);
    exit (1);
  }
  if (! snpsToKeepFile.empty() && ! doesFileExist (snpsToKeepFile))
  {
    fprintf (stderr, "ERROR: can't find '%s'\n\n", snpsToKeepFile.c_str());
    help (argv);
    exit (1);
  }
  if (seed == string::npos)
    seed = getSeed ();
}

struct sFtr
{
  string name; // eg. ENSG00000182816
  string chr; // eg. chr21
  size_t start; // 1-based coordinate
  size_t end; // idem
  vector< vector<double> > vvPhenos; // phenotypes of samples per subgroup
  vector< vector<bool> > vvIsNa; // missing values per subgroup
  double maxL10TrueAbf; // max log10(ABF) over SNPs on non-permuted data
  double permPval; // permutation P-value of the joint analysis
  size_t nbPermsSoFar;
};

void
loadSamples (
  const vector<string> & vGenoPaths,
  const vector<string> & vPhenoPaths,
  vector<string> & vSamples,
  vector<vector<size_t> > vvSampleIdxs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load samples ..." << endl << flush;
  
  // load sample names from all subgroups
  vector<vector<string> > vvSamples;
  ifstream phenoStream;
  string line;
  for (size_t s = 0; s < vPhenoPaths.size(); ++s)
  {
    openFile (vPhenoPaths[s], phenoStream);
    getline (phenoStream, line);
    phenoStream.close();
    if (s == 0)
    {
      split (line, " \t", vSamples);
      if (vSamples[0].compare("Id") == 0)
	vSamples.erase (vSamples.begin());
      vvSamples.push_back (vSamples);
    }
    else
    {
      vector<string> tokens;
      split (line, " \t", tokens);
      if (tokens[0].compare("Id") == 0)
	tokens.erase (tokens.begin());
      vvSamples.push_back (tokens);
      for (size_t i = 0; i < tokens.size(); ++i)
	if (find (vSamples.begin(), vSamples.end(), tokens[i])
	    == vSamples.end())
	  vSamples.push_back (tokens[i]);
    }
  }
  if (verbose > 0)
  {
    cout << "total nb of samples: " << vSamples.size() << endl << flush;
    for (size_t s = 0; s < vPhenoPaths.size(); ++s)
      cout << "s" << (s+1) << " (" << vPhenoPaths[s] << "): "
	   << vvSamples[s].size() << " samples" << endl << flush;
  }
  
  // determine index of each sample in the global vector
  for (size_t s = 0; s < vvSamples.size(); ++s)
  {
    vector<size_t> vSampleIdxs (vvSamples[s].size(), string::npos);
    for (size_t i = 0; i < vvSamples[s].size(); ++i)
      vSampleIdxs[i] = find(vSamples.begin(), vSamples.end(), vvSamples[s][i])
	- vSamples.begin();
    vvSampleIdxs.push_back (vSampleIdxs);
//    for (size_t i = 0; i < vvSamples[s].size(); ++i)
//      cout << (s+1) << " " << (i+1) << " " << (vvSampleIdxs[s][i]+1) << endl;
  }
  
  // check that genotypes are available for each sample
  ifstream genoStream;
  vector<string> tokens;
  size_t nbSamplesWithGeno;
  for (size_t s = 0; s < vGenoPaths.size(); ++s)
  {
    openFile (vGenoPaths[s], genoStream);
    getline (genoStream, line);
    genoStream.close();
    split (line, " \t", tokens);
    nbSamplesWithGeno = 0;
    for (size_t i = 5; i < tokens.size(); ++i)
      if (find(vSamples.begin(), vSamples.end(), split (tokens[i], "_", 0))
	  != vSamples.end())
	++nbSamplesWithGeno;
    if (nbSamplesWithGeno != 3 * vSamples.size())
    {
      cerr << "ERROR: samples have phenotypes but no genotypes" << endl;
      exit (1);
    }
  }
}
/*
void
loadPhenos (
  const vector<string> & vPhenoPaths,
  const vector<string> & vFtrsToKeep,
  map<string, sFtr> & mFtrs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load phenotypes ..." << endl << flush;
  
  ifstream phenoStream;
  string line;
  vector<string> tokens;
  
  for (size_t s = 0; s < vPhenoPaths.size(); ++s)
  {
    openFile (vPhenoPaths[s], phenoStream);
    while (true)
    {
      getline (phenoStream, line);
      if (line.empty())
	break;
      split (line, " \t", tokens);
      if (! vFtrsToKeep.empty()
	  && find (vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0])
	  == vFtrsToKeep.end())
	continue;
      
      if (mFeatures.find(tokens[0]) == mFeatures.end())
      {
	if (s != 0)
	{
	  cerr << "ERROR: all phenotype files don't have the same features"
	       << endl;
	  exit (1);
	}
	FtrMultiS iFtr;
	FtrMultiS_init (iFtr, tokens[0], vPhenoPaths.size(), vSamples.size());
	for (size_t i = 1; i < tokens.size(); ++i)
	{
	  if (tokens[i].compare("NA") == 0)
	    iFtr.vvIsNa[s][i-1] = true;
	  iFtr.vvPhenos[s][i-1] = atof (tokens[i].c_str());
	}
	mFeatures.insert (make_pair (tokens[0], iFtr));
      }
      else
      {
	for (size_t i = 1; i < tokens.size() ; ++i)
	{
	  if (tokens[i].compare("NA") == 0)
	    mFeatures[tokens[0]].vvIsNa[s][i-1] = true;
	  mFeatures[tokens[0]].vvPhenos[s][i-1] = atof (tokens[i].c_str());
	}
      }
    }
    
    phenoStream.close();
  }
  
  if (mFeatures.size() == 0)
  {
    cerr << "ERROR: no feature to analyze" << endl;
    exit (1);
  }
  if (verbose > 0)
    cout << "nb of features: " << mFeatures.size() << endl;
}

void
loadFtrInfo (
  const string & ftrCoordsFile,
  map<string, sFtr> mFtrs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load feature coordinates ..." << endl << flush;
  
  ifstream ftrCoordsStream;
  openFile (ftrCoordsFile, ftrCoordsStream);
  
  string line;
  vector<string> tokens;
  size_t countFtrs = 0;
  while (true)
  {
    getline (ftrCoordsStream, line);
    if (line.empty())
      break;
    split (line, " \t", tokens);
    if (mFeatures.find(tokens[3]) == mFeatures.end())
      continue;
    ++countFtrs;
    mFeatures[tokens[3]].chr = tokens[0];
    mFeatures[tokens[3]].start = atol (tokens[1].c_str()) + 1;
    mFeatures[tokens[3]].end = atol (tokens[2].c_str());
  }
  
  ftrCoordsStream.close();
  
  if (countFtrs < mFeatures.size())
  {
    cerr << "ERROR: " << mFeatures.size() - countFtrs
	 << " feature coordinates are missing" << endl;
    exit (1);
  }
}

void
Snp_init (
  Snp & iSnp,
  const string & name,
  const size_t & nbSubgroups,
  const size_t & nbSamples)
{
  iSnp.name = name;
  iSnp.vGenos.assign (nbSamples, 0.0);
  iSnp.vIsNa.assign (nbSamples, false);
}

void
FtrMultiS_init (
  FtrMultiS & iFtr,
  const string & name,
  const size_t & nbSubgroups,
  const size_t & nbSamples)
{
  iFtr.name = name;
  for (size_t s = 0; s < nbSubgroups; ++s)
  {
    iFtr.vvPhenos.push_back (vector<double> (nbSamples, 0.0));
    iFtr.vvIsNa.push_back (vector<bool> (nbSamples, false));
  }
  iFtr.maxL10TrueAbf = 0.0;
  iFtr.permPval = 1.0;
  iFtr.nbPermsSoFar = 0;
}

void
loadLinksFtrCisSnpsFromFile (
  map<string, FtrMultiS> & mFeatures,
  map<string, Snp> & mSnps,
  const string & linksFile,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load links features - SNPs in cis ..." << endl << flush;
  
  ifstream linksStream;
  openFile (linksFile, linksStream);
  
  string line;
  vector<string> tokens, tokens2;
  size_t countLines = 0;
  while (true)
  {
    getline (linksStream, line);
    if (line.empty())
      break;
    ++countLines;
    split (line, " \t", tokens);
    if (mFeatures.find(tokens[0]) == mFeatures.end())
      continue;
    split (tokens[1], '|', tokens2); // SNPid|SNPcoord
    if (mSnps.find(tokens2[0]) == mSnps.end())
    {
      Snp iSnp;
      Snp_init (iSnp, tokens2[0], mFeatures[tokens[0]].vvPhenos.size(),
		mFeatures[tokens[0]].vvPhenos[0].size());
      mSnps.insert (make_pair (iSnp.name, iSnp));
    }
    mFeatures[tokens[0]].vCisSnps.push_back (&(mSnps[tokens2[0]]));
  }
  
  linksStream.close();
  
  map<string, FtrMultiS>::iterator it = mFeatures.begin();
  while (it != mFeatures.end())
  {
    if (it->second.vCisSnps.size() == 0)
      mFeatures.erase (it++);
    else++it;
  }
  
  if (verbose > 0)
    cout << mSnps.size() << " SNPs in cis of " << mFeatures.size () << " features" << endl;
}

void
loadGenoSnpsFromFile (
  map<string, Snp> & mSnps,
  const string & genoFile,
  const vector<string> & vSamples,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load genotype values ..." << endl << flush;
  
  ifstream genoStream;
  openFile (genoFile, genoStream);
  
  // check header line
  string line;
  vector<string> tokens;
  getline (genoStream, line);
  if (line.empty())
  {
    cerr << "ERROR: file " << genoFile << " is empty" << endl;
    exit (1);
  }
  split (line, " \t", tokens);
  if (tokens.size() -5 != vSamples.size())
  {
    cerr << "ERROR: different number of samples between genotypes ("
	 << (tokens.size() - 5) << ") and phenotypes (" << vSamples.size()
	 << ")" << endl;
    exit (1);
  }
  for (size_t i = 0; i < vSamples.size(); ++i)
    if (tokens[5+i].compare(vSamples[i]) != 0)
    {
      cerr << "ERROR: different samples between genotype and phenotype files"
	   << endl << "or they are not in the same order" << endl;
      exit (1);
    }
  
  // load genotype values
  size_t countSnps = 0;
  double maf, AA, AB, BB;
  while (true)
  {
    getline (genoStream, line);
    if (line.empty())
      break;
    split (line, " \t", tokens);
    if (mSnps.find(tokens[1]) == mSnps.end())
      continue;
    ++countSnps;
    mSnps[tokens[1]].chr = tokens[0];
    mSnps[tokens[1]].coord = atol (tokens[2].c_str());
    maf = 0;
    for (size_t i = 0; i < vSamples.size(); ++i)
    {
      AA = atof(tokens[5+3*i].c_str());
      AB = atof(tokens[5+3*i+1].c_str());
      BB = atof(tokens[5+3*i+2].c_str());
      if (AA == 0 && AB == 0 && BB == 0)
	mSnps[tokens[1]].vIsNa[i] = true;
      else
      {
	mSnps[tokens[1]].vGenos[i] = 0 * AA + 1 * AB + 2 * BB;
	maf += mSnps[tokens[1]].vGenos[i];
      }
    }
    maf /= 2 * (vSamples.size()
		- count (mSnps[tokens[1]].vIsNa.begin(),
			 mSnps[tokens[1]].vIsNa.end(),
			 true));
    mSnps[tokens[1]].maf = maf <= 0.5 ? maf : (1 - maf);
  }
  
  genoStream.close();
  
  for (map<string, Snp>::iterator it = mSnps.begin();
       it != mSnps.end(); ++it)
    if (it->second.chr.empty())
    {
      cerr << "ERROR: SNP " << it->first << " has no genotype" << endl;
      exit (1);
    }
}

void
getAbfsFromTruthFile (
  map<string, FtrMultiS> & mFeatures,
  const string & truthFile,
  const string & whichBf,
  const map<string, Snp> & mSnps,
  const int & verbose)
{  
  if (verbose > 0)
    cout << "load true ABFs ..." << endl << flush;
  
  ifstream truthStream;
  openFile (truthFile, truthStream);
  
  // check header line
  string line;
  getline (truthStream, line);
  if (line.empty())
  {
    cerr << "ERROR: file " << truthFile << " is empty" << endl;
    exit (1);
  }
  
  vector<string> tokens;
  split (line, " \t", tokens);
  if ((whichBf.compare("abf.meta") == 0
       || whichBf.compare("abf.fix") == 0)
      && tokens.size() < 7)
  {
    cerr << "ERROR: file " << truthFile << " should have at least"
	 << " 7 columns in his header" << endl;
    exit (1);
  }
  if ((whichBf.compare("abf.meta.avg.subset") == 0
       || whichBf.compare("abf.meta.avg.all") == 0)
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
  if ((whichBf.compare("abf.meta.avg.subset") == 0
       || whichBf.compare("abf.meta.avg.all") == 0)
      && (tokens[8].compare("l10abf.meta.avg.subset.c") != 0
	  || tokens[9].compare("l10abf.meta.avg.all.c") != 0))
  {
    cerr << "ERROR: header line of file " << truthFile
	 << " should start like this:" << endl
	 << "ftr snp coord nb.subgroups nb.samples l10abf.meta l10abf.fix l10abf.maxh l10abf.meta.avg.subset.c l10abf.meta.avg.all.c" << endl;
    exit (1);
  }
  
  while (true)
  {
    getline (truthStream, line);
    if (line.empty())
      break;
    split (line, " \t", tokens);
    if (mFeatures.find(tokens[0]) == mFeatures.end())
      continue;
    if (mSnps.find(tokens[1]) == mSnps.end())
      continue;
    if (whichBf.compare("abf.meta") == 0
	&& atof(tokens[5].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[5].c_str());
    else if (whichBf.compare("abf.fix") == 0
	     && atof(tokens[6].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[6].c_str());
    else if (whichBf.compare("abf.meta.avg.subset") == 0
	     && atof(tokens[8].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[8].c_str());
    else if (whichBf.compare("abf.meta.avg.all") == 0
	     && atof(tokens[9].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[9].c_str());
  }
  
  truthStream.close();
  
  for (map<string, FtrMultiS>::iterator it = mFeatures.begin();
       it != mFeatures.end(); ++it)
    if (it->second.maxL10TrueAbf == -1)
    {
      cerr << "ERROR: feature " << it->first << " has no true ABF" << endl;
      exit (1);
    }
}

void
getSummaryStatsForAllSubgroups (
  vector< vector<double> > & vvAllSumStatsCorr,
  const FtrMultiS & iFtr,
  const Snp & iSnp,
  const size_t & nbSubgroups,
  const size_t & nbSamples,
  gsl_permutation * perm)
{
  vvAllSumStatsCorr.clear ();
  vector<double> g, y, vSumStats;
  double betaPval, R2;
  for (size_t s = 0; s < nbSubgroups; ++s)
  {
    g.clear();
    y.clear();
    vSumStats.assign (4, 0);
    for (size_t i = 0; i < nbSamples; ++i)
      if (! iFtr.vvIsNa[s][i] && ! iSnp.vIsNa[i])
      {
	g.push_back (iSnp.vGenos[i]);
	y.push_back (iFtr.vvPhenos[s][gsl_permutation_get (perm, i)]);
      }
    vSumStats[0] = g.size(); // n
    ols (iFtr.name, iSnp.name, g, y, &(vSumStats[1]), &(vSumStats[2]),
	 &(vSumStats[3]), &betaPval, &R2, 0);
    vvAllSumStatsCorr.push_back (vector<double> ());
    correctSummaryStats (vvAllSumStatsCorr[s], vSumStats, 0);
  }
}

void
getWeightedAbf (
  const vector< vector<double> > & grid,
  const string & whichBf,
  const vector< vector<double> > & vvAllSumStatsCorr,
  const vector<double> & vWeights,
  double & l10Abf)
{
  double l10Abfs[grid.size()];
  for (size_t gridIdx = 0; gridIdx < grid.size(); ++gridIdx)
  {
    if (whichBf.compare("abf.meta") == 0)
      l10Abfs[gridIdx] = getAbfFromStdSumStats (vvAllSumStatsCorr,
						grid[gridIdx][0],
						grid[gridIdx][1],
						0);
    else if (whichBf.compare("abf.fix") == 0)
      l10Abfs[gridIdx] = getAbfFromStdSumStats (vvAllSumStatsCorr,
						0,
						grid[gridIdx][0]
						+ grid[gridIdx][1],
						0);
  }
  l10Abf = log10_weighted_sum (l10Abfs, &(vWeights[0]), grid.size());
}

void
getAbfOverConfigs (
  const vector< vector<double> > & grid,
  const string & whichBf,
  const vector< vector<double> > & vvAllSumStatsCorr,
  gsl_combination * comb,
  map<string, double> & mAllConfigsAbfs,
  const vector<double> & vWeightsGrid,
  double & l10Abf)
{
  size_t nbSubgroups = vvAllSumStatsCorr.size();
  stringstream ssConfig;
  vector< vector<double> > vvSomeSumStatsCorr;
  double l10Abfs[grid.size()];
  
  mAllConfigsAbfs.clear();
  
  // config "consistent"
  for(size_t s = 0; s < nbSubgroups; ++s)
    ssConfig << (s+1);
  getWeightedAbf (grid, "abf.meta", vvAllSumStatsCorr, vWeightsGrid,
		  mAllConfigsAbfs[ssConfig.str()]);
  
  // configs "single-specific"
  ssConfig.clear();
  comb = gsl_combination_calloc (nbSubgroups, 1);
  if (comb == NULL)
  {
    cerr << "ERROR: can't allocate memory for the combination" << endl;
    exit (1);
  }
  while (true)
  {
    vvSomeSumStatsCorr.clear();
    for (size_t i = 0; i < comb->k; ++i)
    {
      ssConfig << gsl_combination_get (comb, i) + 1;
      vvSomeSumStatsCorr.push_back (vvAllSumStatsCorr[
				      gsl_combination_get (comb, i)]);
    }
    for (size_t gridIdx = 0; gridIdx < grid.size(); ++gridIdx)
      l10Abfs[gridIdx] = getAbfFromStdSumStats (vvSomeSumStatsCorr,
						grid[gridIdx][0],
						grid[gridIdx][1],
						0);
    mAllConfigsAbfs.insert (make_pair(ssConfig.str(),
				      log10_weighted_sum (l10Abfs,
							  &(vWeightsGrid[0]),
							  grid.size())));
    if (gsl_combination_next (comb) != GSL_SUCCESS)
      break;
  }
  gsl_combination_free (comb);
  
  // all other configs
  if (whichBf.compare("abf.meta.avg.all") == 0)
  {
    for (size_t s = 2; s < nbSubgroups; ++s)
    {
      ssConfig.clear();
      comb = gsl_combination_calloc (nbSubgroups, s);
      if (comb == NULL)
      {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit (1);
      }
      while (true)
      {
	vvSomeSumStatsCorr.clear();
	for (size_t i = 0; i < comb->k; ++i)
	{
	  ssConfig << gsl_combination_get (comb, i) + 1;
	  vvSomeSumStatsCorr.push_back (vvAllSumStatsCorr[
					  gsl_combination_get (comb, i)]);
	}
	for (size_t gridIdx = 0; gridIdx < grid.size(); ++gridIdx)
	  l10Abfs[gridIdx] = getAbfFromStdSumStats (vvSomeSumStatsCorr,
						    grid[gridIdx][0],
						    grid[gridIdx][1],
						    0);
	mAllConfigsAbfs.insert (make_pair(ssConfig.str(),
					  log10_weighted_sum (l10Abfs,
							      &(vWeightsGrid[0]),
							      grid.size())));
	if (gsl_combination_next (comb) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free (comb);
    }
  }

  // average over all configs
  vector<double> vL10Abfs,
    vWeightsConfig (mAllConfigsAbfs.size(),
		    1.0 / ((double) mAllConfigsAbfs.size()));
  for (map<string, double>::iterator it = mAllConfigsAbfs.begin();
       it != mAllConfigsAbfs.end(); ++it)
    vL10Abfs.push_back (it->second);
  l10Abf = log10_weighted_sum (&vL10Abfs[0], &(vWeightsConfig[0]),
			       vL10Abfs.size());
}

void
computeJointAnalysisPermPvaluesFtrPerFtr (
  map<string, FtrMultiS> & mFeatures,
  map<string, Snp> & mSnps,
  const vector< vector<double> > & grid,
  const string & whichBf,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const size_t & nbSubgroups,
  const size_t & nbSamples,
  const int & verbose)
{
  if (verbose > 0)
    cout << "compute feature-level permutation P-values"
	 << " for the joint analysis ..." << endl
	 << "testStat=" << whichBf
         << " seed=" << seed
	 << " trick=" << trick
	 << endl << flush;
  
  FtrMultiS iFtr;
  Snp iSnp;
  size_t permId;
  vector<double> y, g, vSumStats, vWeights (grid.size(),
					    1.0 / (double) grid.size());
  double l10Abf, maxL10Abf;
  vector< vector<double> > vvAllSumStatsCorr;
  vector<size_t> vCounters = getCounters (mFeatures.size(), 5);
  
  gsl_rng_env_setup();
  gsl_rng * rngPerm = gsl_rng_alloc (gsl_rng_default);
  if (rngPerm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit (1);
  }
  gsl_rng_set (rngPerm, seed);
  gsl_rng * rngTrick = NULL;
  if (trick != 0)
  {
    rngTrick = gsl_rng_alloc (gsl_rng_default);
    if (rngTrick == NULL)
    {
      cerr << "ERROR: can't allocate memory for the RNG" << endl;
      exit (1);
    }
    gsl_rng_set (rngTrick, seed);
  }
  
  gsl_permutation * perm = gsl_permutation_calloc (nbSamples);
  if (perm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the permutation" << endl;
    exit (1);
  }
  
  gsl_combination * comb = NULL;
  map<string, double> mAllConfigsAbfs;
  
  size_t countFtrs = 0;
  for (map<string, FtrMultiS>::iterator itF = mFeatures.begin();
       itF != mFeatures.end(); ++itF)
  {
    iFtr = itF->second;
    ++countFtrs;
    if (verbose > 0)
      printCounter (countFtrs, vCounters);
    if (verbose > 1)
      cout << setfill('0') << setw((int)floor(log10(mFeatures.size()))+1)
	   << countFtrs << "/" << mFeatures.size() << " " << iFtr.name
	   << " (" << iFtr.vCisSnps.size() << " cis SNPs)" << endl << flush;
    
    bool shuffleOnly = false;
    for (permId = 0; permId < nbPerms; ++permId)
    {
      gsl_ran_shuffle (rngPerm, perm->data, perm->size, sizeof(size_t));
      if (shuffleOnly)
	continue;
      ++itF->second.nbPermsSoFar;
      maxL10Abf = 0;
      if (verbose > 1 && (permId+1) % 500 == 0)
	cout << setfill('0') << setw((int)floor(log10(nbPerms))+1)
	     << (permId+1) << "/" << nbPerms << "\r" << flush;
      
      for (vector<Snp *>::iterator itS = iFtr.vCisSnps.begin();
	   itS != iFtr.vCisSnps.end(); ++itS)
      {
	iSnp = *(*itS);
	
	getSummaryStatsForAllSubgroups (vvAllSumStatsCorr, iFtr, iSnp,
					nbSubgroups, nbSamples, perm);
	
	if (whichBf.find("avg") == string::npos)
	  getWeightedAbf (grid, whichBf, vvAllSumStatsCorr, vWeights, l10Abf);
	else
	  getAbfOverConfigs (grid, whichBf, vvAllSumStatsCorr, comb,
			     mAllConfigsAbfs, vWeights, l10Abf);
	if (l10Abf > maxL10Abf)
	  maxL10Abf = l10Abf;
      }
      if (maxL10Abf >= iFtr.maxL10TrueAbf)
	++iFtr.permPval;
      if (trick != 0 && iFtr.permPval == 11)
      {
	if (trick == 1)
	  break;
	else if (trick == 2)
	  shuffleOnly = true;
      }
    }
    
    if (itF->second.nbPermsSoFar == nbPerms)
      itF->second.permPval /= ((double) (nbPerms + 1));
    else
      itF->second.permPval = gsl_ran_flat (rngTrick,
					   (11 / ((double) (itF->second.nbPermsSoFar + 2))),
					   (11 / ((double) (itF->second.nbPermsSoFar + 1))));
  }
  
  gsl_permutation_free (perm);
  gsl_rng_free (rngPerm);
  if (trick != 0)
    gsl_rng_free (rngTrick);
}

void
writeResults (
  const map<string, FtrMultiS> & mFeatures,
  const string & outFile,
  const int & trick,
  const int & verbose)
{
  if (verbose > 0)
    cout << "write results in " << outFile << " ..." << endl << flush;
  
  ofstream outStream;
  openFile (outFile, outStream);
  outStream << "ftr jointPermPval";
  if (trick != 0)
    outStream << " nbPermsSoFar";
  outStream << endl << flush;
  
  for (map<string, FtrMultiS>::const_iterator it = mFeatures.begin();
       it != mFeatures.end(); ++it)
  {
    if (! outStream.good())
    {
      cerr << "ERROR: problem while writing output file" << endl;
      exit (1);
    }
    outStream << it->first << " " << it->second.permPval;
    if (trick != 0)
      outStream << " " << it->second.nbPermsSoFar;
    outStream << endl << flush;
  }
  
  outStream.close();
}

void getAndCheckSamples (
  const string & genoFile,
  const vector<string> & vPhenoFiles,
  const vector<string> & vSamplesToSkip,
  const vector<size_t> & vPermIndices,
  size_t & nbSamples,
  vector<bool> & vIdxSamplesToSkip,
  gsl_permutation * perm,
  const int & verbose)
{
  string line;
  size_t totNbSamples = 0;
  vector<string> tokens;
  
  // load all samples from phenotype file
  ifstream phenoStream;
  openFile (phenoFile, phenoStream);
  getline (phenoStream, line);
  phenoStream.close();
  split (line, " \t", tokens);
  if (tokens[0].compare("Id") == 0)
    tokens.erase (tokens.begin());
  totNbSamples = tokens.size();
  vIdxSamplesToSkip.assign (totNbSamples, false);
  vector<string> vSamples;
  for(size_t i = 0; i < totNbSamples; ++i)
  {
    if(! vSamplesToSkip.empty() &&
       find(vSamplesToSkip.begin(), vSamplesToSkip.end(), tokens[i])
       != vSamplesToSkip.end())
      vIdxSamplesToSkip[i] = true;
    vSamples.push_back (tokens[i]);
  }
  
  // check they are identical and in same order in genotype file (IMPUTE format)
  ifstream genoStream;
  openFile (genoFile, genoStream);
  getline (genoStream, line);
  genoStream.close();
  split (line, " \t", tokens);
  if (tokens.size() - 5 != totNbSamples)
  {
    cerr << "ERROR: " << tokens.size() - 5 << "samples in genotype file"
	 << endl << "compare to " << totNbSamples << " in phenotype file"
	 << endl;
    exit (1);
  }
  for(size_t i = 0; i < totNbSamples; ++i)
  {
    if (tokens[i+5].compare(vSamples[i]) != 0)
    {
      cerr << "ERROR: samples are not in same order between genotype"
	   << " and phenotype files" << endl;
      exit (1);
    }
  }
  
  nbSamples = count (vIdxSamplesToSkip.begin(), vIdxSamplesToSkip.end(), false);
  if (verbose > 0)
    cout << "nb of samples: " << nbSamples << endl << flush;
  
  // initialize the permutation struct if required
  if (! vPermIndices.empty())
  {
    if (vPermIndices.size() != totNbSamples)
    {
      cerr << "ERROR: wrong number of permutation indices ("
	   << vPermIndices.size() << " versus " << totNbSamples
	   << ")" << endl;
      exit (1);
    }
    perm = gsl_permutation_calloc (vPermIndices.size());
    if (perm == NULL)
    {
      cerr << "ERROR: can't allocate memory for the permutation" << endl;
      exit (1);
    }
    for (size_t i = 0; i < vPermIndices.size(); ++i)
      perm->data[i] = vPermIndices[i];
  }
}

void writeResults (
  const vector<StatsFtr> & vStatsFtrs,
  const string & outFile,
  const bool & calcSpearman,
  const size_t & nbPerms,
  const int & verbose)
{
  ofstream outStream;
  openFile (outFile, outStream);
  outStream << "ftr jointPermPval" << endl;
  
  for (vector<StatsFtr>::const_iterator it = vStatsFtrs.begin();
       it != vStatsFtrs.end(); ++it)
    StatsFtr_write (*it, outStream);
  
  if (verbose > 0)
    cout << "results written in " << outFile << endl;
}
*/
void
run (
  const string & genoPathsFile,
  const string & phenoPathsFile,
  const string & ftrCoordsFile,
  const string & gridFile,
  const string & whichBf,
  const string & outPrefix,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const string & ftrsToKeepFile,
  const string & snpsToKeepFile,
  const bool & needQnorm,
  const int & verbose)
{
  vector<string> vGenoPaths = loadOneColumnFile (genoPathsFile, verbose);
  vector<string> vPhenoPaths = loadOneColumnFile (phenoPathsFile, verbose);
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsToKeepFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsToKeepFile, verbose);
  vector<vector<double> > grid = loadGrid (gridFile, verbose);
  
  vector<string> vSamples;
  vector<vector<size_t> > vvSampleIdxs;
  loadSamples (vGenoPaths, vPhenoPaths, vSamples, vvSampleIdxs, verbose);
/*  
  map<string, sFtr> mFtrs;
  loadPhenos (vPhenoPaths, vFtrsToKeep, mFtrs, verbose);
  loadFtrInfo (ftrCoordsFile, mFtrs, verbose);
  
  map<string, sSnp> mSnps;
  loadSnpInfoAndGenos (vGenoPaths, mSnps, verbose);
  
  testAssos (mFtrs, mSnps, verbose, grid, whichBf, nbPerms, seed, trick,
	     needQnorm);
  
	     writeRes (outPrefix, verbose);
*/
}

int main (int argc, char ** argv)
{
  int verbose = 1, trick = 0;
  string genoPathsFile, phenoPathsFile, ftrCoordsFile, gridFile,
    whichBf = "abf.const", outPrefix, ftrsToKeepFile, snpsToKeepFile;
  size_t nbPerms = 10000, seed = string::npos;
  bool needQnorm = false;
  
  parseArgs (argc, argv, genoPathsFile, phenoPathsFile, ftrCoordsFile,
	     gridFile, whichBf, outPrefix, nbPerms, seed,
	     trick, ftrsToKeepFile, snpsToKeepFile, needQnorm, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl;
  }
  
  run (genoPathsFile, phenoPathsFile, ftrCoordsFile,
       gridFile, whichBf, outPrefix, nbPerms, seed,
       trick, ftrsToKeepFile, snpsToKeepFile, needQnorm, verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime)
	 << ": elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << ")" << endl;
  }
  
  return EXIT_SUCCESS;
}
