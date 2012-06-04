/** \file get_joint-analysis_ftr-level_perm-pval.cpp
 *
 *  `get_joint-analysis_ftr-level_perm-pval' computes feature-level permutation P-values.
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
 *  g++ -Wall -O3 -fopenmp utils.cpp get_joint-analysis_ftr-level_perm-pval.cpp -lgsl -lgslcblas -o get_joint-analysis_ftr-level_perm-pval
 *  help2man -o get_joint-analysis_ftr-level_perm-pval.man ./get_joint-analysis_ftr-level_perm-pval
 *  groff -mandoc get_joint-analysis_ftr-level_perm-pval.man > get_joint-analysis_ftr-level_perm-pval.ps
*/

#include <cmath>
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
#include "get_summary_stats.cpp"
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
       << "  -g, --geno\tfile with genotypes in IMPUTE format" << endl
       << "\t\ta header line with sample names is required" << endl
       << "\t\tsamples in columns should be in same order as in phenotype file" << endl
       << "  -p, --pheno\trelative path to directory with one phenotypes file per subgroup" << endl
       << "\t\trow 1 for sample names, column 1 for feature names" << endl
       << "\t\teach subgroup should have the same set of features" << endl
       << "\t\tcan also be a single file (ie. ABFs for a single subgroup)" << endl
       << "      --fcoord\tBED file with the features coordinates" << endl
       << "\t\tfeatures should be in same order than in phenotypes files" << endl
       << "\t\tfield delimiters are tabs" << endl
       << "      --grid\tfile with the grid of values for phi2 and omega2 (ES model)" << endl
       << "\t\tsee GetGridPhiOmega() in package Rquantgen" << endl
       << "      --truth\tfile with the ABFs on the original dataset" << endl
       << "\t\toutput from get_abf_meta" << endl
       << "      --abf\twhich ABF to use as the test statistic for the permutations" << endl
       << "\t\tdefault=abf.meta/abf.fix/abf.meta.avg.all/abf.meta.avg.subset" << endl
       << "  -o, --out\tfile that will contain the permutation P-values" << endl
       << "  -P, --perm\tnumber of phenotype permutations at each feature" << endl
       << "\t\tdefault=10000" << endl
       << "      --seed\tseed for the two random number generators" << endl
       << "\t\tone RNG is used for the permutations and another for the trick" << endl
       << "\t\tby default, both are initialized via microseconds from epoch" << endl
       << "      --trick\tapply trick to speed-up permutations" << endl
       << "\t\tstop after the tenth permutation for which the test statistic" << endl
       << "\t\tis better than or equal to the true value, and sample from" << endl
       << "\t\ta uniform between 11/(nbPermsSoFar+2) and 11/(nbPermsSoFar+1)" << endl
       << "\t\tif '1', the permutations really stops" << endl
       << "\t\tif '2', all permutations are done but the test statistics are not computed" << endl
       << "\t\tallow to compare different test statistics on the same permutations" << endl
       << "  -f, --ftr\tfile with a list of features to analyze" << endl
       << "\t\tone feature name per line" << endl
       << "  -l, --links\tfile with links between genes and SNPs" << endl
       << "\t\tcustom format: feature<space/tab>SNP|coord" << endl
       << "\t\tfeatures should be in same order than in phenotype file" << endl
       << "\t\tuseful to focus on genetic variants in cis (use windowBed)" << endl
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
  string & gridFile,
  string & truthFile,
  string & whichAbf,
  string & outFile,
  size_t & nbPermutations,
  size_t & seed,
  int & trick,
  string & ftrsToKeepFile,
  string & linksFile,
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
	{"truth", required_argument, 0, 0},
	{"abf", required_argument, 0, 0},
	{"out", required_argument, 0, 'o'},
	{"perm", required_argument, 0, 'P'},
	{"seed", required_argument, 0, 0},
	{"trick", required_argument, 0, 0},
	{"ftr", required_argument, 0, 'f'},
	{"links", required_argument, 0, 'l'},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:o:P:f:l:",
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
      if (strcmp(long_options[option_index].name, "seed") == 0)
      {
	seed = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "trick") == 0)
      {
	trick = atoi (optarg);
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
      genoFile = optarg;
      break;
    case 'p':
      phenoDir = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'P':
      nbPermutations = atol (optarg);
      break;
    case 'f':
      ftrsToKeepFile = optarg;
      break;
    case 'l':
      linksFile = optarg;
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
  ofstream outStream;
  openFile (outFile, outStream);
  outStream.close();
  if (whichAbf.compare("abf.meta") == 0 && whichAbf.compare("abf.fix") == 0
      && whichAbf.compare("abf.meta.avg.subset") == 0
      && whichAbf.compare("abf.meta.avg.all") == 0)
  {
    fprintf (stderr, "ERROR: --abf should be abf.meta, abf.fix, abf.meta.avg.subset or abf.meta.avg.all.\n\n");
    help (argv);
    exit (1);
  }
  if (trick != 0 && trick != 1 && trick != 2)
  {
    fprintf (stderr, "ERROR: unrecognized trick %d (--trick).\n\n", trick);
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
  if (seed == string::npos)
    seed = getSeed ();
}

struct Snp
{
  string name; // eg. rs2205177
  string chr; // eg. chr21
  size_t coord; // 1-based coordinate
  vector<double> vGenos; // genotypes of samples
  vector<bool> vIsNa; // missing values
  double maf; // minor allele frequency
};

struct Feature
{
  string name; // eg. ENSG00000182816
  string chr; // eg. chr21
  size_t start; // 1-based coordinate
  size_t end; // idem
  int bin; // UCSC binning system
  vector< vector<double> > vvPhenos; // phenotypes of samples per subgroup
  vector< vector<bool> > vvIsNa; // missing values per subgroup
  vector <Snp *> vCisSnps;
  double maxL10TrueAbf; // max log10(ABF) over SNPs on non-permuted data
  double permPval; // permutation P-value of the joint analysis
  size_t nbPermsSoFar;
};

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
Feature_init (
  Feature & iFtr,
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

void loadFeaturesFromPhenoFiles(
  map<string, Feature> & mFeatures,
  vector<string> & vSamples,
  const vector<string> & vPhenoFiles,
  const vector<string> & vFtrsToKeep,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load phenotype values from each subgroup ..." << endl << flush;
  
  ifstream phenoStream;
  string line;
  vector<string> tokens;
  
  for (size_t s = 0; s < vPhenoFiles.size(); ++s)
  {
    openFile (vPhenoFiles[s], phenoStream);
    
    // check header line
    getline (phenoStream, line);
    if (line.empty())
    {
      cerr << "ERROR: file " << vPhenoFiles[s] << " is empty" << endl;
      exit (1);
    }
    if (vSamples.empty())
      split (line, " \t", vSamples);
    else
    {
      split (line, " \t", tokens);
      if (tokens.size() != vSamples.size())
      {
	cerr << "ERROR: all phenotype files don't have same number of samples"
	     << endl;
	exit (1);
      }
      for (size_t i = 0; i < vSamples.size(); ++i)
	if (tokens[i].compare(vSamples[i]) != 0)
	{
	  cerr << "ERROR: all phenotype files don't have the same samples"
	       << endl << "or they are not in the same order" << endl;
	  exit (1);
	}
    }
    if (verbose > 0 && s == 0)
      cout << "nb of samples: " << vSamples.size() << endl << flush;
    
    // load phenotype values
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
	Feature iFtr;
	Feature_init (iFtr, tokens[0], vPhenoFiles.size(), vSamples.size());
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
loadFtrCoordsFromBedFile (
  map<string, Feature> & mFeatures,
  const string & ftrCoordsFile,
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
loadLinksFtrCisSnpsFromFile (
  map<string, Feature> & mFeatures,
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
  
  map<string, Feature>::iterator it = mFeatures.begin();
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
  map<string, Feature> & mFeatures,
  const string & truthFile,
  const string & whichAbf,
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
  if ((whichAbf.compare("abf.meta") == 0
       || whichAbf.compare("abf.fix") == 0)
      && tokens.size() < 7)
  {
    cerr << "ERROR: file " << truthFile << " should have at least"
	 << " 7 columns in his header" << endl;
    exit (1);
  }
  if ((whichAbf.compare("abf.meta.avg.subset") == 0
       || whichAbf.compare("abf.meta.avg.all") == 0)
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
  if ((whichAbf.compare("abf.meta.avg.subset") == 0
       || whichAbf.compare("abf.meta.avg.all") == 0)
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
    if (whichAbf.compare("abf.meta") == 0
	&& atof(tokens[5].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[5].c_str());
    else if (whichAbf.compare("abf.fix") == 0
	     && atof(tokens[6].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[6].c_str());
    else if (whichAbf.compare("abf.meta.avg.subset") == 0
	     && atof(tokens[8].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[8].c_str());
    else if (whichAbf.compare("abf.meta.avg.all") == 0
	     && atof(tokens[9].c_str()) > mFeatures[tokens[0]].maxL10TrueAbf)
      mFeatures[tokens[0]].maxL10TrueAbf = atof(tokens[9].c_str());
  }
  
  truthStream.close();
  
  for (map<string, Feature>::iterator it = mFeatures.begin();
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
  const Feature & iFtr,
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
  const string & whichAbf,
  const vector< vector<double> > & vvAllSumStatsCorr,
  const vector<double> & vWeights,
  double & l10Abf)
{
  double l10Abfs[grid.size()];
  for (size_t gridIdx = 0; gridIdx < grid.size(); ++gridIdx)
  {
    if (whichAbf.compare("abf.meta") == 0)
      l10Abfs[gridIdx] = getAbfFromStdSumStats (vvAllSumStatsCorr,
						grid[gridIdx][0],
						grid[gridIdx][1],
						0);
    else if (whichAbf.compare("abf.fix") == 0)
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
  const string & whichAbf,
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
  if (whichAbf.compare("abf.meta.avg.all") == 0)
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
  map<string, Feature> & mFeatures,
  map<string, Snp> & mSnps,
  const vector< vector<double> > & grid,
  const string & whichAbf,
  const size_t & nbPermutations,
  const size_t & seed,
  const int & trick,
  const size_t & nbSubgroups,
  const size_t & nbSamples,
  const int & verbose)
{
  if (verbose > 0)
    cout << "compute feature-level permutation P-values"
	 << " for the joint analysis ..." << endl
	 << "testStat=" << whichAbf
         << " seed=" << seed
	 << " trick=" << trick
	 << endl << flush;
  
  Feature iFtr;
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
  for (map<string, Feature>::iterator itF = mFeatures.begin();
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
    for (permId = 0; permId < nbPermutations; ++permId)
    {
      gsl_ran_shuffle (rngPerm, perm->data, perm->size, sizeof(size_t));
      if (shuffleOnly)
	continue;
      ++itF->second.nbPermsSoFar;
      maxL10Abf = 0;
      if (verbose > 1 && (permId+1) % 500 == 0)
	cout << setfill('0') << setw((int)floor(log10(nbPermutations))+1)
	     << (permId+1) << "/" << nbPermutations << "\r" << flush;
      
      for (vector<Snp *>::iterator itS = iFtr.vCisSnps.begin();
	   itS != iFtr.vCisSnps.end(); ++itS)
      {
	iSnp = *(*itS);
	
	getSummaryStatsForAllSubgroups (vvAllSumStatsCorr, iFtr, iSnp,
					nbSubgroups, nbSamples, perm);
	
	if (whichAbf.find("avg") == string::npos)
	  getWeightedAbf (grid, whichAbf, vvAllSumStatsCorr, vWeights, l10Abf);
	else
	  getAbfOverConfigs (grid, whichAbf, vvAllSumStatsCorr, comb,
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
    
    if (itF->second.nbPermsSoFar == nbPermutations)
      itF->second.permPval /= ((double) (nbPermutations + 1));
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
  const map<string, Feature> & mFeatures,
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
  
  for (map<string, Feature>::const_iterator it = mFeatures.begin();
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

int main (int argc, char ** argv)
{
  string genoFile, phenoDir, ftrCoordsFile, gridFile, truthFile,
    whichAbf = "abf.meta", outFile, ftrsToKeepFile, linksFile;
  size_t nbPermutations = 10000, seed = string::npos;
  int verbose = 1, trick = 0;
  
  parse_args (argc, argv, genoFile, phenoDir, ftrCoordsFile, gridFile,
	      truthFile, whichAbf, outFile, nbPermutations, seed,
	      trick, ftrsToKeepFile, linksFile, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl;
  }
  
  vector<string> vPhenoFiles;
  if (isDirectory(phenoDir.c_str()))
    vPhenoFiles = scanInputDirectory (phenoDir, verbose);
  else
    vPhenoFiles.push_back (phenoDir);
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsToKeepFile, verbose);
  map<string, Feature> mFeatures;
  vector<string> vSamples;
  loadFeaturesFromPhenoFiles (mFeatures, vSamples, vPhenoFiles, vFtrsToKeep,
			      verbose);
  loadFtrCoordsFromBedFile (mFeatures, ftrCoordsFile, verbose);
  map<string, Snp> mSnps;
  loadLinksFtrCisSnpsFromFile (mFeatures, mSnps, linksFile, verbose);
  loadGenoSnpsFromFile (mSnps, genoFile, vSamples, verbose);
  getAbfsFromTruthFile (mFeatures, truthFile, whichAbf, mSnps, verbose);
  vector< vector<double> > grid = loadGrid (gridFile, verbose);
  computeJointAnalysisPermPvaluesFtrPerFtr (mFeatures,
  					    mSnps,
  					    grid,
  					    whichAbf,
  					    nbPermutations,
					    seed,
					    trick,
					    vPhenoFiles.size(),
					    vSamples.size(),
  					    verbose);
  writeResults (mFeatures, outFile, trick, verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime)
	 << ": elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << ")" << endl;
  }
  
  return EXIT_SUCCESS;
}
