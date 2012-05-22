/** \file get_abf_meta.cpp
  *
  *  `get_abf_meta' performs Bayesian meta-analysis of eQTL studies.
  *  Copyright (C) 2011  T. Flutre
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
  *  g++ -Wall -O3 get_abf_meta.cpp -lgsl -lgslcblas -o get_abf_meta
  *  help2man -o get_abf_meta.man ./get_abf_meta
  *  groff -mandoc get_abf_meta.man > get_abf_meta.ps
  */

#include <cmath>
#include <ctime>
#include <getopt.h>

#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <map>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm>
#include <utility>
#include <functional>
using namespace std;

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_combination.h>

#include "utils.cpp"

/** \brief Display the help on stdout.
 */
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       <<" performs the meta-analysis by computing three Bayes Factors," << endl
       << "one per model (meta-analysis, fixed-effects, max heterogeneity)," << endl
       << "for each pair feature-SNP." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS]..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (default=1)" << endl
       << "  -i, --input\tdirectory with the summary stats for each subgroup" << endl
       << "\t\tone file per subgroup, with at least 8 columns:" << endl
       << "\t\tftr chrF start end snp chrS coord maf n betahat sebetahat sigmahat" << endl
       << "\t\tsee output from other program 'get_summary_stats'" << endl
       << "\t\tcan also be a single file (ie. ABFs for a single subgroup)" << endl
       << "  -g, --grid\tfile with the grid of values for phi2 and omega2 (ES model)" << endl
       << "\t\tsee GetGridPhiOmega() in package Rquantgen" << endl
       << "  -o, --out\toutput file with the Bayes factors for each pair feature-SNP" << endl
       << "\t\t3 BFs are reported: meta-analysis, fixed-effects, max-heterogeneity" << endl
       << "  -f, --ftr\tfile with a list of features to analyze" << endl
       << "  -s, --snp\tfile with a list of SNP identifiers to analyze" << endl
       << "  -x, --index\twhich SNP info to use to index pairs feature-SNP" << endl
       << "\t\tdefault=id+coord, but can be only 'id' or 'coord'" << endl
       << "  -b, --beta\tsave the summary statistics of each subgroup in the output file" << endl
       << "  -c, --config\tcalculate the BF 'meta' for each configuration" << endl
       << "\t\tthe total number of configurations is 2^nb.subgroups" << endl
       << "\t\t2 average BFs are reported (over all configs, or only subgroup-" << endl
       << "\t\tconsistent and each single-subgroup-specific)" << endl
       << endl
       << "Examples:" << endl
       << "  " << argv[0] << " -i <input> -g <grid> -o <output>" << endl
       << endl
       << "References:" << endl
       << "  Wen and Stephens (http://arxiv.org/abs/1111.1210)" << endl;
}

/** \brief Display version and license information on stdout.
 */
void version (char ** argv)
{
  cout << argv[0] << " 0.1" << endl
       << endl
       << "Copyright (C) 2011 T. Flutre." << endl
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
  string & inDir,
  string & gridFile,
  string & outFile,
  string & ftrsFile,
  string & snpsFile,
  string & snpIdx,
  bool & saveSumStats,
  bool & calcAllConfigs,
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
	{"input", required_argument, 0, 'i'},
	{"grid", required_argument, 0, 'g'},
	{"output", required_argument, 0, 'o'},
	{"ftr", required_argument, 0, 'f'},
	{"snp", required_argument, 0, 's'},
	{"index", required_argument, 0, 'x'},
	{"beta", no_argument, 0, 'b'},
	{"config",no_argument, 0, 'c'},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:i:g:o:f:s:x:bc",
		     long_options, &option_index);
    if (c == -1)
      break;
    switch (c)
    {
    case 0:
      if (long_options[option_index].flag != 0)
	break;
      printf ("option %s", long_options[option_index].name);
      if (optarg)
	printf (" with arg %s", optarg);
      printf ("\n");
      break;
    case 'h':
      help (argv);
      exit (0);
    case 'V':
      version (argv);
      exit (0);
    case 'v':
      verbose = atoi(optarg);
      break;
    case 'i':
      inDir = optarg;
      break;
    case 'g':
      gridFile = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'f':
      ftrsFile = optarg;
      break;
    case 's':
      snpsFile = optarg;
      break;
    case 'x':
      snpIdx = optarg;
      break;
    case 'b':
      saveSumStats = true;
      break;
    case 'c':
      calcAllConfigs = true;
      break;
    case '?':
      break;
    default:
      abort ();
    }
  }
  if (inDir.empty())
  {
    fprintf (stderr, "ERROR: missing input directory with summary statistics (-i).\n\n");
    help (argv);
    exit (1);
  }
  if (gridFile.empty())
  {
    fprintf (stderr, "ERROR: missing grid file (-g).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (gridFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", gridFile.c_str());
    help (argv);
    exit (1);
  }
  if (outFile.empty())
  {
    fprintf (stderr, "ERROR: missing output file (-o).\n\n");
    help (argv);
    exit (1);
  }
  if (snpIdx.empty())
    snpIdx = "id+coord";
}

/** \brief Index the files with a map to know which pairs feature-SNP are 
 *  in common among subgroups.
 */
void
indexInputFiles (
  const vector<string> & vInFiles,
  const vector<string> & vFtrsToKeep,
  const vector<string> & vSnpsToKeep,
  const string & snpIdx,
  map<string, vector<size_t> > & mPairs2Positions,
  const int & verbose)
{
  size_t subgroup = 0, line_id = 0, pos = 0;
  string line;
  vector<string> tokens;
  map<string, vector<size_t> >::iterator it_mP2P;
  
  if (verbose > 0)
  {
    cout << "index input files ..." << endl;
  }
  
  for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
  {
    ifstream inStream;
    inStream.open (vInFiles[subgroup].c_str());
    if (! inStream.good())
    {
      cerr << "ERROR: can't open file " << vInFiles[subgroup] << endl;
      exit (1);
    }
    if (verbose > 0)
    {
      cout << "s" << subgroup+1 << " " << vInFiles[subgroup] << " ...";
      fflush (stdout);
    }
    
    // check first line is proper header
    getline (inStream, line);
    if (line.empty())
    {
      cerr << "ERROR: file " << vInFiles[subgroup] << " is empty" << endl;
      exit (1);
    }
    split (line, ' ', tokens);
    if (tokens.size() < 12)
    {
      cerr << "ERROR: file " << vInFiles[subgroup]
	   << " has less than 12 columns" << endl;
      exit (1);
    }
    if (tokens[0].compare("ftr") != 0
	|| tokens[1].compare("chrF") != 0
	|| tokens[2].compare("start") != 0
	|| tokens[3].compare("end") != 0
	|| tokens[4].compare("snp") != 0
	|| tokens[5].compare("chrS") != 0
	|| tokens[6].compare("coord") != 0
	|| tokens[7].compare("maf") != 0
	|| tokens[8].compare("n") != 0
	|| tokens[9].compare("betahat") != 0
	|| tokens[10].compare("sebetahat") != 0
	|| tokens[11].compare("sigmahat") != 0)
    {
      cerr << "ERROR: header line of file " << vInFiles[subgroup]
	   << " should start like this:" << endl
	   << "ftr chrF start end snp chrS coord maf n betahat sebetahat sigmahat" << endl;
      exit (1);
    }
    
    line_id = 1;
    while (inStream.good())
    {
      pos = inStream.tellg();
      getline (inStream, line);
      if (line.empty())
      {
	break;
      }
      ++line_id;
      split (line, ' ', tokens);
      if (tokens.size() < 12)
      {
	cerr << endl << "ERROR: file " << vInFiles[subgroup] 
	     << " at line " << line_id 
	     << " should have at least 12 fields:" << endl
	     << "ftr chrF start end snp chrS coord maf n betahat sebetahat sigmahat"
	     << endl;
	exit (1);
      }
      if (! vFtrsToKeep.empty() &&
	  find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0]) ==
	  vFtrsToKeep.end())
      {
	continue;  // skip line, based on feature name
      }
      if (! vSnpsToKeep.empty() &&
	  find(vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[4]) ==
	  vSnpsToKeep.end())
      {
	continue;  // skip line, based on SNP identifier
      }
      stringstream ss; // ftr|snp-id or ftr|snp-coord or ftr|snp-id|snp-coord
      ss << tokens[0];
      if (snpIdx.compare("id") == 0)
	ss << "|" << tokens[4];
      else if (snpIdx.compare("coord") == 0)
	ss << "|" << tokens[6];
      else if (snpIdx.compare("id+coord") == 0)
	ss << "|" << tokens[4] << "|" << tokens[6];
      if (mPairs2Positions.find(ss.str()) == mPairs2Positions.end())
      {
	vector<size_t> vPositions (vInFiles.size(), string::npos);
	mPairs2Positions.insert ( make_pair (ss.str(), vPositions) );
      }
      mPairs2Positions[ss.str()][subgroup] = pos;
    }
    inStream.close ();
    if (verbose > 0)
    {
      cout << " (" << line_id << " lines)" << endl;
      fflush (stdout);
    }
  }
  
  if (verbose > 0)
  {
    cout << "nb of pairs feature-SNP: " << mPairs2Positions.size() << endl;
  }
}

/** \brief Load the file with the grid values.
 */
vector< vector<double> >
loadGrid (
  const string & gridFile,
  const int & verbose)
{
  if (verbose > 0)
  {
    cout << "load grid ..." << endl;
  }
  vector< vector<double> > grid;
  ifstream gridStream;
  vector<string> tokens;
  string line;
  gridStream.open (gridFile.c_str());
  if (! gridStream.is_open())
  {
    cerr << "ERROR: can't open file " << gridFile << endl;
    exit (1);
  }
  while (gridStream.good())
  {
    getline (gridStream, line);
    if (line.empty())
    {
      break;
    }
    split (line, ' ', tokens);
    if (tokens.size() != 2)
    {
      cerr << "ERROR: format of file " << gridFile
	   << " should be phi2<space>oma2" << endl;
      exit (1);
    }
    vector<double> grid_values;
    grid_values.push_back (atof (tokens[0].c_str()));
    grid_values.push_back (atof (tokens[1].c_str()));
    grid.push_back (grid_values);
  }
  gridStream.close();
  if (verbose > 0)
  {
    cout << "grid size: " << grid.size() << endl;
  }
  return grid;
}

void
writeOutputHeader (
  ofstream & outStream,
  const bool & calcAllConfigs,
  gsl_combination * c,
  const size_t & nbSubgroups,
  const bool & saveSumStats)
{
  size_t i, j;
  outStream << "ftr"
	    << " snp"
	    << " coord"
	    << " nb.subgroups"
	    << " nb.samples"
	    << " l10abf.meta"
	    << " l10abf.fix"
	    << " l10abf.maxh";
  
  if (calcAllConfigs)
  {
    outStream << " l10abf.meta.avg.all.c"
	      << " l10abf.meta.avg.subset.c";
    for(i=1; i<nbSubgroups; ++i)
    {
      c = gsl_combination_calloc (nbSubgroups, i);
      if (c == NULL)
      {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit (1);
      }
      while (true)
      {
	outStream << " l10abf.meta.c";
	for (j=0; j<c->k; ++j)
	  outStream << gsl_combination_get (c, j) + 1;
	if (gsl_combination_next (c) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free (c);
    }
  }
  
  if (saveSumStats)
    for (i = 0; i < nbSubgroups; ++i)
      outStream << " betahat.s" << i+1
		<< " sebetahat.s" << i+1
		<< " sigmahat.s" << i+1;
  
  outStream << endl;
}

/** \brief Retrieve the summary stats (n,betahat,sebetahat,sigmahat) for a 
 *  given pair feature-SNP in each  subgroup in which this pair is present.
 *
 *  \note If the pair is not present in the given subgroup (missing value),
 *  n will be 0 and the other summary stats will be NaN. If the pair is present
 *  but has infinite summary stats (should be coded as Inf in input file), only
 *  n will be an integer, the other summary stats remaining Inf.
 */
void
getSummaryStatsForPair (
  const string & pairId,
  const string & snpIdx,
  const vector<size_t> & vPositions,
  vector<ifstream *> & vPtInStreams,
  map<size_t, vector<double> > & mAllSumStats,
  size_t & nbSubgroupsUsed,
  size_t & nbSamplesUsed,
  const int & verbose)
{
  size_t subgroup = 0, colIdx = 0, nbSubgroups = vPositions.size();
  string line;
  vector<string> tokens;
  
  mAllSumStats.clear();
  nbSubgroupsUsed = 0;
  nbSamplesUsed = 0;
  
  for (subgroup = 0; subgroup < nbSubgroups; ++subgroup)
  {
    vector<double> vSumStats (4, numeric_limits<double>::quiet_NaN());
    
    if (vPositions[subgroup] == string::npos)
    {
      vSumStats[0] = 0;
      if (verbose > 0)
	cout << "subgroup " << subgroup+1
	     << ": skip because pair not present in this subgroup"
	     << endl;
    }
    
    else
    {
      ++nbSubgroupsUsed;
      vPtInStreams[subgroup]->seekg (vPositions[subgroup]);
      getline (*(vPtInStreams[subgroup]), line);
      if (line.find('\t') != string::npos)
	split (line, '\t', tokens);
      else
	split (line, ' ', tokens);
      stringstream ss;
      ss << tokens[0];
      if (snpIdx.compare("id") == 0)
	ss << "|" << tokens[4];
      else if (snpIdx.compare("coord") == 0)
	ss << "|" << tokens[6];
      else if (snpIdx.compare("id+coord") == 0)
	ss << "|" << tokens[4] << "|" << tokens[6];
      if (ss.str() != pairId)
      {
	cerr << "ERROR: stream position is wrong for pair " << pairId
	     << " in subgroup " << subgroup+1 << endl;
	cerr << line << endl;
	exit (1);
      }
      nbSamplesUsed += atoi (tokens[8].c_str());
      for (colIdx = 8; colIdx < 12; ++colIdx)
      {
	vSumStats[colIdx-8] = atof (tokens[colIdx].c_str());
      }
    }
    
    mAllSumStats.insert (make_pair(subgroup, vSumStats));
  }
}

/** \brief Apply the small sample size correction and return the standardized
 *  summary stats (bhat,sebhat,t).
 *
 *  \note If there is no variation in the phenotypes of the given subgroup,
 *  t is practically 0, and bhat and sebhat will be 0 too. If there is no 
 *  variation in the genotype (betahat=0, sebetahat=Inf), t will be 0, whereas
 *  bhat will be 0 and sebhat will be Inf.
 */
vector<double>
correctSummaryStats (
  const vector<double> & vSumStats,
  const int & verbose)
{
  vector<double> stdSumStatsCorr;
  double n = 0, bhat = 0, sebhat = 0, sigmahat = 0, t = 0;
  
  if (vSumStats[0] != 0)  // if=0 -> pair is not present in the given subgroup
  {
    n = vSumStats[0];
    bhat = vSumStats[1] / vSumStats[3];      // before correction
    sebhat = vSumStats[2] / vSumStats[3];    // before correction
#ifdef DEBUG
    if (verbose > 0)
      printf ("n=%i betahat=%e sebetahat=%e sigmahat=%e bhat=%e sebhat=%e\n",
	      (int) n, vSumStats[1], vSumStats[2], vSumStats[3], bhat, sebhat);
#endif
    t = gsl_cdf_gaussian_Pinv(gsl_cdf_tdist_P(-fabs(bhat/sebhat), n-2), 1.0);
#ifdef DEBUG
    if (verbose > 0)
      printf ("t=%.8f\n", t);
#endif
    if (fabs(t) > 1e-8)
    {
      sigmahat = fabs(vSumStats[1]) / (fabs(t) * sebhat);
      bhat = vSumStats[1] / sigmahat;
      sebhat = bhat / t;
    }
    else
    {
      sigmahat = numeric_limits<double>::quiet_NaN();
      bhat = 0;
      sebhat = numeric_limits<double>::infinity();
    }
#ifdef DEBUG
    if (verbose > 0)
      printf ("after: sigmahat=%e bhat=%e sebhat=%e t=%e\n",
	      sigmahat, bhat, sebhat, t);
#endif
    stdSumStatsCorr.push_back (bhat);
    stdSumStatsCorr.push_back (sebhat);
    stdSumStatsCorr.push_back (t);
  }
  
  return stdSumStatsCorr;
}

/** \brief Return the log10 approximate Bayes Factor given standardized 
 *  summary statistics and hyperparameters.
 *  \note If bhat=0 and sebhat=Inf for a given subgroup, its ABF will be 1.
 *  If this is the case for all subgroups, the overall ABF will also be 1.
 */
double
getAbfFromStdSumStats (
  const vector< vector<double> > & stdSumStats,
  const double & phi2,
  const double & oma2,
  const int & verbose)
{
  double lABF_all = 0;
  size_t subgroup = 0;
  double bhat = 0, varbhat = 0, t = 0, bbarhat_num = 0, bbarhat_denom = 0,
    varbbarhat = 0;
  vector<double> lABF_single_s;
#ifdef DEBUG
  if (verbose > 0)
    printf ("phi2=%f oma2=%f\n", phi2, oma2);
#endif
  
  for (subgroup = 0; subgroup < stdSumStats.size(); ++subgroup)
  {
    bhat = stdSumStats[subgroup][0];
    varbhat = pow (stdSumStats[subgroup][1], 2);
    t = stdSumStats[subgroup][2];
    double lABF_single;
    if (fabs(t) < 1e-8)
    {
      lABF_single = 0;
    }
    else
    {
      bbarhat_num += bhat / (varbhat + phi2);
      bbarhat_denom += 1 / (varbhat + phi2);
      varbbarhat += 1 / (varbhat + phi2);
      lABF_single = 0.5 * log10(varbhat)
	- 0.5 * log10(varbhat + phi2)
	+ (0.5 * pow(t,2) * phi2 / (varbhat + phi2)) / log(10);
  }
  lABF_single_s.push_back (lABF_single);
#ifdef DEBUG
    if (verbose > 0)
      printf ("lABF_single_s[%ld]=%e\n", subgroup+1, lABF_single_s[subgroup]);
#endif
  }
  
  double bbarhat = (bbarhat_denom != 0) ?
    bbarhat_num / bbarhat_denom
    : 0;
  varbbarhat = (varbbarhat != 0) ?
    1 / varbbarhat
    : numeric_limits<double>::infinity();
  double T2 = pow(bbarhat, 2.0) / varbbarhat;
  double lABF_bbar = (T2 != 0) ?
    0.5 * log10(varbbarhat) - 0.5 * log10(varbbarhat + oma2)
    + (0.5 * T2 * oma2 / (varbbarhat + oma2)) / log(10)
    : 0;
#ifdef DEBUG
  if (verbose > 0)
    printf ("bbarhat=%e varbbarhat=%e T2=%e lABF_bbar=%e\n",
	    bbarhat, varbbarhat, T2, lABF_bbar);
#endif
  
  lABF_all = lABF_bbar;
  for (subgroup = 0; subgroup < lABF_single_s.size(); ++subgroup)
    lABF_all += lABF_single_s[subgroup];
#ifdef DEBUG
  if (verbose > 0)
    printf ("lABF_all=%e\n", lABF_all);
#endif
  
  return lABF_all;
}

/** \brief Return log_{10}(\sum_i w_i 10^vec_i)
 */
double
log10_weighted_sum (
  double * vec,
  double * weights,
  size_t size)
{
  size_t i = 0;
  double max = vec[0];
  for (i = 0; i < size; i++)
    if (vec[i] > max)
      max = vec[i];
  double sum = 0;
  for (i = 0; i < size; i++)
    sum += weights[i] * pow(10, vec[i] - max);
  return max + log10(sum);
}

void
computeAbfMetaOfAllConfigsForPairOverGrid (
  map<string, double> & mAllConfigsAbfs,
  const size_t & nbSubgroupsUsed,
  gsl_combination * c,
  const vector< vector<double> > & mAllSumStatsCorr,
  const vector< vector<double> > & grid)
{
  size_t i, j, gridIdx;
  vector< vector<double> > mSomeSumStatsCorr;
  double l10_abf_metas_c[grid.size()], weights[grid.size()];
  
  mAllConfigsAbfs.clear();
  
  for(i=1; i<nbSubgroupsUsed; ++i)
  {
    c = gsl_combination_calloc (nbSubgroupsUsed, i);
    if (c == NULL)
    {
      cerr << "ERROR: can't allocate memory for the combination" << endl;
      exit (1);
    }
    
    while (true)
    {
      stringstream ss; // will contain the identifier of the configuration
      mSomeSumStatsCorr.clear();
      for (j=0; j<c->k; ++j)
      {
	ss << gsl_combination_get (c, j) + 1;
	mSomeSumStatsCorr.push_back (mAllSumStatsCorr[
				       gsl_combination_get (c, j)]);
      }
      
      for (gridIdx = 0; gridIdx < grid.size(); ++gridIdx)
      {
	l10_abf_metas_c[gridIdx] = getAbfFromStdSumStats (mSomeSumStatsCorr,
							  grid[gridIdx][0],
							  grid[gridIdx][1],
							  0);
	weights[gridIdx] = 1.0 / (double) grid.size();
      }
      
      mAllConfigsAbfs.insert (make_pair(ss.str(),
					log10_weighted_sum (l10_abf_metas_c,
							    weights,
							    grid.size())));
      if (gsl_combination_next (c) != GSL_SUCCESS)
	break;
    }
    
    gsl_combination_free (c);
  }
}

void
computeAveragesAbfMetaOverConfigs (
  map<string, double> & mAllConfigsAbfs,
  const double & l10_abf_meta)
{
  vector<double> vL10Abfs, vWeights;
  
  // compute the average ABF over all configurations
  for (map<string, double>::iterator it = mAllConfigsAbfs.begin();
       it != mAllConfigsAbfs.end();
       ++it)
    vL10Abfs.push_back (it->second);
  vL10Abfs.push_back (l10_abf_meta);
  vWeights.assign (vL10Abfs.size(), 1.0 / (double) vL10Abfs.size());
  mAllConfigsAbfs["avg.all"] = log10_weighted_sum (&vL10Abfs[0],
						   &vWeights[0],
						   vL10Abfs.size());
  
  // compute the average ABF over the consistent configuration
  // and each one-subgroup-eQTL configuration
  vL10Abfs.clear();
  vWeights.clear();
  for (map<string, double>::iterator it = mAllConfigsAbfs.begin();
       it != mAllConfigsAbfs.end();
       ++it)
    if (it->first.length() == 1)
      vL10Abfs.push_back (it->second);
  vL10Abfs.push_back (l10_abf_meta);
  vWeights.assign (vL10Abfs.size(), 1.0 / (double) vL10Abfs.size());
  mAllConfigsAbfs["avg.subset"] = log10_weighted_sum (&vL10Abfs[0],
						      &vWeights[0],
						      vL10Abfs.size());
}

/** \brief Compute the log10 of the three approximate Bayes factors for 
 *  a given pair feature-SNP, weighted over a grid of hyperparameter values.
 */
void
computeAbfsForPairOverGrid (
  const vector< vector<double> > & grid, 
  const map<size_t, vector<double> > & mAllSumStats,
  double & l10_abf_meta,
  double & l10_abf_fix,
  double & l10_abf_maxh,
  size_t & nbSubgroupsUsed,
  size_t & nbSamplesUsed,
  const size_t & nbSubgroups,
  const bool & calcAllConfigs,
  gsl_combination * c,
  map<string, double> & mAllConfigsAbfs,
  const int & verbose)
{
  size_t subgroup = 0, gridIdx = 0;
  vector< vector<double> > mAllSumStatsCorr;
  double l10_abf_metas[grid.size()], l10_abf_fixs[grid.size()], 
    l10_abf_maxhs[grid.size()], weights[grid.size()];
  
  // apply small sample size correction for the summary stats of each study,
  // and keep the one being defined
  nbSamplesUsed = 0;
  for (subgroup = 0; subgroup < mAllSumStats.size(); ++subgroup)
  {
    vector<double> vSumStats = mAllSumStats.find(subgroup)->second;  // can't use [] because map is const
    vector<double> vSumStatsCorr = correctSummaryStats (vSumStats, verbose-1);
    if (vSumStatsCorr.size() > 0)  // otherwise, pair not present in subgroup
    {
      mAllSumStatsCorr.push_back (vSumStatsCorr);
      nbSamplesUsed += (int) vSumStats[0];
    }
  }
  nbSubgroupsUsed = mAllSumStatsCorr.size();
  
  // calculate the 3 ABFs over the grid
  for (gridIdx = 0; gridIdx < grid.size(); ++gridIdx)
  {
    l10_abf_metas[gridIdx] = getAbfFromStdSumStats (mAllSumStatsCorr,
						    grid[gridIdx][0],
						    grid[gridIdx][1],
						    verbose-1);
    l10_abf_fixs[gridIdx] = getAbfFromStdSumStats (mAllSumStatsCorr,
						   0,
						   grid[gridIdx][0]
						   + grid[gridIdx][1],
						   verbose-1);
    l10_abf_maxhs[gridIdx] = getAbfFromStdSumStats (mAllSumStatsCorr,
						    grid[gridIdx][0]
						    + grid[gridIdx][1],
						    0,
						    verbose-1);
    weights[gridIdx] = 1.0 / (double) grid.size();
  }
  
  // calculate the 3 weighted ABFs
  l10_abf_meta = log10_weighted_sum (l10_abf_metas, weights, grid.size());
  l10_abf_fix = log10_weighted_sum (l10_abf_fixs, weights, grid.size());
  l10_abf_maxh = log10_weighted_sum (l10_abf_maxhs, weights, grid.size());
  if (verbose > 0)
    printf("subgroups=%zu samples=%zu l10_abf_meta=%e l10_abf_fix=%e l10_abf_maxh=%e\n",
	   nbSubgroupsUsed, nbSamplesUsed, l10_abf_meta, l10_abf_fix,
	   l10_abf_maxh);
  
  if (calcAllConfigs)
  {
    computeAbfMetaOfAllConfigsForPairOverGrid (mAllConfigsAbfs,
					       nbSubgroupsUsed,
					       c,
					       mAllSumStatsCorr,
					       grid);
    computeAveragesAbfMetaOverConfigs (mAllConfigsAbfs, l10_abf_meta);
    if (verbose > 0)
      printf ("configs=%zu avg.all=%f avg.subset=%f\n", mAllConfigsAbfs.size(),
	      mAllConfigsAbfs["avg.all"], mAllConfigsAbfs["avg.subset"]);
  }
}

void
writeOutputResultsForPair (
  ofstream & outStream,
  const string & pairId,
  const string & snpIdx,
  const size_t & nbSubgroupsUsed,
  const size_t & nbSamplesUsed,
  const double & l10_abf_meta,
  const double & l10_abf_fix,
  const double & l10_abf_maxh,
  const size_t & nbSubgroups,
  const bool & calcAllConfigs,
  gsl_combination * c,
  const map<string, double> & mAllConfigsAbfs,
  const bool & saveSumStats,
  const map<size_t, vector<double> > & mAllSumStats)
{
  size_t subgroup, i, j;
  vector<string> tokens;
  stringstream ss;
  map<string, double>::const_iterator it;
  
  split (pairId, '|', tokens);
  outStream << tokens[0];
  
  if (snpIdx.compare("id") == 0)
    outStream << " " << tokens[1] << " --";
  else if (snpIdx.compare("coord") == 0)
    outStream << " -- " << tokens[1];
  else if (snpIdx.compare("id+coord") == 0)
    outStream << " " << tokens[1] << " " << tokens[2];
  
  outStream << " " << nbSubgroupsUsed
	    << " " << nbSamplesUsed
	    << " " << l10_abf_meta
	    << " " << l10_abf_fix
	    << " " << l10_abf_maxh;
  
  if (calcAllConfigs)
  {
    outStream << " " << mAllConfigsAbfs.find("avg.all")->second
	      << " " << mAllConfigsAbfs.find("avg.subset")->second;
    for(i=1; i<nbSubgroups; ++i)
    {
      c = gsl_combination_calloc (nbSubgroups, i);
      if (c == NULL)
      {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit (1);
      }
      while (true)
      {
	ss.str("");
	for (j=0; j<c->k; ++j)
	  ss << gsl_combination_get (c, j) + 1;
	it = mAllConfigsAbfs.find(ss.str());
	if (it != mAllConfigsAbfs.end())
	  outStream << " " << it->second;
	else
	  outStream << " " << numeric_limits<double>::quiet_NaN();
	if (gsl_combination_next (c) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free (c);
    }
    
  }
  
  if (saveSumStats)
    for (subgroup = 0; subgroup < nbSubgroups; ++subgroup)
      outStream << " " << mAllSumStats.find(subgroup)->second[1]  // betahat
		<< " " << mAllSumStats.find(subgroup)->second[2]  // sebetahat
		<< " " << mAllSumStats.find(subgroup)->second[3]; // sigmahat
  
  outStream << endl;
}

/** \brief Loop over all pairs to compute ABFs and write results.
 *  \note Use a string feature|coord for the pair id.
 */
void
computeAndWriteAbfsForAllPairs (
  const vector<string> & vInFiles,
  const map<string, vector<size_t> > & mPairs2Positions,
  const vector< vector<double> > & grid,
  const string & outFile,
  const string & snpIdx,
  const bool & saveSumStats,
  const bool & calcAllConfigs,
  const int & verbose)
{
  size_t subgroup = 0, i = 0, nbPairs = 0, nbSubgroups = vInFiles.size(),
    nbSubgroupsUsed = 0, nbSamplesUsed = 0;
  vector<ifstream *> vPtInStreams;
  map<string, vector<size_t> >::const_iterator it_mP2P;
  string pairId;
  vector<string> tokens;
  vector<size_t> vPositions;
  map<size_t, vector<double> > mAllSumStats;  // 1 vector of doubles per subgroup
  double l10_abf_meta = 0, l10_abf_fix = 0, l10_abf_maxh = 0;
  ofstream outStream;
  vector<size_t> vNbSubgroups2Occurrences (nbSubgroups+1, 0);
  vector<size_t> vCounters = getCounters (mPairs2Positions.size());
  gsl_combination * c = NULL;
  map<string, double> mAllConfigsAbfs;
  
  if (verbose > 0)
  {
    cout << "compute the approximate Bayes factors ..." << endl;
    fflush (stdout);
  }
  
  // open the input files once for all
  for (subgroup = 0; subgroup < nbSubgroups; ++subgroup)
  {
    ifstream * pt_inStream = new ifstream;
    pt_inStream->open (vInFiles[subgroup].c_str());
    if (! pt_inStream->good())
    {
      cerr << "ERROR: can't open file '" << vInFiles[subgroup] << "' for reading." << endl;
      exit (1);
    }
    vPtInStreams.push_back (pt_inStream);
  }
  
  // prepare the output file
  outStream.open (outFile.c_str());
  if (! outStream.good())
  {
    cerr << "ERROR: can't open file '" << outFile << "' for writing." << endl;
    exit (1);
  }
  writeOutputHeader (outStream, calcAllConfigs, c, nbSubgroups, saveSumStats);
  
  // for each pair feature-SNP
  for (it_mP2P = mPairs2Positions.begin();
       it_mP2P != mPairs2Positions.end();
       ++it_mP2P)
  {
    ++nbPairs;
    if (verbose > 0)
      printCounter (nbPairs, vCounters);
    pairId = it_mP2P->first;
    vPositions = it_mP2P->second;
    if (verbose > 1)
      cout << "#" << nbPairs << " " << pairId << endl;
    
    // retrieve its summary statistics
    getSummaryStatsForPair (pairId, snpIdx, vPositions, vPtInStreams,
			    mAllSumStats, nbSubgroupsUsed, nbSamplesUsed,
			    verbose-1);
    ++vNbSubgroups2Occurrences[nbSubgroupsUsed];
    
    // compute the Bayes Factors for the 3 models over the grid
    computeAbfsForPairOverGrid (grid, mAllSumStats, l10_abf_meta,
				l10_abf_fix, l10_abf_maxh,
				nbSubgroupsUsed, nbSamplesUsed,
				nbSubgroups, calcAllConfigs, c,
				mAllConfigsAbfs, verbose-1);
    
    // write the results
    writeOutputResultsForPair (outStream, pairId, snpIdx, nbSubgroupsUsed,
			       nbSamplesUsed, l10_abf_meta, l10_abf_fix,
			       l10_abf_maxh, nbSubgroups, calcAllConfigs, c,
			       mAllConfigsAbfs, saveSumStats, mAllSumStats);
  }
  
  // close all files
  outStream.close();
  for (subgroup = 0; subgroup < nbSubgroups; ++subgroup)
  {
    vPtInStreams[subgroup]->close();
    delete vPtInStreams[subgroup];
  }
  vPtInStreams.clear();
  
  if (verbose > 0)
  {
    cout << "results written in " << outFile << endl;
    int sumPerc = 0;
    int perc = 0;
    for(i = 0; i < vNbSubgroups2Occurrences.size()-1; ++i)
    {
      if (vNbSubgroups2Occurrences[i] != 0)
      {
	perc = (int) round (100 * vNbSubgroups2Occurrences[i] / mPairs2Positions.size());
	sumPerc += perc;
	cout << "nb of pairs gene-SNP with " << i
	     << " subgroup" << ((i < 2) ? "" : "s")
	     << ": " << vNbSubgroups2Occurrences[i]
	     << " (" << perc << "%)" << endl;
      }
    }
    i = vNbSubgroups2Occurrences.size() - 1;
    cout << "nb of pairs gene-SNP with " << i
	 << " subgroups: " << vNbSubgroups2Occurrences[i]
	 << " (" << (int) 100 - sumPerc
	 << "%)" << endl;
  }
}

int main (int argc, char ** argv)
{
  string inDir, gridFile, outFile, ftrsFile, snpsFile, snpIdx;
  bool saveSumStats = false, calcAllConfigs = false;
  int verbose = 1;
  parse_args (argc, argv, inDir, gridFile, outFile,
	      ftrsFile, snpsFile, snpIdx, saveSumStats, calcAllConfigs,
	      verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")" << endl;
  }
  
  vector<string> vInFiles;
  if (isDirectory(inDir.c_str()))
    vInFiles = scanInputDirectory (inDir, verbose);
  else
    vInFiles.push_back (inDir);
  vector< vector<double> > grid = loadGrid (gridFile, verbose);
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsFile, verbose);
  
  map<string, vector<size_t> > mPairs2Positions;
  indexInputFiles (vInFiles,
		   vFtrsToKeep,
		   vSnpsToKeep,
		   snpIdx,
		   mPairs2Positions,
		   verbose);
  
  computeAndWriteAbfsForAllPairs (vInFiles,
				  mPairs2Positions,
				  grid,
				  outFile,
				  snpIdx,
				  saveSumStats,
				  calcAllConfigs,
				  verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime)
	 << ": elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << ")" << endl;
  }
  
  return 0;
}
