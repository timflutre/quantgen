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
  *  gcc -Wall -g -DDEBUG get_abf_meta.cpp -lstdc++ -lgzstream -lz -lgsl -lgslcblas -o get_abf_meta
  *  help2man -o get_abf_meta.man ./get_abf_meta
  *  groff -mandoc get_abf_meta.man > get_abf_meta.ps
  */

#include <cmath>
#include <ctime>
#include <dirent.h>
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

#include <gzstream.h>
#include <gsl/gsl_cdf.h>

#include "utils.cpp"

//-----------------------------------------------------------------------------

struct Pair
{
  string feature;
  string snp;
  string coord;
  string id;
  bool operator< (const Pair & rhs) const 
  {return id < rhs.id;}
};

void Pair_set (Pair & iPair, string feature, string snp, string coord)
{
  iPair.feature = feature;
  iPair.snp = snp;
  iPair.coord = coord;
  stringstream ss;
  ss << iPair.feature << "|" << iPair.coord;
  iPair.id = ss.str();
}

// http://stackoverflow.com/q/8234253/597069
struct checkPairId : unary_function<pair<Pair, vector<size_t> >, bool>
{
private:
  string _exp;
public:
  checkPairId (Pair iPair) : _exp(iPair.id) {}
  bool operator() (const pair<Pair, vector<size_t> > & x) const
  {
    return x.first.id.compare(_exp) == 0;
  }
};

//-----------------------------------------------------------------------------

/** \brief Display the help on stdout.
 */
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       <<" performs the meta-analysis by computing three Bayes Factors," << endl
       << "one per model (average, fixed-effects, max heterogeneity)," << endl
       << "for each pair feature-SNP (use SNP coord rather than name)." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS]..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (default=1)" << endl
       << "  -i, --input\tdirectory with the summary stats for each subgroup" << endl
       << "\t\t(one file per subgroup, with at least 7 columns:\n\t\tftr snp coord n betahat sebetahat sigmahat)" << endl
       << "  -g, --grid\tfile with the grid of values for psi2 and omega2 (ES model)" << endl
       << "  -o, --out\tgzipped output file" << endl
       << "  -f, --ftr\tgzipped file with a list of features to analyze" << endl
       << "  -s, --snp\tgzipped file with a list of SNP coordinates to analyze" << endl
       << "  "
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
void parse_args (int argc, char ** argv,
		 string * pt_inDir,
		 string * pt_gridFile,
		 string * pt_outFile,
		 string * pt_ftrsFile,
		 string * pt_snpsFile,
		 int * pt_verbose)
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
	{"snp", required_argument, 0, 's'}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:i:g:o:f:s:",
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
      *pt_verbose = atoi(optarg);
      break;
    case 'i':
      *pt_inDir = optarg;
      break;
    case 'g':
      *pt_gridFile = optarg;
      break;
    case 'o':
      *pt_outFile = optarg;
      break;
    case 'f':
      *pt_ftrsFile = optarg;
      break;
    case 's':
      *pt_snpsFile = optarg;
      break;
    case '?':
      break;
    default:
      abort ();
    }
  }
  if ((*pt_inDir).empty())
  {
    fprintf (stderr, "ERROR: missing input directory with summary statistics (-i).\n\n");
    help (argv);
    exit (1);
  }
  if ((*pt_gridFile).empty())
  {
    fprintf (stderr, "ERROR: missing grid file (-g).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (*pt_gridFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", pt_gridFile->c_str());
    help (argv);
    exit (1);
  }
  if ((*pt_outFile).empty())
  {
    fprintf (stderr, "ERROR: missing output file (-o).\n\n");
    help (argv);
    exit (1);
  }
}

/** \brief List the input directory with the OLS summary stats.
 */
vector<string> scanInputDirectory (string inDir,
				   int verbose)
{
  vector<string> vInFiles;
  struct dirent ** inFiles = NULL;
  int nbInFiles;
  if (verbose > 0)
  {
    cout << "scan directory " << inDir << "..." << endl;
  }
  nbInFiles = scandir(inDir.c_str(), &inFiles, dummy_selector, alphasort);
  if (nbInFiles == -1)
  {
    cerr << "ERROR: can't scan " << inDir << endl;
    exit (1);
  }
  else if (nbInFiles == 0)
  {
    cerr << "ERROR: " << inDir << " contains no file" << endl;
    exit (1);
  }
  else
  {
    for (int s = 0; s < nbInFiles; ++s)
    {
      if (string(inFiles[s]->d_name) == "." ||
	  string(inFiles[s]->d_name) == "..")
      {
	continue;
      }
      char path[500];
      int nbChar;
      if (inDir[inDir.size()-1] != '/')
	nbChar = sprintf (path, "%s/%s", inDir.c_str(), inFiles[s]->d_name);
      else
	nbChar = sprintf (path, "%s%s", inDir.c_str(), inFiles[s]->d_name);
      if (nbChar < 0)
      {
	cerr << "ERROR: variable 'path' is not big enough" << endl;
      }
      vInFiles.push_back (string(path));
      free (inFiles[s]);
    }
    if (verbose > 0)
      cout << "nb of subgroups: " << vInFiles.size() << endl;
  }
  free (inFiles);
  return vInFiles;
}

/** \brief Index the files with a map to know which pairs feature-SNP are 
 *  in common among subgroups.
 *  \note Use a struct for the pair.
 */
void indexInputFiles (vector<string> vInFiles,
		      vector<string> vFtrsToKeep,
		      vector<string> vSnpsToKeep,
		      map<Pair, vector<size_t> > & mPairs2Positions,
		      int verbose)
{
  size_t subgroup = 0, line_id = 0, pos = 0;
  string line;
  vector<string> tokens;
  map<Pair, vector<size_t> >::iterator it_mP2P;
  
  if (verbose > 0)
  {
    cout << "index input files..." << endl;
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
    
    line_id = 0;
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
      if (tokens.size() < 7)
      {
	cerr << endl << "ERROR: file " << vInFiles[subgroup] 
	     << " at line " << line_id 
	     << " should have at least 7 fields:" << endl
	     << "feature<space>SNP<space>coord<space>n<space>betahat<space>sebetahat<space>sigmahat"
	     << endl;
	exit (1);
      }
      if (! vFtrsToKeep.empty() &&
	  find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0]) == vFtrsToKeep.end())
      {
	continue;  // skip line, based on feature name
      }
      if (! vSnpsToKeep.empty() &&
	  find(vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[2]) == vSnpsToKeep.end())
      {
	continue;  // skip line, based on SNP coordinate
      }
      Pair iPair_obs;
      Pair_set (iPair_obs, tokens[0], tokens[1], tokens[2]);
      it_mP2P = find_if(mPairs2Positions.begin(), mPairs2Positions.end(), checkPairId(iPair_obs));
      if (it_mP2P == mPairs2Positions.end())
      {
	vector<size_t> vPositions (vInFiles.size(), string::npos);
	mPairs2Positions.insert ( make_pair (iPair_obs, vPositions) );
	mPairs2Positions[iPair_obs][subgroup] = pos;
      }
      else
      {
	it_mP2P->second[subgroup] = pos;
	// keep "rs" identifier as SNP name, if available
	if (it_mP2P->first.snp.find("rs") == string::npos
	    && iPair_obs.snp.find("rs") != string::npos)
	{
	  vector<size_t> vPositions;
	  vPositions.swap (it_mP2P->second);
	  mPairs2Positions.erase (it_mP2P);
	  mPairs2Positions.insert ( make_pair (iPair_obs, vPositions) );
	}
      }
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

/** \brief Index the files with a map to know which pairs feature-SNP are 
 *  in common among subgroups.
 *  \note Use a string feature|coord for the pair id.
 */
void indexInputFiles (vector<string> vInFiles,
		      vector<string> vFtrsToKeep,
		      vector<string> vSnpsToKeep,
		      map<string, vector<size_t> > & mPairs2Positions,
		      int verbose)
{
  size_t subgroup = 0, line_id = 0, pos = 0;
  string line;
  vector<string> tokens;
  map<string, vector<size_t> >::iterator it_mP2P;
  
  if (verbose > 0)
  {
    cout << "index input files..." << endl;
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
    
    line_id = 0;
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
      if (tokens.size() < 7)
      {
	cerr << endl << "ERROR: file " << vInFiles[subgroup] 
	     << " at line " << line_id 
	     << " should have at least 7 fields:" << endl
	     << "ftr<space>snp<space>coord<space>n<space>betahat<space>sebetahat<space>sigmahat"
	     << endl;
	exit (1);
      }
      if (tokens[0].compare("ftr") == 0)
      {
	continue;  // skip header line
      }
      if (! vFtrsToKeep.empty() &&
	  find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0]) ==
	  vFtrsToKeep.end())
      {
	continue;  // skip line, based on feature name
      }
      if (! vSnpsToKeep.empty() &&
	  find(vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[2]) ==
	  vSnpsToKeep.end())
      {
	continue;  // skip line, based on SNP coordinate
      }
      stringstream ss;
      ss << tokens[0] << "|" << tokens[2];  // feature|coord
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
vector< vector<double> > loadGrid (string gridFile, int verbose)
{
  if (verbose > 0)
  {
    cout << "load grid..." << endl;
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

/** \brief Retrieve the summary stats (n,betahat,sebetahat,sigmahat) for a 
 *  given pair feature-SNP in each  subgroup in which this pair is present.
 *
 *  \note If the pair is not present in the given subgroup (missing value),
 *  n will be 0 and the other summary stats will be NaN. If the pair is present
 *  but has infinite summary stats (should be coded as Inf in input file), only
 *  n will be an integer, the other summary stats remaining Inf.
 */
void getSummaryStatsForPair (Pair iPair,
			     vector<size_t> vPositions,
			     vector<ifstream *> vPtInStreams,
			     map<size_t, vector<double> > & allSumStats,
			     int * pt_nbSubgroupsUsed,
			     int * pt_nbSamplesUsed,
			     int verbose)
{
  size_t subgroup = 0, colIdx = 0;
  string line;
  vector<string> tokens;
  
  allSumStats.clear();
  *pt_nbSubgroupsUsed = 0;
  *pt_nbSamplesUsed = 0;
  
  for (subgroup = 0; subgroup < vPositions.size(); ++subgroup)
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
      ++(*pt_nbSubgroupsUsed);
      vPtInStreams[subgroup]->seekg (vPositions[subgroup]);
      getline (*(vPtInStreams[subgroup]), line);
      if (line.find('\t') != string::npos)
	split (line, '\t', tokens);
      else
	split (line, ' ', tokens);
      stringstream pairId;
      pairId << tokens[0] << "|" << tokens[2];  // feature|coord
      if (pairId.str() != iPair.id)
      {
	cerr << "ERROR: stream position is wrong for feature " << iPair.feature
	     << " and SNP " << iPair.snp << " (" << iPair.coord << ") "
	     << " in subgroup " << subgroup+1 << endl;
	cerr << line << endl;
	exit (1);
      }
      *pt_nbSamplesUsed += atoi (tokens[3].c_str());
      for (colIdx = 3; colIdx < 7; ++colIdx)
      {
	vSumStats[colIdx-3] = atof (tokens[colIdx].c_str());
      }
    }
    
    allSumStats.insert (make_pair(subgroup, vSumStats) );
  }
}

/** \brief Retrieve the summary stats (n,betahat,sebetahat,sigmahat) for a 
 *  given pair feature-SNP in each  subgroup in which this pair is present.
 *
 *  \note If the pair is not present in the given subgroup (missing value),
 *  n will be 0 and the other summary stats will be NaN. If the pair is present
 *  but has infinite summary stats (should be coded as Inf in input file), only
 *  n will be an integer, the other summary stats remaining Inf.
 */
void getSummaryStatsForPair (string pairId,
			     vector<size_t> vPositions,
			     vector<ifstream *> vPtInStreams,
			     map<size_t, vector<double> > & allSumStats,
			     int * pt_nbSubgroupsUsed,
			     int * pt_nbSamplesUsed,
			     int verbose)
{
  size_t subgroup = 0, colIdx = 0;
  string line;
  vector<string> tokens;
  
  allSumStats.clear();
  *pt_nbSubgroupsUsed = 0;
  *pt_nbSamplesUsed = 0;
  
  for (subgroup = 0; subgroup < vPositions.size(); ++subgroup)
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
      ++(*pt_nbSubgroupsUsed);
      vPtInStreams[subgroup]->seekg (vPositions[subgroup]);
      getline (*(vPtInStreams[subgroup]), line);
      if (line.find('\t') != string::npos)
	split (line, '\t', tokens);
      else
	split (line, ' ', tokens);
      stringstream ss;
      ss << tokens[0] << "|" << tokens[2];  // feature|coord
      if (ss.str() != pairId)
      {
	cerr << "ERROR: stream position is wrong for pair " << pairId
	     << " in subgroup " << subgroup+1 << endl;
	cerr << line << endl;
	exit (1);
      }
      *pt_nbSamplesUsed += atoi (tokens[3].c_str());
      for (colIdx = 3; colIdx < 7; ++colIdx)
      {
	vSumStats[colIdx-3] = atof (tokens[colIdx].c_str());
      }
    }
    
    allSumStats.insert (make_pair(subgroup, vSumStats) );
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
vector<double> correctSummaryStats (const vector<double> sumStats, int verbose)
{
  vector<double> stdSumStatsCorr;
  double n = 0, bhat = 0, sebhat = 0, sigmahat = 0, t = 0;
  
  if (sumStats[0] != 0)  // if=0 -> pair is not present in the given subgroup
  {
    n = sumStats[0];
    bhat = sumStats[1] / sumStats[3];      // before correction
    sebhat = sumStats[2] / sumStats[3];    // before correction
#ifdef DEBUG
    if (verbose > 0)
      printf ("n=%i betahat=%e sebetahat=%e sigmahat=%e bhat=%e sebhat=%e\n",
	      (int) n, sumStats[1], sumStats[2], sumStats[3], bhat, sebhat);
#endif
    t = gsl_cdf_gaussian_Pinv(gsl_cdf_tdist_P(-fabs(bhat/sebhat), n-2), 1.0);
#ifdef DEBUG
    if (verbose > 0)
      printf ("t=%.8f\n", t);
#endif
    if (fabs(t) > 1e-8)
    {
      sigmahat = fabs(sumStats[1]) / (fabs(t) * sebhat);
      bhat = sumStats[1] / sigmahat;
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
double getAbfFromStdSumStats (vector< vector<double> > stdSumStats,
			      double phi2,
			      double oma2,
			      int verbose)
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
double log10_weighted_sum (double * vec, double * weights, size_t size)
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

/** \brief Compute the log10 of the three approximate Bayes factors for 
 *  a given pair feature-SNP, weighted over a grid of hyperparameter values.
 */
void computeAbfsForPairOverGrid (vector< vector<double> > grid, 
				 const map<size_t, vector<double> > & allSumStats,
				 double * pt_l10_abf_meta,
				 double * pt_l10_abf_fix,
				 double * pt_l10_abf_maxh,
				 int * pt_nbSubgroupsUsed,
				 int * pt_nbSamplesUsed,
				 int verbose)
{
  size_t subgroup = 0, gridIdx = 0;
  vector< vector<double> > allSumStatsCorr;
  double l10_abf_metas[grid.size()], l10_abf_fixs[grid.size()], 
    l10_abf_maxhs[grid.size()], weights[grid.size()];
  
  *pt_l10_abf_meta = 0;
  *pt_l10_abf_fix = 0;
  *pt_l10_abf_maxh = 0;
  
  // apply small sample size correction for the summary stats of each study,
  // and keep the one being defined
  *pt_nbSamplesUsed = 0;
  for (subgroup = 0; subgroup < allSumStats.size(); ++subgroup)
  {
    vector<double> sumStats = allSumStats.find(subgroup)->second;  // can't use [] because map is const
    vector<double> sumStatsCorr = correctSummaryStats (sumStats, verbose-1);
    if (sumStatsCorr.size() > 0)  // otherwise, pair not present in subgroup
    {
      allSumStatsCorr.push_back (sumStatsCorr);
      *pt_nbSamplesUsed += (int) sumStats[0];
    }
  }
  *pt_nbSubgroupsUsed = allSumStatsCorr.size();
  
  // calculate the 3 ABFs over the grid
  for (gridIdx = 0; gridIdx < grid.size(); ++gridIdx)
  {
    l10_abf_metas[gridIdx] = getAbfFromStdSumStats (allSumStatsCorr,
						grid[gridIdx][0],
						grid[gridIdx][1],
						verbose-1);
    l10_abf_fixs[gridIdx] = getAbfFromStdSumStats (allSumStatsCorr,
					       0,
					       grid[gridIdx][0]
					       + grid[gridIdx][1],
					       verbose-1);
    l10_abf_maxhs[gridIdx] = getAbfFromStdSumStats (allSumStatsCorr,
						grid[gridIdx][0]
						+ grid[gridIdx][1],
						0,
						verbose-1);
    weights[gridIdx] = 1.0 / (double) grid.size();
  }
  
  // calculate the 3 weighted ABFs
  *pt_l10_abf_meta = log10_weighted_sum (l10_abf_metas, weights, grid.size());
  *pt_l10_abf_fix = log10_weighted_sum (l10_abf_fixs, weights, grid.size());
  *pt_l10_abf_maxh = log10_weighted_sum (l10_abf_maxhs, weights, grid.size());
  if (verbose > 0)
    printf("subgroups=%i samples=%i l10_abf_meta=%e l10_abf_fix=%e l10_abf_maxh=%e\n",
	   *pt_nbSubgroupsUsed, *pt_nbSamplesUsed, *pt_l10_abf_meta, *pt_l10_abf_fix,
	   *pt_l10_abf_maxh);
}

/** \brief Loop over all pairs to compute ABFs and write results.
 */
void computeAndWriteAbfsForAllPairs (vector<string> vInFiles,
				     const map<Pair, vector<size_t> > & mPairs2Positions,
				     string gridFile,
				     string outFile,
				     int verbose)
{
  size_t subgroup = 0, i = 0;
  vector<ifstream *> vPtInStreams;
  map<Pair, vector<size_t> >::const_iterator it;
  Pair iPair;
  vector<size_t> vPositions;
  map<size_t, vector<double> > allSumStats;  // 1 vector of doubles per subgroup
  double l10_abf_meta = 0, l10_abf_fix = 0, l10_abf_maxh = 0,
    nbPairs = 0, step = floor(mPairs2Positions.size() / 5);
  int nbSamplesUsed = 0, nbSubgroupsUsed = 0;
  ogzstream outStream;
  vector<size_t> vNbSubgroups2Occurrences (vInFiles.size()+1, 0);
  
  if (verbose > 0)
  {
    cout << "compute the approximate Bayes factors..." << endl;
    fflush (stdout);
  }
  
  // open the input files once for all
  for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
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
  outStream << "ftr"
	    << " snp"
	    << " coord"
	    << " nb.subgroups"
	    << " nb.samples"
	    << " l10abf.avg"
	    << " l10abf.fix"
	    << " l10abf.maxh";
  for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
  {
    outStream << " betahat.s" << subgroup+1
	      << " sebetahat.s" << subgroup+1
	      << " sigmahat.s" << subgroup+1;
  }
  outStream << endl;
  
  vector< vector<double> > grid = loadGrid (gridFile, verbose);
  
  // for each pair feature-SNP
  for (it = mPairs2Positions.begin();
       it != mPairs2Positions.end();
       ++it)
  {
    ++nbPairs;
    if (verbose > 0 && fmod (nbPairs, step) == 0.0)
    {
      printf ("%.0f%%\n", 100 * nbPairs / mPairs2Positions.size());
      fflush (stdout);
    }
    iPair = it->first;
    vPositions = it->second;
    
    getSummaryStatsForPair (iPair, vPositions, vPtInStreams, allSumStats,
			    &nbSubgroupsUsed, &nbSamplesUsed, verbose-1);
    ++vNbSubgroups2Occurrences[nbSubgroupsUsed];
    if (verbose > 1)
      cout << "#" << nbPairs << " " << iPair.feature
	   << " " << iPair.snp << " (" << iPair.coord << "):"
	   << " " << nbSubgroupsUsed << " subgroups" << endl;
    
    computeAbfsForPairOverGrid (grid, allSumStats, &l10_abf_meta, &l10_abf_fix, &l10_abf_maxh,
				&nbSubgroupsUsed, &nbSamplesUsed, verbose-1);
    
    outStream << iPair.feature
	      << " " << iPair.snp
	      << " " << iPair.coord
	      << " " << nbSubgroupsUsed
	      << " " << nbSamplesUsed
	      << " " << l10_abf_meta
	      << " " << l10_abf_fix
	      << " " << l10_abf_maxh;
    for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
    {
      outStream << " " << allSumStats[subgroup][1]  // betahat
		<< " " << allSumStats[subgroup][2]  // sebetahat
		<< " " << allSumStats[subgroup][3]; // sigmahat
    }
    outStream << endl;
  }
  
  outStream.close();
  for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
  {
    vPtInStreams[subgroup]->close();
    delete vPtInStreams[subgroup];
  }
  vPtInStreams.clear();
  if (verbose > 0)
  {
    cout << "results written in " << outFile << endl;
    for(i = 0; i < 2; ++i)
    {
      cout << "nb of pairs gene-SNP with " << i
	   << " subgroup: " << vNbSubgroups2Occurrences[i] << endl;
    }
    for(i = 2; i < vNbSubgroups2Occurrences.size(); ++i)
    {
      cout << "nb of pairs gene-SNP with " << i
	   << " subgroups: " << vNbSubgroups2Occurrences[i] << endl;
    }
  }
}

/** \brief Loop over all pairs to compute ABFs and write results.
 *  \note Use a string feature|coord for the pair id.
 */
void computeAndWriteAbfsForAllPairs (vector<string> vInFiles,
				     const map<string, vector<size_t> > & mPairs2Positions,
				     string gridFile,
				     string outFile,
				     int verbose)
{
  size_t subgroup = 0, i = 0, nbPairs = 0;
  vector<ifstream *> vPtInStreams;
  map<string, vector<size_t> >::const_iterator it_mP2P;
  string pairId;
  vector<string> tokens;
  vector<size_t> vPositions;
  map<size_t, vector<double> > allSumStats;  // 1 vector of doubles per subgroup
  double l10_abf_meta = 0, l10_abf_fix = 0, l10_abf_maxh = 0;
  int nbSamplesUsed = 0, nbSubgroupsUsed = 0;
  ogzstream outStream;
  vector<size_t> vNbSubgroups2Occurrences (vInFiles.size()+1, 0);
  vector<size_t> vCounters = getCounters (mPairs2Positions.size());
  
  if (verbose > 0)
  {
    cout << "compute the approximate Bayes factors..." << endl;
    fflush (stdout);
  }
  
  // open the input files once for all
  for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
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
  outStream << "ftr"
	    << " snp"
	    << " coord"
	    << " nb.subgroups"
	    << " nb.samples"
	    << " l10abf.avg"
	    << " l10abf.fix"
	    << " l10abf.maxh";
  for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
  {
    outStream << " betahat.s" << subgroup+1
	      << " sebetahat.s" << subgroup+1
	      << " sigmahat.s" << subgroup+1;
  }
  outStream << endl;
  
  vector< vector<double> > grid = loadGrid (gridFile, verbose);
  
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
    
    // retrieve its summary statistics
    getSummaryStatsForPair (pairId, vPositions, vPtInStreams, allSumStats,
			    &nbSubgroupsUsed, &nbSamplesUsed, verbose-1);
    ++vNbSubgroups2Occurrences[nbSubgroupsUsed];
    if (verbose > 1)
      cout << "#" << nbPairs << " " << pairId
	   << " " << nbSubgroupsUsed << " subgroups" << endl;
    
    // compute the Bayes Factors for the 3 models over the grid
    computeAbfsForPairOverGrid (grid, allSumStats, &l10_abf_meta, &l10_abf_fix, &l10_abf_maxh,
				&nbSubgroupsUsed, &nbSamplesUsed, verbose-1);
    
    // write the results
    split (pairId, '|', tokens);
    outStream << tokens[0]
	      << " " << "--"
	      << " " << tokens[1]
	      << " " << nbSubgroupsUsed
	      << " " << nbSamplesUsed
	      << " " << l10_abf_meta
	      << " " << l10_abf_fix
	      << " " << l10_abf_maxh;
    for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
    {
      outStream << " " << allSumStats[subgroup][1]  // betahat
		<< " " << allSumStats[subgroup][2]  // sebetahat
		<< " " << allSumStats[subgroup][3]; // sigmahat
    }
    outStream << endl;
  }
  
  // close all files
  outStream.close();
  for (subgroup = 0; subgroup < vInFiles.size(); ++subgroup)
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
      perc = (int) round (100 * vNbSubgroups2Occurrences[i] / mPairs2Positions.size());
      sumPerc += perc;
      cout << "nb of pairs gene-SNP with " << i
	   << " subgroup" << ((i < 2) ? "" : "s")
	   << ": " << vNbSubgroups2Occurrences[i]
	   << " (" << perc << "%)" << endl;
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
  string inDir, gridFile, outFile, ftrsFile, snpsFile;
  int verbose = 1;
  parse_args (argc, argv, &inDir, &gridFile, &outFile,
	      &ftrsFile, &snpsFile, &verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")" << endl;
  }
  
  vector<string> vInFiles = scanInputDirectory (inDir, verbose);
  
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsFile, verbose);
  
//  map<Pair, vector<size_t> > mPairs2Positions;
  map<string, vector<size_t> > mPairs2Positions;
  indexInputFiles (vInFiles,
		   vFtrsToKeep,
		   vSnpsToKeep,
		   mPairs2Positions,
		   verbose);
  
  computeAndWriteAbfsForAllPairs (vInFiles,
				  mPairs2Positions,
				  gridFile,
				  outFile,
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
