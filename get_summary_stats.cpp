/** \file get_summary_stats.cpp
 *
 *  `get_summary_stats' computes summary statistics of association between genotypes and phenotypes.
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
 *  gcc -Wall -g -DDEBUG get_summary_stats.cpp -lstdc++ -lgzstream -lz -lgsl -lgslcblas -o get_summary_stats
 *  help2man -o get_summary_stats.man ./get_summary_stats
 *  groff -mandoc get_summary_stats.man > get_summary_stats.ps
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
#include <algorithm>
#include <limits>
using namespace std;

#include <gzstream.h>
#include <gsl/gsl_cdf.h>

#include "utils.cpp"

//-----------------------------------------------------------------------------
// OOP-like definitions of two structures and their associated functions

struct SnpStats
{
  string name;  // eg. rs2205177
  size_t coord;
  double betahat;
  double sebetahat;
  double sigmahat;
  double log10_pval;
};

struct FtrStats
{
  string name;  // eg. ENSG00000182816
  vector <SnpStats> vSnpStats;
};

void FtrStats_reset (FtrStats * pt_iFtrStats)
{
  pt_iFtrStats->name.clear();
  pt_iFtrStats->vSnpStats.clear();
}

void FtrStats_write (FtrStats iFtrStats, long int n, ostream & outStream)
{
  SnpStats iSnpStats;
  for (size_t snp_id = 0; snp_id < iFtrStats.vSnpStats.size(); snp_id++)
  {
    iSnpStats = iFtrStats.vSnpStats[snp_id];
    if (iSnpStats.name.empty())
      continue;
    outStream << iFtrStats.name
	      << " " << iSnpStats.name
	      << " " << iSnpStats.coord
	      << " " << n
	      << " " << iSnpStats.betahat
	      << " " << iSnpStats.sebetahat
	      << " " << iSnpStats.sigmahat
	      << " " << iSnpStats.log10_pval
	      << endl;
  }
}

//-----------------------------------------------------------------------------
// list of all functions

void help (char ** argv);
void version (char ** argv);
void parse_args (int argc, char ** argv,
		 string * pt_linksFile,
		 string * pt_genoFile,
		 string * pt_phenoFile,
		 string * pt_outFile,
		 bool * pt_shouldOutputFileBeGzipped,
		 string * pt_ftrsFile,
		 string * pt_snpsFile,
		 int * pt_verbose);
map<string, vector<string> >
loadLinksFtr2Snps (const string linksFile,
		   vector<string> vFtrsToKeep,
		   vector<string> vSnpsToKeep,
		   int verbose);
map<string, long int>
indexFtrs (const string phenoFile,
	   const map<string, vector<string> > ftrName2CisSnpNameCoords,
	   int verbose);
map<string, long int>
indexSnps (const string genoFile,
	   vector<string> vSnpsToKeep,
	   int verbose);
void getPhenoValues (istream & phenoStream,
		     long int ftrPos,
		     double ** y,
		     size_t * n_y,
		     int ** idx);
void getGenoValues (istream & genoStream,
		    long int snpPos,
		    double ** g,
		    size_t * n_g,
		    const int * idx);
void ols (const size_t n, const double * g, const double * y,
	  double * betahat, double * sebetahat, double * sigmahat,
	  double * log10_pval, int verbose);
void
computeSummaryStatsForOneFeature (
  const map<string, vector<string> >::iterator ftrName2CisSnpNameCoords_it,
  map<string, long int> snpNameCoord2Pos,
  ifstream & genoStream,
  const double * y,
  const size_t n_y,
  const int * idx,
  FtrStats * pt_iFtrStats,
  const int verbose);
void computeAndWriteSummaryStatsFtrPerFtr (
  map<string, vector<string> > ftrName2CisSnpNameCoords,
  map<string, long int> ftrName2Pos,
  map<string, long int> snpNameCoord2Pos,
  string phenoFile,
  string genoFile,
  string outFile,
  int verbose);
int main (int argc, char ** argv);

//-----------------------------------------------------------------------------
// functions

/** \brief Display the usage on stdout.
*/
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " assesses associations between genotypes and phenotypes" << endl
       << "(one genetic variant per phenotype) and returns the summary"
       << " statistics betahat," << endl
       << "se(betahat) and sigmahat, as well as the P-value"
       << " for H0:\"beta=0\"." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS]..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (default=1)" << endl
       << "  -l, --links\tgziped file with links <gene><space/tab><SNP|coord>" << endl
       << "  -g, --geno\tfile with genotypes in IMPUTE format" << endl
       << "  -p, --pheno\tfile with phenotypes (row 1 for sample names,"
       << " column 1" << endl
       << "\t\tfor gene names, delimiter=<space/tab>)" << endl
       << "  -o, --output\toutput file (will be gzipped if '.gz' in file name)" << endl
       << "  -f, --ftr\tgzipped file with a list of features to analyze" << endl
       << "\t\t(one feature name per line)" << endl
       << "  -s, --snp\tgzipped file with a list of SNPs to analyze" << endl
       << "\t\t(one SNP coordinate per line)" << endl
       << endl
       << "Examples:" << endl
       << "  " << argv[0] << " -l <links> -g <genotypes> -p <phenotypes> -o <output>" << endl;
}

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
		 string * pt_linksFile,
		 string * pt_genoFile,
		 string * pt_phenoFile,
		 string * pt_outFile,
		 bool * pt_shouldOutputFileBeGzipped,
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
	{"links", required_argument, 0, 'l'},
	{"geno", required_argument, 0, 'g'},
	{"pheno", required_argument, 0, 'p'},
	{"output", required_argument, 0, 'o'},
	{"ftr", required_argument, 0, 'f'},
	{"snp", required_argument, 0, 's'}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:l:g:p:o:f:s:",
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
    case 'l':
      *pt_linksFile = optarg;
      break;
    case 'g':
      *pt_genoFile = optarg;
      break;
    case 'p':
      *pt_phenoFile = optarg;
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
  if ((*pt_linksFile).empty())
  {
    fprintf (stderr, "ERROR: missing file with links feature-SNPs.\n");
    help (argv);
    exit (1);
  }
  if ((*pt_genoFile).empty())
  {
    fprintf (stderr, "ERROR: missing file with genotypes.\n");
    help (argv);
    exit (1);
  }
  if ((*pt_phenoFile).empty())
  {
    fprintf (stderr, "ERROR: missing file with phenotypes.\n");
    help (argv);
    exit (1);
  }
  if ((*pt_outFile).empty())
  {
    fprintf (stderr, "ERROR: missing output file.\n");
    help (argv);
    exit (1);
  }
  if (pt_outFile->find(".gz") != string::npos)
  {
    *pt_shouldOutputFileBeGzipped = true;
  }
  else
  {
    *pt_shouldOutputFileBeGzipped = false;
  }
}

/** \brief Load the 3-column file feature<tab>snp<tab>coord into a map 
 *  which keys are feature names and values are vectors of snp-coord.
 */
map<string, vector<string> >
loadLinksFtr2Snps (const string linksFile,
		   vector<string> vFtrsToKeep,
		   vector<string> vSnpsToKeep,
		   int verbose)
{
  map<string, vector<string> > ftrName2CisSnpNameCoords;
  igzstream linksStream;
  string line;
  vector<string> tokens;
  size_t line_id = 0;
  
  linksStream.open(linksFile.c_str());
  if (! linksStream.good())
  {
    cerr << "ERROR: can't open file " << linksFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout <<"load file " << linksFile << "..." << endl;
  }
  
  while (linksStream.good())
  {
    getline (linksStream, line);
    if (line.empty())
    {
      break;
    }
    if (line.find('\t') != string::npos)
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    if (tokens.size() != 2)
    {
      cerr << line << endl;
      cerr << "ERROR: format of file " << linksFile
	   << " should be feature<space/tab>snp" << endl;
      exit (1);
    }
    if (! vFtrsToKeep.empty() &&
	find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0]) == vFtrsToKeep.end())
    {
      continue;
    }
    if (! vSnpsToKeep.empty() &&
	find(vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[1]) == vSnpsToKeep.end())
    {
      continue;
    }
    line_id++;
    if (ftrName2CisSnpNameCoords.find(tokens[0]) == ftrName2CisSnpNameCoords.end())
    {
      ftrName2CisSnpNameCoords.insert( make_pair(tokens[0], vector<string>()) );
    }
    ftrName2CisSnpNameCoords[tokens[0]].push_back (tokens[1]);
  }
  
  linksStream.close();
  
  if (verbose > 0)
  {
    cout << "features with cis SNPs: " << ftrName2CisSnpNameCoords.size() << endl
	 << "pairs feature-SNP: " << line_id << endl;
  }
  
  return ftrName2CisSnpNameCoords;
}

/** \brief Index the file with the phenotype values, ie. a big matrix with
 *  samples in columns and features in row.
 *  Only for the features having SNPs in cis.
 */
map<string, long int>
indexFtrs (const string phenoFile,
	   const map<string, vector<string> > ftrName2CisSnpNameCoords,
	   int verbose)
{
  map<string, long int> ftrName2Pos;
  ifstream phenoStream;
  string line;
  vector<string> tokens;
  size_t nbSamples = 0;
  long int ftrPos;
  
  phenoStream.open(phenoFile.c_str());
  if (! phenoStream.is_open())
  {
    cerr << "ERROR: can't open file " << phenoFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout << "index file " << phenoFile << "..." << endl;
  }
  
  // header line (sample names)
  getline (phenoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  nbSamples = tokens.size();
  
  while (phenoStream.good())
  {
    ftrPos = phenoStream.tellg();
    getline (phenoStream, line);
    if (line.empty())
    {
      break;
    }
    if (line.find('\t') != string::npos)
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    if (tokens.size() != nbSamples+1)
    {
      cerr << "ERROR: different number of samples for feature "
	   << tokens[0] << " (" << tokens.size() << " vs " << nbSamples+1 << ")" << endl;
      exit (1);
    }
    if (ftrName2CisSnpNameCoords.find(tokens[0]) == ftrName2CisSnpNameCoords.end())
    {
      continue;
    }
    ftrName2Pos.insert( make_pair(tokens[0], ftrPos) );
  }
  
  phenoStream.close();
  
  if (ftrName2Pos.size() == 0)
  {
    cerr << "WARNING: no feature has phenotypic value" << endl;
    exit (0);
  }
  if (verbose > 0)
  {
    cout << "features with values: " << ftrName2Pos.size() << endl;
  }
  
  return ftrName2Pos;
}

/** \brief Index the file with the genotype values (IMPUTE2 format).
 */
map<string, long int>
indexSnps (const string genoFile,
	   vector<string> vSnpsToKeep,
	   int verbose)
{
  map<string, long int> snpNameCoord2Pos;
  ifstream genoStream;
  string line;
  vector<string> tokensLine;
  vector<string> tokensSnp;
  long int snpPos;
  
  genoStream.open(genoFile.c_str());
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout << "index file " << genoFile << "..." << endl;
  }
  
  while (genoStream.good())
  {
    snpPos = genoStream.tellg();
    getline (genoStream, line);
    if (line.empty())
    {
      break;
    }
    if (line.find('\t') != string::npos)
      split (line, '\t', tokensLine);
    else
      split (line, ' ', tokensLine);
    stringstream ss;
    ss << tokensLine[1] << "|" << tokensLine[2];  // SNP_name|SNP_coord
    if (! vSnpsToKeep.empty() &&
	find(vSnpsToKeep.begin(), vSnpsToKeep.end(), ss.str()) == vSnpsToKeep.end())
    {
      continue;
    }
    snpNameCoord2Pos.insert( make_pair(ss.str(), snpPos) );
  }
  
  genoStream.close();
  
  if (verbose > 0)
  {
    cout << "SNPs with values: " << snpNameCoord2Pos.size() << endl;
  }
  
  return snpNameCoord2Pos;
}

/** \brief Retrieve phenotype values for a given feature
 *  but assume no missing values.
 */
void getPhenoValues (istream & phenoStream,
		     long int ftrPos,
		     double ** y,
		     size_t * n_y)
{
  string line;
  vector<string> tokens;
  phenoStream.seekg(ftrPos);
  getline (phenoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  *n_y = tokens.size() - 1;
  *y = (double *) malloc (*n_y * sizeof(double *));
  for (size_t i = 0; i < *n_y; i++)
  {
    (*y)[i] = atof(tokens[i+1].c_str());
  }
}

/** \brief Retrieve phenotype values for a given feature
 *  and skip missing values encoded as NA.
 */
void getPhenoValues (istream & phenoStream,
		     long int ftrPos,
		     double ** y,
		     size_t * n_y,
		     int ** idx)
{
  size_t i = 0, j = 0;
  string line, tok;
  vector<string> tokens;
  phenoStream.seekg(ftrPos);
  getline (phenoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  *n_y = tokens.size() - 1;
  *y = (double *) malloc (*n_y * sizeof(double *));
  *idx = (int *) malloc (*n_y * sizeof(int *));
  for (i = 0; i < *n_y; ++i)
  {
    tok = tokens[i+1];
    if (tok.compare("NA") == 0)
    {
      (*idx)[i] = 0;
    }
    else
    {
      (*idx)[i] = 1;
      (*y)[j] = atof(tok.c_str());
      ++j;
    }
  }
  if (j != *n_y)
  {
    *y = (double *) realloc (*y, j);
    *n_y = j;
  }
}

/** \brief Retrieve genotype values for a given SNP
 *  and assume no sample has missing phenotype.
 */
void getGenoValues (istream & genoStream,
		    long int snpPos,
		    double ** g,
		    size_t * n_g)
{
  string line;
  vector<string> tokens;
  genoStream.seekg(snpPos);
  getline (genoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  *n_g = (tokens.size() - 5) / 3;
  *g = (double *) malloc (*n_g * sizeof(double *));
  for (size_t i = 0; i < *n_g; i++)
  {
    (*g)[i] = 0 * atof(tokens[5+3*i].c_str())
      + 1 * atof(tokens[5+3*i+1].c_str())
      + 2 * atof(tokens[5+3*i+2].c_str());
  }
}

/** \brief Retrieve genotype values for a given SNP
 *  and skip samples having a missing phenotype.
 */
void getGenoValues (istream & genoStream,
		    long int snpPos,
		    double ** g,
		    size_t * n_g,
		    const int * idx)
{
  string line;
  vector<string> tokens;
  genoStream.seekg(snpPos);
  getline (genoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  *n_g = (tokens.size() - 5) / 3;
  *g = (double *) malloc (*n_g * sizeof(double *));
  for (size_t i = 0; i < *n_g; i++)
  {
    if (idx[i] == 1)
    {
      (*g)[i] = 0 * atof(tokens[5+3*i].c_str())
	+ 1 * atof(tokens[5+3*i+1].c_str())
	+ 2 * atof(tokens[5+3*i+2].c_str());
    }
  }
}

/** \brief Compute the summary statistics of the linear regression.
 */
void ols (const size_t n, const double * g, const double * y,
	  double * betahat, double * sebetahat, double * sigmahat,
	  double * log10_pval, int verbose)
{
  size_t i = 0;
  double ym = 0, gm = 0, yty = 0, gtg = 0, gty = 0;
  for(i=0; i<n; ++i){
    ym += y[i];
    gm += g[i];
    yty += y[i] * y[i];
    gtg += g[i] * g[i];
    gty += g[i] * y[i];
  }
  ym /= n;
  gm /= n;
  double denom = gtg - n * gm * gm;
#ifdef DEBUG
  if (verbose > 0)
    printf ("ym=%f gm=%f yty=%f gtg=%f gty=%f denom=%f\n", ym, gm, yty, gtg,
	    gty, denom);
#endif
  if(denom > 1e-8)
  {
    *betahat = (gty - n * gm * ym) / denom;
    double inv_xtx[2][2];
    inv_xtx[0][0] = gtg;
    inv_xtx[0][1] = - n * gm;
    inv_xtx[1][0] = - n * gm;
    inv_xtx[1][1] = n;
    double xty[2];
    xty[0] = n * ym;
    xty[1] = gty;
    double rss1 = yty - 1/denom * (n*ym*(gtg*ym - gm*gty) - gty*(n*gm*ym - gty));
    if (fabs(*betahat) > 1e-8)
      *sigmahat = sqrt(rss1 / (n-2));
    else  // case where phenotypes are not variable enough
      *sigmahat = sqrt((yty - n * ym * ym) / (n-2));  // sqrt(rss0/(n-2))
    *sebetahat = *sigmahat / sqrt(gtg - n*gm*gm);
    double muhat = (ym*gtg - gm*gty) / (gtg - n*gm*gm);
    double mss = 0;
    for(i=0; i<n; ++i)
      mss += pow(muhat + *betahat * g[i] - ym, 2);
    *log10_pval = gsl_cdf_fdist_Q (mss/pow(*sigmahat,2), 1, n-2);
  }
  else
  {
    if (verbose > 0)
      cout << "genotypes are not variable enough" << endl;
    *betahat = numeric_limits<double>::infinity();  // or use quiet_NaN?
    *sebetahat = numeric_limits<double>::infinity();
    *sigmahat = numeric_limits<double>::infinity();
    *log10_pval = numeric_limits<double>::infinity();
  }
}

/** \brief First, read all genotypes of all cis SNPs of the given feature.
 *   Second, compute the summary stats for each pair gene-SNP.
 */
void
computeSummaryStatsForOneFeature (
  const map<string, vector<string> >::iterator ftrName2CisSnpNameCoords_it,
  map<string, long int> snpNameCoord2Pos,
  ifstream & genoStream,
  const double * y,
  const size_t n_y,
  const int * idx,
  FtrStats * pt_iFtrStats,
  const int verbose)
{
  vector<string> cisSnpNameCoords = (*ftrName2CisSnpNameCoords_it).second;
  size_t snp_id;
  map<string, double *> mSnpNameCoord2Genotypes;
  
  // for each SNP in cis of the given feature,
  // retrieve its genotypes and keep them in a map
  for (snp_id = 0; snp_id < cisSnpNameCoords.size(); snp_id++)
  {
    string snpNameCoord = cisSnpNameCoords[snp_id];
    if (verbose > 0)
      cout << "#" << snp_id+1 << "/" << cisSnpNameCoords.size()
	   << " " << snpNameCoord << endl;
    if (snpNameCoord2Pos.find(snpNameCoord) == snpNameCoord2Pos.end())
      continue;
    
    size_t snpPos = snpNameCoord2Pos[snpNameCoord];
    double * g;
    size_t n_g;
    getGenoValues (genoStream, snpPos, &g, &n_g, idx);
    if (n_y != n_g)
    {
      cerr << "ERROR: different number of samples for feature "
	   << pt_iFtrStats->name << " (" << n_y << ")  and SNP "
	   << snpNameCoord << " (" << n_g << ")" << endl;
      exit (1);
    }
    if (verbose > 0)
      cout << n_g << " values: " << g[0] << " " << g[1] << " ..." << endl;
    
    mSnpNameCoord2Genotypes.insert( make_pair(snpNameCoord, g) );
  }
  
  // for each SNP in cis of the given feature
  // compute the OLS summary stats and keep them
  map<string, double *>::iterator mSnpNameCoord2Genotypes_it;
  for (mSnpNameCoord2Genotypes_it = mSnpNameCoord2Genotypes.begin();
       mSnpNameCoord2Genotypes_it != mSnpNameCoord2Genotypes.end();
       mSnpNameCoord2Genotypes_it++)
  {
    SnpStats iSnpStats;
    
    string snpNameCoord = mSnpNameCoord2Genotypes_it->first;
    vector<string> tokens;
    split (snpNameCoord, '|', tokens);
    iSnpStats.name = tokens[0];
    iSnpStats.coord = atol(tokens[1].c_str());
    
    double * g = mSnpNameCoord2Genotypes_it->second;
    ols (n_y, g, y, &iSnpStats.betahat, &iSnpStats.sebetahat,
	 &iSnpStats.sigmahat, &iSnpStats.log10_pval, verbose-1);
    if (verbose > 0)
      printf ("betahat=%.8f sebetahat=%.8f sigmahat=%.8f P-value=%.8f\n",
	      iSnpStats.betahat, iSnpStats.sebetahat, iSnpStats.sigmahat,
	      iSnpStats.log10_pval);
    
    pt_iFtrStats->vSnpStats.push_back (iSnpStats);
    free (g);
  }
}

/** \brief Compute OLS summary statistics for each pair feature-SNP
 *  and write the results feature by feature in a gzipped file.
 */
void computeAndGzWriteSummaryStatsFtrPerFtr (
  map<string, vector<string> > ftrName2CisSnpNameCoords,
  map<string, long int> ftrName2Pos,
  map<string, long int> snpNameCoord2Pos,
  string phenoFile,
  string genoFile,
  string outFile,
  int verbose)
{
  FtrStats iFtrStats;
  map<string, vector<string> >::iterator ftrName2CisSnpNameCoords_it;
  vector<string> cisSnpNameCoords, tokens;
  string ftrName, snpNameCoord, snpName;
  size_t ftrPos, n_y = 0;
  double nbFtrs = 0, step = floor(ftrName2CisSnpNameCoords.size() / 5);
  double * y = NULL;
  int * idx = NULL;
  ifstream phenoStream, genoStream;
  ogzstream outStream;
  
  // open all files
  phenoStream.open(phenoFile.c_str(), ifstream::in);
  if (! phenoStream.is_open())
  {
    cerr << "ERROR: can't open file " << phenoFile << endl;
    exit (1);
  }
  genoStream.open(genoFile.c_str(), ifstream::in);
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  outStream.open(outFile.c_str());
  if (! outStream.good())
  {
    cerr << "ERROR: can't open file " << outFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout << "compute summary statistics for each pair feature-SNP..." << endl;
  }
  
  // for each feature
  for (ftrName2CisSnpNameCoords_it = ftrName2CisSnpNameCoords.begin();
       ftrName2CisSnpNameCoords_it != ftrName2CisSnpNameCoords.end();
       ftrName2CisSnpNameCoords_it++)
  {
    nbFtrs++;
    if (verbose > 0 && fmod(nbFtrs, step) == 0.0)
      cout << nbFtrs << "/" << ftrName2CisSnpNameCoords.size() << endl;
    ftrName = (*ftrName2CisSnpNameCoords_it).first;
    if (verbose > 1)
      cout << "#" << nbFtrs << "/" << ftrName2CisSnpNameCoords.size()
	   << " " << ftrName << endl;
    FtrStats_reset (&iFtrStats);
    iFtrStats.name = ftrName;
    
    // retrieve its values
    ftrPos = ftrName2Pos[ftrName];
    getPhenoValues (phenoStream, ftrPos, &y, &n_y, &idx);
    if (verbose > 1)
      cout << n_y << " values: " << y[0] << " " << y[1] << " ..." << endl;
    
    // loop over SNPs in cis
    computeSummaryStatsForOneFeature (ftrName2CisSnpNameCoords_it,
				      snpNameCoord2Pos,
				      genoStream,
				      y,
				      n_y,
				      idx,
				      &iFtrStats,
				      verbose-2);
    
    // write the results (one line per SNP)
    FtrStats_write (iFtrStats, n_y, outStream);
  }
  
  phenoStream.close();
  genoStream.close();
  outStream.close();
}

/** \brief Compute OLS summary statistics for each pair feature-SNP
 *  and write the results feature by feature in an uncompressed file.
 */
void computeAndWriteSummaryStatsFtrPerFtr (
  map<string, vector<string> > ftrName2CisSnpNameCoords,
  map<string, long int> ftrName2Pos,
  map<string, long int> snpNameCoord2Pos,
  string phenoFile,
  string genoFile,
  string outFile,
  int verbose)
{
  FtrStats iFtrStats;
  map<string, vector<string> >::iterator ftrName2CisSnpNameCoords_it;
  vector<string> cisSnpNameCoords, tokens;
  string ftrName, snpNameCoord, snpName;
  size_t ftrPos, n_y = 0;
  double nbFtrs = 0, step = floor(ftrName2CisSnpNameCoords.size() / 5);
  double * y = NULL;
  int * idx = NULL;
  ifstream phenoStream, genoStream;
  ofstream outStream;
  
  // open all files
  phenoStream.open(phenoFile.c_str(), ifstream::in);
  if (! phenoStream.is_open())
  {
    cerr << "ERROR: can't open file " << phenoFile << endl;
    exit (1);
  }
  genoStream.open(genoFile.c_str(), ifstream::in);
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  outStream.open(outFile.c_str());
  if (! outStream.good())
  {
    cerr << "ERROR: can't open file " << outFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout << "compute summary statistics for each pair feature-SNP..." << endl;
  }
  
  // for each feature
  for (ftrName2CisSnpNameCoords_it = ftrName2CisSnpNameCoords.begin();
       ftrName2CisSnpNameCoords_it != ftrName2CisSnpNameCoords.end();
       ++ftrName2CisSnpNameCoords_it)
  {
    nbFtrs++;
    if (verbose > 0 && fmod(nbFtrs, step) == 0.0)
      printf ("%.0f%%\n", 100 * nbFtrs / ftrName2CisSnpNameCoords.size());
    ftrName = (*ftrName2CisSnpNameCoords_it).first;
    if (ftrName2Pos.find(ftrName) == ftrName2Pos.end())
      continue;
    if (verbose > 1)
      cout << "#" << nbFtrs << "/" << ftrName2CisSnpNameCoords.size()
	   << " " << ftrName << endl;
    FtrStats_reset (&iFtrStats);
    iFtrStats.name = ftrName;
    
    // retrieve its values
    ftrPos = ftrName2Pos[ftrName];
    getPhenoValues (phenoStream, ftrPos, &y, &n_y, &idx);
    if (verbose > 1)
      cout << n_y << " values: " << y[0] << " " << y[1] << " ..." << endl;
    
    // loop over SNPs in cis
    computeSummaryStatsForOneFeature (ftrName2CisSnpNameCoords_it,
				      snpNameCoord2Pos,
				      genoStream,
				      y,
				      n_y,
				      idx,
				      &iFtrStats,
				      verbose-2);
    
    // write the results (one line per SNP)
    FtrStats_write (iFtrStats, n_y, outStream);
    
    free (y);
    free (idx);
  }
  
  phenoStream.close();
  genoStream.close();
  outStream.close();
}

int main (int argc, char ** argv)
{
  string linksFile, genoFile, phenoFile, outFile, ftrsFile, snpsFile;
  bool shouldOutputFileBeGzipped;
  int verbose = 1;
  parse_args (argc, argv, &linksFile, &genoFile, &phenoFile, &outFile,
	      &shouldOutputFileBeGzipped, &ftrsFile, &snpsFile, &verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")" << endl;
  }
  
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsFile, verbose);
  
  map<string, vector<string> > ftrName2CisSnpNameCoords
    = loadLinksFtr2Snps (linksFile, vFtrsToKeep, vSnpsToKeep, verbose);
  
  map<string, long int> ftrName2Pos
    = indexFtrs (phenoFile, ftrName2CisSnpNameCoords, verbose);
  map<string, long int> snpNameCoord2Pos
    = indexSnps (genoFile, vSnpsToKeep, verbose);
  
  if (shouldOutputFileBeGzipped)
    computeAndGzWriteSummaryStatsFtrPerFtr (ftrName2CisSnpNameCoords,
					    ftrName2Pos,
					    snpNameCoord2Pos,
					    phenoFile,
					    genoFile,
					    outFile,
					    verbose);
  else
    computeAndWriteSummaryStatsFtrPerFtr (ftrName2CisSnpNameCoords,
					  ftrName2Pos,
					  snpNameCoord2Pos,
					  phenoFile,
					  genoFile,
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
