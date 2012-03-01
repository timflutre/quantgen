/** \file get_summary_stats.cpp
 *
 *  `get_summary_stats' computes summary statistics of association between genotypes and phenotypes.
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "utils.cpp"

const gsl_rng_type * T;
gsl_rng * r;

//-----------------------------------------------------------------------------

struct SnpStats
{
  string name;  // eg. rs2205177
  size_t coord;
  double maf;
  double betahat;
  double sebetahat;
  double sigmahat;
  double pval;
  double R2;  // coef of determination, proportion of variance explained
  double rs;
  double rsZscore;  // using Fisher transformation
  double rsPvalZtest;
  double rsPvalPerms;
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
	      << " " << iSnpStats.maf
	      << " " << n
	      << " " << iSnpStats.betahat
	      << " " << iSnpStats.sebetahat
	      << " " << iSnpStats.sigmahat
	      << " " << iSnpStats.pval
	      << " " << iSnpStats.R2;
    if (iSnpStats.rs != numeric_limits<double>::quiet_NaN())
    {
      outStream << " " << iSnpStats.rs
		<< " " << iSnpStats.rsZscore
		<< " " << iSnpStats.rsPvalZtest
		<< " " << iSnpStats.rsPvalPerms;
    }
    outStream << endl;
  }
}

//-----------------------------------------------------------------------------

/** \brief Display the help on stdout.
*/
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " assesses associations between genotypes and phenotypes" << endl
       << "(one genetic variant per phenotype) and returns the summary"
       << " statistics betahat," << endl
       << "se(betahat) and sigmahat, as well as the P-value"
       << " for H0:\"beta=0\" and the PVE." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS]..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (default=1)" << endl
       << "  -l, --links\tgzipped file with links <gene><space/tab><SNP|coord>" << endl
       << "\t\t(especially useful to focus on genetic variants in cis)" << endl
       << "  -g, --geno\tfile with genotypes in IMPUTE format (delimiter=<space/tab>)" << endl
       << "\t\tsamples in columns should be in same order as the phenotype file" << endl
       << "  -p, --pheno\tfile with phenotypes (row 1 for sample names,"
       << " column 1" << endl
       << "\t\tfor feature names, delimiter=<space/tab>)" << endl
       << "  -o, --output\toutput file for the summary stats" << endl
       << "  -f, --ftr\tgzipped file with a list of features to analyze" << endl
       << "\t\t(one feature name per line)" << endl
       << "  -s, --snp\tgzipped file with a list of SNPs to analyze" << endl
       << "\t\t(one SNP coordinate per line)" << endl
       << "  -d, --discard\tgzipped file with a list of individuals to discard" << endl
       << "\t\t(one individual per line, should match header of phenotype file)" << endl
       << "  -m, --maf\tthreshold for the minor allele frequency (default=0.0)" << endl
       << "\t\t(whatever the option, the MAF will still be computed and saved)" << endl
       << "  -c, --cor\tnumber of permutations to assess significance of the Spearman" << endl
       << "\t\trank correlation coefficient" << endl
       << "\t\t(default c<0, Spearman coef is not computed;" << endl
       << "\t\tif c=0, add Z-score and P-value for Z-test;" << endl
       << "\t\tif c>0, add P-value for c permutations of phenotype labels)" << endl
       << endl
       << "Examples:" << endl
       << "  " << argv[0] << " -l <links> -g <genotypes> -p <phenotypes> -o <output>" << endl
       << endl
       << "Remarks:" << endl
       << "  Samples with missing phenotypes (NA) are skipped, but missing genotypes are" << endl
       << "forbidden: consider imputing them first. For non-variable genotypes, a zero" << endl
       << "effect size is returned, along with an infinite std error and a P-value of 1." << endl;
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
		 string * pt_linksFile,
		 string * pt_genoFile,
		 string * pt_phenoFile,
		 string * pt_outFile,
		 string * pt_ftrsFile,
		 string * pt_snpsFile,
		 string * pt_indsFile,
		 double * pt_minMaf,
		 int * pt_nbPermutations,
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
	{"snp", required_argument, 0, 's'},
	{"discard", required_argument, 0, 'd'},
	{"maf", required_argument, 0, 'm'},
	{"cor", no_argument, 0, 'c'}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:l:g:p:o:f:s:d:m:c:",
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
    case 'd':
      *pt_indsFile = optarg;
      break;
    case 'm':
      *pt_minMaf = atof(optarg);
      break;
    case 'c':
      *pt_nbPermutations = atoi(optarg);
      break;
    case '?':
      printf ("\n"); help (argv);
      abort ();
    default:
      printf ("\n"); help (argv);
      abort ();
    }
  }
  if ((*pt_linksFile).empty())
  {
    fprintf (stderr, "ERROR: missing file with links feature-SNPs (-l).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (*pt_linksFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", pt_linksFile->c_str());
    help (argv);
    exit (1);
  }
  if ((*pt_genoFile).empty())
  {
    fprintf (stderr, "ERROR: missing file with genotypes (-g).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (*pt_genoFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", pt_genoFile->c_str());
    help (argv);
    exit (1);
  }
  if ((*pt_phenoFile).empty())
  {
    fprintf (stderr, "ERROR: missing file with phenotypes (-p).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (*pt_phenoFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", pt_phenoFile->c_str());
    help (argv);
    exit (1);
  }
  if ((*pt_outFile).empty())
  {
    fprintf (stderr, "ERROR: missing output file (-o).\n\n");
    help (argv);
    exit (1);
  }
  if (*pt_minMaf < 0 || *pt_minMaf > 1)
  {
    fprintf (stderr, "ERROR: min MAF should be between 0 and 1 (-m).\n\n");
    help (argv);
    exit (1);
  }
}

/** \brief Load the 2-column file feature<tab>snp|coord into a map 
 *  which keys are feature names and values are vectors of snp-coord.
 */
map<string, vector<string> >
loadLinksFtr2Snps (const string linksFile,
		   const vector<string> vFtrsToKeep,
		   const vector<string> vSnpsToKeep,
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
	   << " should be feature<space/tab>snp|coord" << endl;
      exit (1);
    }
    if (! vFtrsToKeep.empty() &&
	find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0]) ==
	vFtrsToKeep.end())
    {
      continue;
    }
    if (! vSnpsToKeep.empty())
    {
      vector<string> tokens2 = split (tokens[1], '|');
      if (find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens2[1]) ==
	  vSnpsToKeep.end())
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
 *  Only for the individuals not discarded.
 */
map<string, long int>
indexFtrs (const string phenoFile,
	   const map<string, vector<string> > ftrName2CisSnpNameCoords,
	   const vector<string> vIndsToSkip,
	   vector<size_t> & vIdxIndsToSkip,
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
  for(size_t i = 0; i < tokens.size(); ++i)
    if(vIndsToSkip.size() > 0 & find(vIndsToSkip.begin(), vIndsToSkip.end(),
				     tokens[i]) != vIndsToSkip.end())
      vIdxIndsToSkip.push_back(i);
  
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
	   << tokens[0] << " (" << tokens.size() << " vs " << nbSamples+1
	   << ")" << endl;
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
    if(vIdxIndsToSkip.size() > 0)
      cout << "individuals to discard: " << vIdxIndsToSkip.size() << endl;
  }
  
  return ftrName2Pos;
}

/** \brief Return the minor allele frequency.
 *  \note The input comes from a line in the IMPUTE format that was splitted.
 */
double getMaf (const vector<string> & tokens,
	       const vector<size_t> & vIdxIndsToSkip)
{
  double maf = 0;
  for (size_t i = 0; i < (tokens.size() - 5) / 3; ++i)
  {
    if(vIdxIndsToSkip.size() == 0 | find(vIdxIndsToSkip.begin(),
					 vIdxIndsToSkip.end(),
					 i) == vIdxIndsToSkip.end())
      maf += 1 * atof(tokens[5+3*i+1].c_str())
	+ 2 * atof(tokens[5+3*i+2].c_str());
  }
  maf /= 2 * ((tokens.size() - 5) / 3 - vIdxIndsToSkip.size());
  return maf <= 0.5 ? maf : (1 - maf);
}

/** \brief Index the file with the genotype values (IMPUTE format).
 */
map<string, long int>
indexSnps (const string genoFile,
	   vector<string> vSnpsToKeep,
	   const vector<size_t> vIdxIndsToSkip,
	   const double minMaf,
	   int verbose)
{
  map<string, long int> snpNameCoord2Pos;
  ifstream genoStream;
  string line;
  vector<string> tokens;
  long int snpPos;
  double maf;
  
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
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    if (! vSnpsToKeep.empty() &&
	find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[2])
	== vSnpsToKeep.end())
    {
      continue;
    }
    if (minMaf > 0)
    {
      maf = getMaf (tokens, vIdxIndsToSkip);
#ifdef DEBUG
      if (verbose > 1)
	cout << "SNP " << tokens[1] << "|" << tokens[2] << ": MAF="
	     << maf << endl;
#endif
      if (maf < minMaf)
      {
#ifdef DEBUG
	if (verbose > 1)
	  cout << "SNP " << tokens[1] << "|" << tokens[2] <<
	    ": skip because MAF < " << minMaf << endl;
#endif
	continue;
      }
    }
    stringstream ss;
    ss << tokens[1] << "|" << tokens[2];  // SNP_name|SNP_coord
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
 *  and skip missing values encoded as NA.
 *  \note isNa[i] == true if sample i has no phenotype
 */
void getPhenoValues (const string ftrName,
		     istream & phenoStream,
		     long int ftrPos,
		     vector<double> & y,
		     const vector<size_t> vIdxIndsToSkip,
		     vector<bool> & isNa,
		     int verbose)
{
  string line, tok;
  vector<string> tokens;
  isNa.clear();
  phenoStream.seekg(ftrPos);
  getline (phenoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  for (size_t i = 1; i < tokens.size(); ++i)
  {
    if(vIdxIndsToSkip.size() > 0 &
       find(vIdxIndsToSkip.begin(), vIdxIndsToSkip.end(), i-1) !=
       vIdxIndsToSkip.end())
    {
      isNa.push_back (false);
      continue;
    }
    tok = tokens[i];
    if (tok.compare("NA") == 0)
    {
      isNa.push_back (true);
    }
    else
    {
      isNa.push_back (false);
      y.push_back (atof(tok.c_str()));
    }
  }
#ifdef DEBUG
  if (verbose > 0)
    printf ("%s samples=%zu to-skip=%zu missing=%zu remaining=%zu %.4f %.4f ...\n",
	    ftrName.c_str(), tokens.size()-1, vIdxIndsToSkip.size(),
	    (size_t) count(isNa.begin(), isNa.end(), true), y.size(), y[0], y[1]);
#endif
}

/** \brief Retrieve genotype values for a given SNP
 *  and skip samples having a missing phenotype.
 *  \note isNa[i] == true if sample i has no phenotype
 */
void getGenoValues (const string ftrName,
		    const string snpNameCoord,
		    istream & genoStream,
		    long int snpPos,
		    vector<double> & g,
		    const vector<size_t> vIdxIndsToSkip,
		    const vector<bool> & isNa,
		    double * pt_maf,
		    int verbose)
{
  size_t i = 0;
  string line;
  vector<string> tokens;
  genoStream.seekg(snpPos);
  getline (genoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  *pt_maf = getMaf (tokens, vIdxIndsToSkip);
  for (i = 0; i < (tokens.size() - 5) / 3; ++i)
  {
    if (vIdxIndsToSkip.size() > 0 &
	find(vIdxIndsToSkip.begin(), vIdxIndsToSkip.end(), i) !=
	vIdxIndsToSkip.end())
      continue;
    if (! isNa[i])
    {
      g.push_back (0 * atof(tokens[5+3*i].c_str())
		   + 1 * atof(tokens[5+3*i+1].c_str())
		   + 2 * atof(tokens[5+3*i+2].c_str()));
    }
    else
      cout << "isNa[" << i << "]=true" << endl;
  }
#ifdef DEBUG
  if (verbose > 0)
    printf ("%s %s samples=%zu to-skip=%zu remaining=%zu %.4f %.4f ...\n",
	    ftrName.c_str(), snpNameCoord.c_str(), (tokens.size()-5)/3,
	    vIdxIndsToSkip.size(), g.size(), g[0], g[1]);
#endif
}

/** \brief Compute the summary statistics of the linear regression.
 *  \note phenotype = mu + genotype * beta + error
 *  \note missing values should have been already filtered out
 */
void ols (const string yName, const string xName,
	  const vector<double> & g, const vector<double> & y,
	  double * betahat, double * sebetahat, double * sigmahat,
	  double * pval, double * R2, int verbose)
{
  size_t i = 0, n = g.size();
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
  double vg = gtg - n * gm * gm;  // variance of the genotypes
#ifdef DEBUG
  if (verbose > 0)
    printf ("%s %s n=%zu ym=%f gm=%f yty=%f gtg=%f gty=%f vg=%f\n",
	    yName.c_str(), xName.c_str(), n, ym, gm, yty, gtg,
	    gty, vg);
#endif
  if(vg > 1e-8)
  {
    *betahat = (gty - n * gm * ym) / vg;
    double inv_xtx[2][2];
    inv_xtx[0][0] = gtg;
    inv_xtx[0][1] = - n * gm;
    inv_xtx[1][0] = - n * gm;
    inv_xtx[1][1] = n;
    double xty[2];
    xty[0] = n * ym;
    xty[1] = gty;
    double rss1 = yty - 1/vg * (n*ym*(gtg*ym - gm*gty) - gty*(n*gm*ym - gty));
    if (fabs(*betahat) > 1e-8)
      *sigmahat = sqrt(rss1 / (n-2));
    else  // case where phenotypes are not variable enough
      *sigmahat = sqrt((yty - n * ym * ym) / (n-2));  // sqrt(rss0/(n-2))
    *sebetahat = *sigmahat / sqrt(gtg - n*gm*gm);
    double muhat = (ym*gtg - gm*gty) / (gtg - n*gm*gm);
    double mss = 0;
    for(i=0; i<n; ++i)
      mss += pow(muhat + *betahat * g[i] - ym, 2);
    *pval = gsl_cdf_fdist_Q (mss/pow(*sigmahat,2), 1, n-2);
    *R2 = mss / (mss + rss1);
  }
  else
  {
    if (verbose > 0)
      cout << "genotypes are not variable enough" << endl;
    *betahat = 0;
    *sebetahat = numeric_limits<double>::infinity();
    *sigmahat = sqrt((yty - n * ym * ym) / (n-2));  // sqrt(rss0/(n-2))
    *pval = 1;
    *R2 = 0;
  }
}

/** \brief Resolve a sequence of ties.
 *  The input ranks array is expected to take the same value for all indices in
 *  tiesTrace. The common value is recoded with the average of the indices. For
 *  example, if ranks = <5,8,2,6,2,7,1,2> and tiesTrace = <2,4,7>, the result
 *  will be <5,8,3,6,3,7,1,3>.
 *
 *  http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.312
 * 
 *  @param ranks array of ranks
 *  @param tiesTrace vector of indices where ranks is constant, that is,
 *  for any i and j in tiesTrace, ranks[i] == ranks[j].
 */
void
resolveTies (double * ranks, vector<size_t> tiesTrace)
{
  // constant value of ranks over tiesTrace
  double c = ranks[tiesTrace[0]];
  
  // length of sequence of tied ranks
  size_t length = tiesTrace.size();
  
  // new rank (ie. the average of the current indices)
  double mean = (2*c + length - 1) / 2;
  
  for(size_t i=0; i<tiesTrace.size(); ++i)
    ranks[tiesTrace[i]] = mean;
}

/** \brief Rank an array using the natural ordering on doubles, ties being 
 *  resolved by taking their average.
 *  \note See http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.190
*/
void
rank (double * ranks, const double data[], const size_t stride, const size_t n)
{
  // copy the input data and sort them
  double * ds = (double *) calloc (n, sizeof(double));
  for(size_t i=0; i<n; ++i)
    ds[i] = data[i];
  gsl_sort (ds, 1, n);
  
  // get the index of the input data as if they were sorted
  size_t * p = (size_t*) calloc (n, sizeof(size_t));
  gsl_sort_index (p, data, stride, n);
  
  // walk the sorted array, filling output array using sorted positions,
  // resolving ties as we go
  size_t pos = 1;
  ranks[p[0]] = pos;
  vector<size_t> tiesTrace;
  tiesTrace.push_back(p[0]);
  for(size_t i=1; i<n; ++i)
  {
    if(ds[i] - ds[i-1] > 0)
    {
      pos = i + 1;
      if(tiesTrace.size() > 1)
	resolveTies (ranks, tiesTrace);
      tiesTrace.clear();
      tiesTrace.push_back(p[i]);
    }
    else
      tiesTrace.push_back(p[i]);
    ranks[p[i]] = pos;
  }
  if(tiesTrace.size() > 1)
    resolveTies (ranks, tiesTrace);
  
  free (p);
  free (ds);
}

/** \brief Compute the Spearman rank correlation coefficient between two arrays
 *  by computing the Pearson correlation coefficient of their ranks.
 *  \note Ties are resolved by taking the mean of their ranks.
*/
double
my_stats_correlation_spearman (const double data1[], const size_t stride1,
			       const double data2[], const size_t stride2,
			       const size_t n)
{
  double rs = 0.0;
  
  double * ranks1 = (double*) calloc (n, sizeof(size_t));
  rank (ranks1, data1, stride1, n);
  
  double * ranks2 = (double*) calloc (n, sizeof(size_t));
  rank (ranks2, data2, stride2, n);
  
  rs = gsl_stats_correlation((double*) ranks1, 1, (double*) ranks2, 1, n );
  
  free (ranks1);
  free (ranks2);
  
  return rs;
}

/** \brief Compute the Spearman rank correlation coefficient between g and y,
 *  and assess its significance with a t test or permutations of y (all other
 *  summary stats are nan).
 *  \note missing values should have been already filtered out
 */
void spearman (const string yName, const string xName,
	       const vector<double> & g, const vector<double> & y,
	       double * rs, double * z, double * pvalZtest,
	       int nbPermutations, double * pvalPerms, int verbose)
{
  gsl_vector_const_view gsl_g = gsl_vector_const_view_array (&g[0],
							     g.size());
  gsl_vector_const_view gsl_y = gsl_vector_const_view_array (&y[0],
							     y.size());
  *rs = my_stats_correlation_spearman (gsl_g.vector.data, 1,
				       gsl_y.vector.data, 1,
				       g.size());
  
  *z = sqrt((g.size() - 3) / 1.06) * 1/2 * (log(1 + *rs) - log(1 - *rs));
  *pvalZtest = gsl_cdf_ugaussian_Q (*z);
  
  if (nbPermutations > 0)
  {
    size_t i;
    double yPerm[y.size()];
    double rsPerm;
    *pvalPerms = 1;
    for(i=0; i<y.size(); ++i)
      yPerm[i] = y[i];
    for(int p=0; p<nbPermutations; ++p)
    {
      gsl_ran_shuffle (r, yPerm, y.size(), sizeof(double));
      rsPerm = my_stats_correlation_spearman (gsl_g.vector.data, 1,
					      yPerm, 1, g.size());
      if (abs(rsPerm) >= abs(*rs))
	++(*pvalPerms);
    }
    *pvalPerms /= (nbPermutations + 1);
  }
  else
    *pvalPerms = numeric_limits<double>::quiet_NaN();
  
#ifdef DEBUG
  if (verbose > 0)
    printf ("%s %s n=%zu rs=%f z=%f pvalZtest=%f pvalPerms=%f\n",
	    yName.c_str(), xName.c_str(), g.size(),
	    *rs, *z, *pvalZtest, *pvalPerms);
#endif
}

/** \brief First, read all genotypes of all cis SNPs of the given feature.
 *   Second, compute the summary stats for each pair gene-SNP.
 *  \note isNa[i] == true if sample i has no phenotype
 */
void
computeSummaryStatsForOneFeature (
  const map<string, vector<string> >::iterator ftrName2CisSnpNameCoords_it,
  map<string, long int> snpNameCoord2Pos,
  ifstream & genoStream,
  const vector<double> & y,
  const vector<size_t> vIdxIndsToSkip,
  const vector<bool> isNa,
  int nbPermutations,
  FtrStats * pt_iFtrStats,
  size_t * pt_nbFtrsLowPval,
  const int verbose)
{
  vector<string> cisSnpNameCoords = (*ftrName2CisSnpNameCoords_it).second;
  size_t snp_id;
  double thresholdLowPval = 1e-7;
  bool isPvalLow = false;
  
  // for each SNP in cis of the given feature,
  // retrieve its genotypes, compute the OLS summary stats and keep them
  for (snp_id = 0; snp_id < cisSnpNameCoords.size(); snp_id++)
  {
    string snpNameCoord = cisSnpNameCoords[snp_id];
    if (verbose > 0)
      printf ("%s %zu/%zu\n", snpNameCoord.c_str(), snp_id+1,
	      cisSnpNameCoords.size());
    if (snpNameCoord2Pos.find(snpNameCoord) == snpNameCoord2Pos.end())
      continue;
    
    SnpStats iSnpStats;
    vector<string> tokens;
    split (snpNameCoord, '|', tokens);
    iSnpStats.name = tokens[0];
    iSnpStats.coord = atol(tokens[1].c_str());
    
    size_t snpPos = snpNameCoord2Pos[snpNameCoord];
    vector<double> g;
    getGenoValues (pt_iFtrStats->name, snpNameCoord, genoStream, snpPos, g,
		   vIdxIndsToSkip, isNa, &iSnpStats.maf, verbose);
    if (y.size() != g.size())
    {
      cerr << "ERROR: different number of samples for feature "
	   << pt_iFtrStats->name << " (" << y.size() << ") and SNP "
	   << snpNameCoord << " (" << g.size() << ")" << endl;
      exit (1);
    }
    
    ols (pt_iFtrStats->name, snpNameCoord, g, y,
	 &iSnpStats.betahat, &iSnpStats.sebetahat, &iSnpStats.sigmahat,
	 &iSnpStats.pval, &iSnpStats.R2, verbose-1);
    if (verbose > 0)
      printf ("%s %s betahat=%.8f sebetahat=%.8f sigmahat=%.8f P-value=%.8f R2=%.8f\n",
	      pt_iFtrStats->name.c_str(), snpNameCoord.c_str(),
	      iSnpStats.betahat, iSnpStats.sebetahat, iSnpStats.sigmahat,
	      iSnpStats.pval, iSnpStats.R2);
    if (nbPermutations >= 0)
    {
      spearman (pt_iFtrStats->name, snpNameCoord, g, y,
		&iSnpStats.rs, &iSnpStats.rsZscore,
		&iSnpStats.rsPvalZtest, nbPermutations,
		&iSnpStats.rsPvalPerms, verbose-1);
      if (verbose > 0)
	printf ("spearman=%.8f PvalZtest=%.8f permutations=%i PvalPerms=%.8f\n",
		iSnpStats.rs, iSnpStats.rsPvalZtest,
		nbPermutations, iSnpStats.rsPvalPerms);
    }
    else
    {
      iSnpStats.rs = numeric_limits<double>::quiet_NaN();
      iSnpStats.rsZscore = numeric_limits<double>::quiet_NaN();
      iSnpStats.rsPvalZtest = numeric_limits<double>::quiet_NaN();
      iSnpStats.rsPvalPerms = numeric_limits<double>::quiet_NaN();
    }
    if (iSnpStats.pval <= thresholdLowPval)
      isPvalLow = true;
    
    pt_iFtrStats->vSnpStats.push_back (iSnpStats);
  }
  if (isPvalLow)
    ++(*pt_nbFtrsLowPval);
}

/** \brief Compute summary statistics for each pair feature-SNP
 *  and write the results feature by feature in an uncompressed file.
 */
void computeAndWriteSummaryStatsFtrPerFtr (
  map<string, vector<string> > ftrName2CisSnpNameCoords,
  map<string, long int> ftrName2Pos,
  map<string, long int> snpNameCoord2Pos,
  vector<size_t> vIdxIndsToSkip,
  string phenoFile,
  string genoFile,
  string outFile,
  int nbPermutations,
  int verbose)
{
  FtrStats iFtrStats;
  map<string, vector<string> >::iterator ftrName2CisSnpNameCoords_it;
  vector<string> cisSnpNameCoords, tokens;
  string ftrName, snpNameCoord, snpName;
  size_t ftrPos;
  size_t nbFtrs = 0;
  vector<size_t> vCounters = getCounters (ftrName2CisSnpNameCoords.size());
  vector<double> y;
  vector<bool> isNa;
  ifstream phenoStream, genoStream;
  ofstream outStream;
  size_t nbFtrsLowPval = 0;
  
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
  outStream << "ftr snp coord maf n betahat sebetahat sigmahat pval R2";
  outStream << " rs rsZscore rsPvalZtest rsPvalPerms";
  outStream << endl;
  if (verbose > 0)
  {
    cout << "compute summary statistics for each pair feature-SNP..." << endl;
  }
  
  // for each feature
  for (ftrName2CisSnpNameCoords_it = ftrName2CisSnpNameCoords.begin();
       ftrName2CisSnpNameCoords_it != ftrName2CisSnpNameCoords.end();
       ++ftrName2CisSnpNameCoords_it)
  {
    ++nbFtrs;
    if (verbose > 0)
      printCounter (nbFtrs, vCounters);
    ftrName = (*ftrName2CisSnpNameCoords_it).first;
    if (ftrName2Pos.find(ftrName) == ftrName2Pos.end())
      continue;
    if (verbose > 1)
      printf ("%s %zu/%zu\n", ftrName.c_str(), nbFtrs,
	      ftrName2CisSnpNameCoords.size());
    FtrStats_reset (&iFtrStats);
    iFtrStats.name = ftrName;
    
    // retrieve its values
    ftrPos = ftrName2Pos[ftrName];
    getPhenoValues (ftrName, phenoStream, ftrPos, y, vIdxIndsToSkip, isNa,
		    verbose-1);
    
    // loop over SNPs in cis
    computeSummaryStatsForOneFeature (ftrName2CisSnpNameCoords_it,
				      snpNameCoord2Pos,
				      genoStream,
				      y,
				      vIdxIndsToSkip,
				      isNa,
				      nbPermutations,
				      &iFtrStats,
				      &nbFtrsLowPval,
				      verbose-2);
    
    // write the results (one line per SNP)
    FtrStats_write (iFtrStats, y.size(), outStream);
    
    y.clear();
    isNa.clear();
  }
  
  phenoStream.close();
  genoStream.close();
  outStream.close();
  if (verbose > 0)
  {
    cout << "results written in " << outFile << endl;
    cout << "features with at least one SNP with P-value <= 1e-7: " 
	 << nbFtrsLowPval << " / " << ftrName2Pos.size() << endl;
  }
}

int main (int argc, char ** argv)
{
  string linksFile, genoFile, phenoFile, outFile, ftrsFile, snpsFile, indsFile;
  vector<size_t> vIdxIndsToSkip;
  double minMaf = 0.0;
  int nbPermutations = -1;
  int verbose = 1;
  parse_args (argc, argv, &linksFile, &genoFile, &phenoFile, &outFile,
	      &ftrsFile, &snpsFile, &indsFile, &minMaf, &nbPermutations,
	      &verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl;
  }
  
  gsl_rng_env_setup();   // used for permutations (Pval Spearman coef)
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsFile, verbose);
  vector<string> vIndsToSkip = loadOneColumnFile (indsFile, verbose);
  
  map<string, vector<string> > ftrName2CisSnpNameCoords
    = loadLinksFtr2Snps (linksFile, vFtrsToKeep, vSnpsToKeep, verbose);
  
  map<string, long int> ftrName2Pos
    = indexFtrs (phenoFile, ftrName2CisSnpNameCoords, vIndsToSkip,
		 vIdxIndsToSkip, verbose);
  map<string, long int> snpNameCoord2Pos
    = indexSnps (genoFile, vSnpsToKeep, vIdxIndsToSkip, minMaf, verbose);
  
  computeAndWriteSummaryStatsFtrPerFtr (ftrName2CisSnpNameCoords,
					ftrName2Pos,
					snpNameCoord2Pos,
					vIdxIndsToSkip,
					phenoFile,
					genoFile,
					outFile,
					nbPermutations,
					verbose);
  
  gsl_rng_free (r);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime)
	 << ": elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << ")" << endl;
  }
  
  return 0;
}
