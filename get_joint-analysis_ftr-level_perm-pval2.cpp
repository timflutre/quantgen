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
#include <limits>
#include <sstream>
#include <numeric>
using namespace std;

#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_combination.h>

#include "utils.h"

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
  string & anchor,
  size_t & lenCis,
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
	{"anchor", required_argument, 0, 'a'},
	{"cis", required_argument, 0, 0},
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
    c = getopt_long (argc, argv, "hVv:g:p:a:o:P:t:f:s:q",
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
      if (strcmp(long_options[option_index].name, "cis") == 0)
      {
	lenCis = atol (optarg);
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
    case 'a':
      anchor = optarg;
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
  if (anchor.empty())
  {
    fprintf (stderr, "ERROR: SNPs in trans not yet implemented, see --anchor.\n\n");
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

vector< vector<double> >
loadGrid (
  const string & gridFile,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load grid ..." << endl << flush;
  vector< vector<double> > grid;
  ifstream gridStream;
  vector<string> tokens;
  string line;
  openFile (gridFile, gridStream);
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
    cout << "grid size: " << grid.size() << endl;
  return grid;
}

/** \brief Compute the summary statistics of the linear regression.
 *  \note phenotype = mu + genotype * beta + error
 *  \note missing values should have been already filtered out
 */
void
ols (
  const vector<double> & g,
  const vector<double> & y,
  double & betahat,
  double & sebetahat,
  double & sigmahat,
  double & pval,
  double & R2)
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
  printf ("n=%zu ym=%f gm=%f yty=%f gtg=%f gty=%f vg=%f\n",
	  n, ym, gm, yty, gtg, gty, vg);
#endif
  if(vg > 1e-8)
  {
    betahat = (gty - n * gm * ym) / vg;
#ifdef DEBUG
    printf ("betahat=%f\n", betahat);
#endif
    double rss1 = yty - 1/vg * (n*ym*(gtg*ym - gm*gty) - gty*(n*gm*ym - gty));
    if (fabs(betahat) > 1e-8)
      sigmahat = sqrt(rss1 / (n-2));
    else  // case where phenotypes are not variable enough
      sigmahat = sqrt((yty - n * ym * ym) / (n-2));  // sqrt(rss0/(n-2))
#ifdef DEBUG
    printf ("sigmahat=%f\n", sigmahat);
#endif
    sebetahat = sigmahat / sqrt(gtg - n*gm*gm);
#ifdef DEBUG
    printf ("sebetahat=%f\n", sebetahat);
#endif
    double muhat = (ym*gtg - gm*gty) / (gtg - n*gm*gm);
    double mss = 0;
    for(i=0; i<n; ++i)
      mss += pow(muhat + betahat * g[i] - ym, 2);
    pval = gsl_cdf_fdist_Q (mss/pow(sigmahat,2), 1, n-2);
    R2 = mss / (mss + rss1);
  }
  else
  {
#ifdef DEBUG
    cout << "genotypes are not variable enough" << endl << flush;
#endif
    betahat = 0;
    sebetahat = numeric_limits<double>::infinity();
    sigmahat = sqrt((yty - n * ym * ym) / (n-2));  // sqrt(rss0/(n-2))
    pval = 1;
    R2 = 0;
  }
}

struct Snp
{
  string name; // eg. rs7263524
  string chr; // eg. chr21
  size_t coord; // 1-based coordinate
  vector<vector<double> > vvGenos; // genotypes of samples per subgroup
  vector<vector<bool> > vvIsNa; // missing values per subgroup
  vector<double> vMafs; // minor allele frequencies per subgroup
};

struct ResFtrSnp
{
  string snp; // name of the SNP
  vector<size_t> vNs; // sample sizes per subgroup
  vector<double> vBetahats; // MLE of the beta per subgroup
  vector<double> vSebetahats; // standard errors
  vector<double> vSigmahats; // MLEs of the sigmas
  vector<double> vBetaPvals; // P-values of H0:"beta=0" and H1:"beta!=0"
  vector<double> vPves; // proportions of variance explained (R2)
  vector<vector<double> > vvStdSstatsCorr;
  map<string, double> mAbfs;
};

struct Ftr
{
  string name; // eg. ENSG00000182816
  string chr; // eg. chr21
  size_t start; // 1-based coordinate
  size_t end; // idem
  vector<vector<double> > vvPhenos; // phenotypes of samples per subgroup
  vector<vector<bool> > vvIsNa; // missing values per subgroup
  vector<Snp*> vPtCisSnps;
  vector<ResFtrSnp> vResFtrSnps;
  double maxL10TrueAbf;
  double permPval;
  size_t nbPermsSoFar;
};

void
Snp_init (
  Snp & iSnp,
  const string & name,
  const size_t & nbSubgroups,
  const size_t & nbSamplesS1)
{
  iSnp.name = name;
  iSnp.chr.clear();
  iSnp.coord = string::npos;
  iSnp.vvGenos.resize (nbSubgroups);
  iSnp.vvIsNa.resize (nbSubgroups);
  iSnp.vMafs = (vector<double> (nbSubgroups, 0.0));
  iSnp.vvGenos[0] = (vector<double> (nbSamplesS1, 0.0));
  iSnp.vvIsNa[0] = (vector<bool> (nbSamplesS1, false));
}

// assume both features are on the same chromosome
bool Snp_compByCoord (
  const Snp* pt_iSnp1,
  const Snp* pt_iSnp2)
{
  bool res = false;
  if (pt_iSnp1->coord < pt_iSnp2->coord)
    res = true;
  return res;
}

int Snp_isInCis (
  const Snp & iSnp,
  const size_t & ftrStart,
  const size_t & ftrEnd,
  const string & anchor,
  const size_t & lenCis)
{
  int res = -1;
  if (anchor.compare("FSS+FES") == 0)
  {
    if (((ftrStart >= lenCis &&
	 iSnp.coord >= ftrStart - lenCis) ||
	(ftrStart < lenCis &&
	 iSnp.coord >= 0)) &&
	iSnp.coord <= ftrEnd + lenCis)
      res = 0;
    else if (iSnp.coord > ftrEnd + lenCis)
      res = 1;
  }
  else if (anchor.compare("FSS") == 0)
  {
    if (((ftrStart >= lenCis &&
	 iSnp.coord >= ftrStart - lenCis) ||
	(ftrStart < lenCis &&
	 iSnp.coord >= 0)) &&
	iSnp.coord <= ftrStart + lenCis)
      res = 0;
    else if (iSnp.coord > ftrStart + lenCis)
      res = 1;
  }
  return res;
}

void
ResFtrSnp_init (
  ResFtrSnp & iResFtrSnp,
  const string & snpName,
  const size_t & nbSubgroups)
{
  iResFtrSnp.snp = snpName;
  iResFtrSnp.vNs.assign (nbSubgroups, string::npos);
  iResFtrSnp.vBetahats.assign (nbSubgroups, 0.0);
  iResFtrSnp.vSebetahats.assign (nbSubgroups, 0.0);
  iResFtrSnp.vSigmahats.assign (nbSubgroups, 0.0);
  iResFtrSnp.vBetaPvals.assign (nbSubgroups, 0.0);
  iResFtrSnp.vPves.assign (nbSubgroups, 0.0);
}

void
ResFtrSnp_getSstatsOneSbgrp (
  ResFtrSnp & iResFtrSnp,
  const Ftr & iFtr,
  const Snp & iSnp,
  const size_t & s,
  const bool & needQnorm)
{
  vector<double> y, g;
  for (size_t i = 0; i < iFtr.vvPhenos[s].size(); ++i)
    if (! iFtr.vvIsNa[s][i] && ! iSnp.vvIsNa[0][i])
    {
      y.push_back (iFtr.vvPhenos[s][i]);
      g.push_back (iSnp.vvGenos[0][i]);
    }
  if (needQnorm)
    qqnorm (&y[0], y.size());
  
  iResFtrSnp.vNs[s] = y.size();
  
  ols (g, y, iResFtrSnp.vBetahats[s], iResFtrSnp.vSebetahats[s],
       iResFtrSnp.vSigmahats[s], iResFtrSnp.vBetaPvals[s],
       iResFtrSnp.vPves[s]);
}

void
ResFtrSnp_corrSmallSampleSize (
  ResFtrSnp & iResFtrSnp)
{
  for (size_t s = 0; s < iResFtrSnp.vBetahats.size(); ++s)
  {
    double n = iResFtrSnp.vNs[s],
      bhat = iResFtrSnp.vBetahats[s] / iResFtrSnp.vSigmahats[s],
      sebhat = iResFtrSnp.vSebetahats[s] / iResFtrSnp.vSigmahats[s],
      t = gsl_cdf_gaussian_Pinv (gsl_cdf_tdist_P (-fabs(bhat/sebhat),
						  n-2), 1.0),
      sigmahat = 0;
    vector<double> vStdSstatsCorr;
    if (fabs(t) > 1e-8)
    {
      sigmahat = fabs (iResFtrSnp.vBetahats[s]) / (fabs (t) * sebhat);
      bhat = iResFtrSnp.vBetahats[s] / sigmahat;
      sebhat = bhat / t;
    }
    else
    {
      sigmahat = numeric_limits<double>::quiet_NaN();
      bhat = 0;
      sebhat = numeric_limits<double>::infinity();
    }
    vStdSstatsCorr.push_back (bhat);
    vStdSstatsCorr.push_back (sebhat);
    vStdSstatsCorr.push_back (t);
    iResFtrSnp.vvStdSstatsCorr.push_back (vStdSstatsCorr);
  }
}

double
getAbfFromStdSumStats (
  const vector< vector<double> > & stdSumStats,
  const double & phi2,
  const double & oma2)
{
  double lABF_all = 0;
  size_t subgroup = 0;
  double bhat = 0, varbhat = 0, t = 0, bbarhat_num = 0, bbarhat_denom = 0,
    varbbarhat = 0;
  vector<double> lABF_single_s;
#ifdef DEBUG
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
  printf ("bbarhat=%e varbbarhat=%e T2=%e lABF_bbar=%e\n",
	  bbarhat, varbbarhat, T2, lABF_bbar);
#endif
  
  lABF_all = lABF_bbar;
  for (subgroup = 0; subgroup < lABF_single_s.size(); ++subgroup)
    lABF_all += lABF_single_s[subgroup];
#ifdef DEBUG
  printf ("lABF_all=%e\n", lABF_all);
#endif
  
  return lABF_all;
}

void
ResFtrSnp_calcAbfs (
  ResFtrSnp & iResFtrSnp,
  const vector< vector<double> > & grid)
{
  ResFtrSnp_corrSmallSampleSize (iResFtrSnp);
  
  vector<double> l10AbfsConst, l10AbfsFix, l10AbfsMaxh;
  for (size_t gridId = 0; gridId < grid.size(); ++gridId)
  {
    l10AbfsConst.push_back (getAbfFromStdSumStats (
			      iResFtrSnp.vvStdSstatsCorr,
			      grid[gridId][0],
			      grid[gridId][1]));
    l10AbfsFix.push_back (getAbfFromStdSumStats (
			    iResFtrSnp.vvStdSstatsCorr,
			    0,
			    grid[gridId][0]
			    + grid[gridId][1]));
    l10AbfsMaxh.push_back (getAbfFromStdSumStats (
			     iResFtrSnp.vvStdSstatsCorr,
			     grid[gridId][0]
			     + grid[gridId][1],
			     0));
  }
  
  vector<double> vWeights (grid.size(), 1.0 / (double) grid.size());
  iResFtrSnp.mAbfs.insert (make_pair ("const", log10_weighted_sum (
					&(l10AbfsConst[0]),
					&(vWeights[0]),
					grid.size())));
  iResFtrSnp.mAbfs.insert (make_pair ("fix", log10_weighted_sum (
					&(l10AbfsFix[0]),
					&(vWeights[0]),
					grid.size())));
  iResFtrSnp.mAbfs.insert (make_pair ("maxh", log10_weighted_sum (
					&(l10AbfsMaxh[0]),
					&(vWeights[0]),
					grid.size())));
}

void
Ftr_init (
  Ftr & iFtr,
  const string & name,
  const size_t & nbSubgroups,
  const size_t & nbSamplesS1)
{
  iFtr.name = name;
  iFtr.chr.clear();
  iFtr.start = string::npos;
  iFtr.end = string::npos;
  iFtr.vvPhenos.resize (nbSubgroups);
  iFtr.vvIsNa.resize (nbSubgroups);
  iFtr.vvPhenos[0] = (vector<double> (nbSamplesS1, 0.0));
  iFtr.vvIsNa[0] = (vector<bool> (nbSamplesS1, false));
}

// assume both features are on the same chromosome
bool Ftr_compByCoord (
  const Ftr* pt_iFtr1,
  const Ftr* pt_iFtr2)
{
  bool res = false;
  if ((pt_iFtr1->start < pt_iFtr2->start) ||
      (pt_iFtr1->start == pt_iFtr2->start && 
       pt_iFtr1->end < pt_iFtr2->end ))
      res = true;
  return res;
}

void
Ftr_getCisSnps (
  Ftr & iFtr,
  const map<string, vector<Snp*> > & mChr2VecPtSnps,
  const string & anchor,
  const size_t & lenCis)
{
  map<string, vector<Snp*> >::const_iterator itVecPtSnps =
    mChr2VecPtSnps.find(iFtr.chr);
  
  for (size_t snpId = 0; snpId < itVecPtSnps->second.size(); ++snpId)
  {
    Snp * ptSnp = (itVecPtSnps->second)[snpId];
    int inCis = Snp_isInCis (*ptSnp, iFtr.start, iFtr.end,
			     anchor, lenCis);
    if (inCis == 1)
      break;
    else if (inCis == -1)
      continue;
    iFtr.vPtCisSnps.push_back ((itVecPtSnps->second)[snpId]);
  }
}

void
Ftr_inferAssos (
  Ftr & iFtr,
  const vector<vector<double> > & grid,
  const string & whichBf,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const bool & needQnorm,
  const int & verbose)
{
  size_t nbSubgroups = iFtr.vvPhenos.size();
  Snp * ptSnp;
  for (size_t snpId = 0; snpId < iFtr.vPtCisSnps.size(); ++snpId)
  {
    ptSnp = iFtr.vPtCisSnps[snpId];
    ResFtrSnp iResFtrSnp;
    ResFtrSnp_init (iResFtrSnp, ptSnp->name, nbSubgroups);
    for (size_t s = 0; s < nbSubgroups; ++s)
      ResFtrSnp_getSstatsOneSbgrp (iResFtrSnp, iFtr, *ptSnp, s, needQnorm);
    ResFtrSnp_calcAbfs (iResFtrSnp, grid);
    iFtr.vResFtrSnps.push_back (iResFtrSnp);
  }
  
//  Ftr_doPerms ();
}

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
    if (tokens[0].compare("chr") != 0)
    {
      cerr << "ERROR: file " << vGenoPaths[s] << " requires a header" << endl;
      exit (1);
    }
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

void
loadPhenos (
  const vector<string> & vPhenoPaths,
  const vector<string> & vFtrsToKeep,
  map<string, Ftr> & mFtrs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load phenotypes ..." << endl << flush;
  
  ifstream phenoStream;
  string line;
  vector<string> tokens;
  size_t nbSamples, nbLines;
  
  for (size_t s = 0; s < vPhenoPaths.size(); ++s)
  {
    openFile (vPhenoPaths[s], phenoStream);
    getline (phenoStream, line); // header
    split (line, " \t", tokens);
    nbSamples = tokens.size() - 1;
    nbLines = 1;
		
    while (true)
    {
      getline (phenoStream, line);
      if (line.empty())
	break;
      ++nbLines;
      split (line, " \t", tokens);
      if (! vFtrsToKeep.empty()
	  && find (vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0])
	  == vFtrsToKeep.end())
	continue;
      if (tokens.size() != nbSamples + 1)
      {
	cerr << "ERROR: not enough columns on line " << nbLines
	     << " of file " << vPhenoPaths[s] << endl;
	exit (1);
      }
			
      if (mFtrs.find(tokens[0]) == mFtrs.end())
      {
	Ftr iFtr;
	Ftr_init (iFtr, tokens[0], vPhenoPaths.size(), nbSamples);
	for (size_t i = 1; i < tokens.size(); ++i)
	{
	  if (tokens[i].compare("NA") == 0)
	    iFtr.vvIsNa[s][i-1] = true;
	  iFtr.vvPhenos[s][i-1] = atof (tokens[i].c_str());
	}
	mFtrs.insert (make_pair (tokens[0], iFtr));
      }
      else
      {
	mFtrs[tokens[0]].vvPhenos[s] = (vector<double> (nbSamples, 0.0));
	mFtrs[tokens[0]].vvIsNa[s] = (vector<bool> (nbSamples, false));
	for (size_t i = 1; i < tokens.size() ; ++i)
	{
	  if (tokens[i].compare("NA") == 0)
	    mFtrs[tokens[0]].vvIsNa[s][i-1] = true;
	  mFtrs[tokens[0]].vvPhenos[s][i-1] = atof (tokens[i].c_str());
	}
      }
    }
    
    phenoStream.close();
  }
  
  if (mFtrs.size() == 0)
  {
    cerr << "ERROR: no feature to analyze" << endl;
    exit (1);
  }
  if (verbose > 0)
    cout << "nb of features: " << mFtrs.size() << endl;
/*	map<string, Ftr>::iterator it = mFtrs.begin();
	while (it != mFtrs.end())
	{
	for (size_t s = 0; s < it->second.vvPhenos.size(); ++s)
	for (size_t i = 0; i < it->second.vvPhenos[s].size(); ++i)
	cout << it->first << " " << (s+1) << " " << (i+1) << " "
	<< it->second.vvPhenos[s][i] << endl;
	++it;
	}*/
}

void
loadFtrInfo (
  const string & ftrCoordsFile,
  map<string, Ftr> & mFtrs,
  map<string, vector<Ftr*> > & mChr2VecPtFtrs,
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
    if (mFtrs.find(tokens[3]) == mFtrs.end())
      continue;
    ++countFtrs;
    mFtrs[tokens[3]].chr = tokens[0];
    mFtrs[tokens[3]].start = atol (tokens[1].c_str()) + 1;
    mFtrs[tokens[3]].end = atol (tokens[2].c_str());
    
    if (mChr2VecPtFtrs.find(tokens[0]) == mChr2VecPtFtrs.end())
      mChr2VecPtFtrs.insert (make_pair (tokens[0],
					vector<Ftr*> ()));
    mChr2VecPtFtrs[tokens[0]].push_back (&(mFtrs[tokens[3]]));
  }
  
  ftrCoordsStream.close();
  
  if (countFtrs < mFtrs.size())
  {
    cerr << "ERROR: " << mFtrs.size() - countFtrs
	 << " feature coordinates are missing" << endl;
    exit (1);
  }
  
  // sort the features per chr
  for (map<string, vector<Ftr*> >::iterator it = mChr2VecPtFtrs.begin();
       it != mChr2VecPtFtrs.end(); ++it)
    sort (it->second.begin(), it->second.end(), Ftr_compByCoord);
}

void
loadGenosAndSnpInfo (
  const vector<string> & vGenoPaths,
  const vector<string> & vSnpsToKeep,
  map<string, Snp> & mSnps,
  map<string, vector<Snp*> > & mChr2VecPtSnps,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load genotypes and SNP coordinates ..." << endl << flush;
  
  ifstream genoStream;
  string line;
  vector<string> tokens;
  size_t nbSamples, nbLines;
  double maf, AA, AB, BB;
  
  for (size_t s = 0; s < vGenoPaths.size(); ++s)
  {
    openFile (vGenoPaths[s], genoStream);
    getline (genoStream, line); // header
    split (line, " \t", tokens);
    nbSamples = (size_t) (tokens.size() - 5) / 3;
    nbLines = 1;
    
    while (true)
    {
      getline (genoStream, line);
      if (line.empty())
	break;
      ++nbLines;
      split (line, " \t", tokens);
      if (! vSnpsToKeep.empty()
	  && find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[1])
	  == vSnpsToKeep.end())
	continue;
      if (tokens.size() != (size_t) (3 * nbSamples + 5))
      {
	cerr << "ERROR: not enough columns on line " << nbLines
	     << " of file " << vGenoPaths[s] << endl;
	exit (1);
      }
			
      if (mSnps.find(tokens[1]) == mSnps.end())
      {
	Snp iSnp;
	Snp_init (iSnp, tokens[1], vGenoPaths.size(), nbSamples);
	maf = 0;
	for (size_t i = 0; i < nbSamples; ++i)
	{
	  AA = atof(tokens[5+3*i].c_str());
	  AB = atof(tokens[5+3*i+1].c_str());
	  BB = atof(tokens[5+3*i+2].c_str());
	  if (AA == 0 && AB == 0 && BB == 0)
	    iSnp.vvIsNa[s][i] = true;
	  else
	  {
	    iSnp.vvGenos[s][i] = 0 * AA + 1 * AB + 2 * BB;
	    maf += iSnp.vvGenos[s][i];
	  }
	}
	maf /= 2 * (nbSamples,
		    - count (iSnp.vvIsNa[s].begin(),
			     iSnp.vvIsNa[s].end(),
			     true));
	iSnp.vMafs[s] = maf <= 0.5 ? maf : (1 - maf);
	iSnp.chr = tokens[0];
	iSnp.coord = atol (tokens[2].c_str());
	mSnps.insert (make_pair (tokens[1], iSnp));
	
	if (mChr2VecPtSnps.find(tokens[0]) == mChr2VecPtSnps.end())
	  mChr2VecPtSnps.insert (make_pair (tokens[0],
					    vector<Snp*> ()));
	mChr2VecPtSnps[tokens[0]].push_back (&(mSnps[tokens[1]]));
      }
    }
    
    genoStream.close();
  }
  
  for (map<string, Snp>::iterator it = mSnps.begin();
       it != mSnps.end(); ++it)
    if (it->second.vvGenos.size() == 0)
    {
      cerr << "ERROR: SNP " << it->first << " has no genotype" << endl;
      exit (1);
    }
  
  // sort the SNPs per chr
  for (map<string, vector<Snp*> >::iterator it = mChr2VecPtSnps.begin();
       it != mChr2VecPtSnps.end(); ++it)
       sort (it->second.begin(), it->second.end(), Snp_compByCoord);
}

void
inferAssos (
  map<string, Ftr> & mFtrs,
  const map<string, vector<Ftr*> > & mChr2VecPtFtrs,
  const map<string, Snp> & mSnps,
  const map<string, vector<Snp*> > & mChr2VecPtSnps,
  const vector<vector<double> > & grid,
  const string & whichBf,
  const string & anchor,
  const size_t & lenCis,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const bool & needQnorm,
  const int & verbose)
{
  if (verbose > 0)
  {
    cout << "look for association between each pair feature-SNP ..." << endl
	 << "anchor=" << anchor << " lenCis=" << lenCis << endl << flush;
    if (nbPerms > 0)
    {
      cout << "permutation"<< (nbPerms > 1 ? "s=" : "=") << nbPerms
	   << ", seed=" << seed
	   << ", trick=" << trick
	   << endl << flush;
    }
  }
  
  size_t nbAnalyzedPairs = 0;
  
  for (map<string, Ftr>::iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    Ftr_getCisSnps (itF->second, mChr2VecPtSnps, anchor, lenCis);
    if (itF->second.vPtCisSnps.size() > 0)
    {
      if (verbose > 1)
	cout << itF->second.name << ": " << itF->second.vPtCisSnps.size()
	     << " SNPs in cis" << endl << flush;
      Ftr_inferAssos (itF->second, grid, whichBf, nbPerms, seed, trick,
		      needQnorm, verbose-1);
      nbAnalyzedPairs += itF->second.vResFtrSnps.size();
    }
  }
  
  if (verbose > 0)
    cout << "nb of analyzed feature-SNP pairs: " << nbAnalyzedPairs << endl;
}

bool isNonZero (size_t i) { return (i != 0); };

void
writeRes (
  const string & outPrefix,
  const map<string, Ftr> & mFtrs,
  const map<string, Snp> & mSnps,
  const size_t & nbSubgroups,
  const int & verbose)
{
  // 1 file per subgrp with summary stats
  for (size_t s = 0; s < nbSubgroups; ++s)
  {
    stringstream ss;
    ss << outPrefix << "_sumstats_s" << (s+1) << ".txt";
    ofstream outStream;
    openFile (ss.str(), outStream);
    
    outStream << "ftr snp betahat sebetahat sigmahat betaPval pve";
    outStream << endl;
    
    for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
	 itF != mFtrs.end(); ++itF)
    {
      const Ftr * ptF = &(itF->second);
      for (vector<ResFtrSnp>::const_iterator itP = ptF->vResFtrSnps.begin();
	   itP != ptF->vResFtrSnps.end(); ++itP)
	outStream << ptF->name
		  << " " << itP->snp
		  << " " << itP->vBetahats[s]
		  << " " << itP->vSebetahats[s]
		  << " " << itP->vSigmahats[s]
		  << " " << itP->vBetaPvals[s]
		  << " " << itP->vPves[s]
		  << endl;
    }
    
    outStream.close();
  }
  
  // 1 file with ABFs
  stringstream ss;
  ss << outPrefix << "_abfs.txt";
  ofstream outStream;
  openFile (ss.str(), outStream);
  outStream << "ftr snp nb.samples nb.subgroups abf.const abf.fix abf.maxh";
  outStream << endl;
  for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    const Ftr * ptF = &(itF->second);
    for (vector<ResFtrSnp>::const_iterator itP = ptF->vResFtrSnps.begin();
	 itP != ptF->vResFtrSnps.end(); ++itP)
      outStream << ptF->name
		<< " " << itP->snp
		<< " " << accumulate (itP->vNs.begin(), itP->vNs.end(), 0)
		<< " " << count_if (itP->vNs.begin(), itP->vNs.end(), isNonZero)
		<< " " << itP->mAbfs.find("const")->second
		<< " " << itP->mAbfs.find("fix")->second
		<< " " << itP->mAbfs.find("maxh")->second
		<< endl;
  }
  outStream.close();
}

/*
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
*/
void
run (
  const string & genoPathsFile,
  const string & phenoPathsFile,
  const string & ftrCoordsFile,
  const string & gridFile,
  const string & whichBf,
  const string & anchor,
  const size_t & lenCis,
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
  if (vGenoPaths.size() > 1)
  {
    cerr << "ERROR: current version can't handle different genotype files" << endl;
    exit (1);
  }
  vector<string> vPhenoPaths = loadOneColumnFile (phenoPathsFile, verbose);
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsToKeepFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsToKeepFile, verbose);
  vector<vector<double> > grid = loadGrid (gridFile, verbose);
  
  vector<string> vSamples;
  vector<vector<size_t> > vvSampleIdxs;
  loadSamples (vGenoPaths, vPhenoPaths, vSamples, vvSampleIdxs, verbose);
  
  map<string, Ftr> mFtrs;
  map<string, vector<Ftr*> > mChr2VecPtFtrs;
  loadPhenos (vPhenoPaths, vFtrsToKeep, mFtrs, verbose);
  loadFtrInfo (ftrCoordsFile, mFtrs, mChr2VecPtFtrs, verbose);
  
  map<string, Snp> mSnps;
  map<string, vector<Snp*> > mChr2VecPtSnps;
  loadGenosAndSnpInfo (vGenoPaths, vSnpsToKeep, mSnps, mChr2VecPtSnps,
		       verbose);
  
  inferAssos (mFtrs, mChr2VecPtFtrs, mSnps, mChr2VecPtSnps, grid, whichBf,
	      anchor, lenCis, nbPerms, seed, trick, needQnorm, verbose);
  
  writeRes (outPrefix, mFtrs, mSnps, vPhenoPaths.size(), verbose);
}

int main (int argc, char ** argv)
{
  int verbose = 1, trick = 0;
  string genoPathsFile, phenoPathsFile, ftrCoordsFile, gridFile,
    whichBf = "abf.const", anchor = "FSS+FES", outPrefix,
    ftrsToKeepFile, snpsToKeepFile;
  size_t lenCis = 100000, nbPerms = 10000, seed = string::npos;
  bool needQnorm = false;
  
  parseArgs (argc, argv, genoPathsFile, phenoPathsFile, ftrCoordsFile,
	     gridFile, whichBf, anchor, lenCis, outPrefix, nbPerms, seed,
	     trick, ftrsToKeepFile, snpsToKeepFile, needQnorm, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl;
  }
  
  run (genoPathsFile, phenoPathsFile, ftrCoordsFile, gridFile, whichBf,
       anchor, lenCis, outPrefix, nbPerms, seed, trick, ftrsToKeepFile,
       snpsToKeepFile, needQnorm, verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime)
	 << ": elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << ")" << endl;
  }
  
  return EXIT_SUCCESS;
}
