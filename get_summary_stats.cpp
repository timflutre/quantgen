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
 *  g++ -Wall -fopenmp -O3 -DGET_SUMSTATS_MAIN get_summary_stats.cpp utils.cpp -lgsl -lgslcblas -o get_summary_stats
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

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#include "utils.h"

struct SnpStats
{
  string name; // eg. rs2205177
  string chr; // eg. chr21
  size_t coord; // 1-based coordinate
  vector<double> g_init; // genotypes of samples
  vector<bool> vIsNa; // missing values
  double maf; // minor allele frequency
  size_t n; // sample size
  double betahat; // MLE of beta
  double sebetahat; // standard error
  double sigmahat; // MLE of sigma
  double betaPval; // P-value of the test where H0:"beta=0" and H1:"beta!=0"
  double R2; // coef of determination (proportion of variance explained)
  double rs; // Spearman rank correlation coefficient
  double rsZscore; // using Fisher transformation
  double rsPval; // using rsZscore
};

struct FtrStats
{
  string name; // eg. ENSG00000182816
  string chr; // eg. chr21
  size_t start; // 1-based coordinate
  size_t end; // idem
  vector<double> y_init; // phenotypes of samples
  vector<bool> vIsNa; // missing values
  vector <SnpStats> vSnpStats;
  double betaPermPval; // permutation P-value based on beta P-value
  double rsPermPval; // permutation P-value based on Spearman-derived P-value
};

void FtrStats_reset (FtrStats & iFtrStats)
{
  iFtrStats.name.clear();
  iFtrStats.chr.clear();
  iFtrStats.start = string::npos;
  iFtrStats.end = string::npos;
  iFtrStats.y_init.clear();
  iFtrStats.vIsNa.clear();
  iFtrStats.vSnpStats.clear();
  iFtrStats.betaPermPval = numeric_limits<double>::quiet_NaN();
  iFtrStats.rsPermPval = numeric_limits<double>::quiet_NaN();
}

void
FtrStats_init (
  FtrStats & iFtrStats,
  ifstream & phenoStream,
  ifstream & ftrCoordsStream,
  const string chrToKeep,
  const vector<bool> & vIdxSamplesToSkip,
  const vector<string> vFtrsToKeep,
  const gsl_permutation * perm,
  const int verbose)
{
  string linePheno, lineFtrCoords, lineLinks, tok;
  vector<string> tokensPheno, tokensFtrCoords, tokensLinks;
  size_t nbSamples = count (vIdxSamplesToSkip.begin(),
			    vIdxSamplesToSkip.end(),
			    false);
  
  FtrStats_reset (iFtrStats);
  
  while (true)
  {
    getline (phenoStream, linePheno);
    if (linePheno.empty())
      break;
    
    if (linePheno.find('\t') != string::npos)
      split (linePheno, '\t', tokensPheno);
    else
      split (linePheno, ' ', tokensPheno);
    if (tokensPheno.size()-1 != vIdxSamplesToSkip.size())
    {
      cerr << "ERROR: different number of samples for feature "
	   << tokensPheno[0] << " (" << tokensPheno.size()-1
	   << " vs " << vIdxSamplesToSkip.size() << ")" << endl;
      exit (1);
    }
    if (! vFtrsToKeep.empty()
	&& find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokensPheno[0])
	== vFtrsToKeep.end())
      break;
    
    while (true)
    {
      getline (ftrCoordsStream, lineFtrCoords);
      if (lineFtrCoords.empty())
	break;
      
      if (lineFtrCoords.find('\t') != string::npos)
	split (lineFtrCoords, '\t', tokensFtrCoords);
      else
	split (lineFtrCoords, ' ', tokensFtrCoords);
      if (tokensPheno[0].compare(tokensFtrCoords[3]) == 0)
	break;
    }
    if (lineFtrCoords.empty())
    {
      cerr << "ERROR: can't find coordinates of feature "
	   << tokensPheno[0] << endl;
      exit (1);
    }
    if (! chrToKeep.empty() && chrToKeep.compare(tokensFtrCoords[0]) != 0)
      break;
    
    iFtrStats.name = tokensPheno[0];
    iFtrStats.chr = tokensFtrCoords[0];
    iFtrStats.start = atol(tokensFtrCoords[1].c_str()) + 1;
    iFtrStats.end = atol(tokensFtrCoords[2].c_str());
    
    iFtrStats.vIsNa.assign (nbSamples, false);
    iFtrStats.y_init.assign (nbSamples, 0);
    size_t j = 0;
    for (size_t colIdx = 1; colIdx < tokensPheno.size(); ++colIdx)
    {
      if(vIdxSamplesToSkip[colIdx-1])
	continue;
      if (perm != 0 &&
	  vIdxSamplesToSkip[gsl_permutation_get (perm, colIdx-1) + 1])
	continue;
      if (perm == 0)
	tok = tokensPheno[colIdx];
      else
	tok = tokensPheno[gsl_permutation_get (perm, colIdx-1) + 1];
      if (tok.compare("NA") == 0)
	iFtrStats.vIsNa[j] = true;
      else
	iFtrStats.y_init[j] = atof(tok.c_str());
      ++j;
    }
    
    break;
  }
}

void
FtrStats_getCisSnps (
  FtrStats & iFtrStats,
  ifstream & linksStream,
  const map<string, streampos> & mSnpNameCoord2Pos)
{
  string line;
  streampos linksPos;
  vector<string> tokens, tokens2;
  bool hasFtrBeenSeen = false;
  
  iFtrStats.vSnpStats.clear();
  
  while (true)
  {
    linksPos = linksStream.tellg();
    getline (linksStream, line);
    if (line.empty()) // end of file
      break;
    if (line.find('\t') != string::npos)
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    if (tokens.size() != 2)
    {
      cerr << line << endl;
      cerr << "ERROR: format of links file should be"
	   << " feature<space/tab>snp|coord" << endl;
      exit (1);
    }
    if (iFtrStats.name.compare(tokens[0]) != 0)
    {
      if (! hasFtrBeenSeen)
	continue;
      else
      {
	linksStream.seekg (linksPos);
	break;
      }
    }
    hasFtrBeenSeen = true;
    split (tokens[1], '|', tokens2);
    if (mSnpNameCoord2Pos.find(tokens2[0]) == mSnpNameCoord2Pos.end())
      continue; // SNP to skip (see indexSnps)
    SnpStats iSnpStats;
    iSnpStats.name = tokens2[0];
    iFtrStats.vSnpStats.push_back (iSnpStats);
  }
}

void
FtrStats_getCisSnps (
  FtrStats & iFtrStats,
  ifstream & genoStream,
  const map<string, streampos> & mSnpNameCoord2Pos,
  const string & anchor,
  const size_t & lenCis)
{
  string line;
  vector<string> tokens;
  
  iFtrStats.vSnpStats.clear();
  genoStream.clear ();
  genoStream.seekg (0, ios::beg);
  
  while (true)
  {
    getline (genoStream, line);
    if (line.empty()) // end of file
      break;
    if (line.find('\t') != string::npos)
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    if (mSnpNameCoord2Pos.find(tokens[1]) == mSnpNameCoord2Pos.end())
      continue; // SNP to skip (see indexSnps)
    if (tokens[0].compare(iFtrStats.chr) != 0)
      continue; // not on same chr
    if (anchor.compare("TSS") == 0)
    {
      if ((size_t) atol(tokens[2].c_str()) >= iFtrStats.start - lenCis
	  && (size_t) atol(tokens[2].c_str()) <= iFtrStats.start + lenCis)
      {
	SnpStats iSnpStats;
	iSnpStats.name = tokens[1];
	iFtrStats.vSnpStats.push_back (iSnpStats);
      }
    }
    else if (anchor.compare("TSS+TES") == 0)
    {
      if ((size_t) atol(tokens[2].c_str()) >= iFtrStats.start - lenCis
	  && (size_t) atol(tokens[2].c_str()) <= iFtrStats.end + lenCis)
      {
	SnpStats iSnpStats;
	iSnpStats.name = tokens[1];
	iFtrStats.vSnpStats.push_back (iSnpStats);
      }
    }
  }
}

void FtrStats_write (FtrStats iFtrStats, ostream & outStream,
		     const size_t nbPermutations, const bool calcSpearman)
{
  SnpStats iSnpStats;
  for (size_t snp_id = 0; snp_id < iFtrStats.vSnpStats.size(); snp_id++)
  {
    iSnpStats = iFtrStats.vSnpStats[snp_id];
    if (iSnpStats.name.empty())
      continue;
    outStream << iFtrStats.name
	      << " " << iFtrStats.chr
	      << " " << iFtrStats.start
	      << " " << iFtrStats.end
	      << " " << iSnpStats.name
	      << " " << iSnpStats.chr
	      << " " << iSnpStats.coord
	      << " " << iSnpStats.maf
	      << " " << iSnpStats.n;
    if (! calcSpearman)
    {
      outStream << " " << iSnpStats.betahat
		<< " " << iSnpStats.sebetahat
		<< " " << iSnpStats.sigmahat
		<< " " << iSnpStats.betaPval
		<< " " << iSnpStats.R2;
      if (nbPermutations > 0)
	outStream << " " << iFtrStats.betaPermPval;
    }
    else
    {
      outStream << " " << iSnpStats.rs
		<< " " << iSnpStats.rsZscore
		<< " " << iSnpStats.rsPval;
      if (nbPermutations > 0)
	outStream << " " << iFtrStats.rsPermPval;
    }
    outStream << endl;
  }
}

void
SnpStats_init (
  SnpStats & iSnpStats,
  ifstream & genoStream,
  map<string, streampos> & mSnpNameCoord2Pos,
  const vector<bool> & vIdxSamplesToSkip,
  const int verbose)
{
  size_t i, j;
  double AA, AB, BB, maf;
  string line;
  vector<string> tokens;
  size_t nbSamples = count (vIdxSamplesToSkip.begin(),
			    vIdxSamplesToSkip.end(),
			    false);
  
  genoStream.clear ();
  genoStream.seekg (mSnpNameCoord2Pos[iSnpStats.name]);
  getline (genoStream, line);
  if (line.empty())
  {
    cerr << "ERROR: empty genotype line for SNP " << iSnpStats.name << endl;
    exit (1);
  }
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  if (iSnpStats.name.compare(tokens[1]) != 0)
  {
    cerr << "ERROR: wrong genotype line for SNP " << iSnpStats.name << endl;
    exit (1);
  }
  
  iSnpStats.chr = tokens[0];
  iSnpStats.coord = atol(tokens[2].c_str());
  
  maf = 0;
  iSnpStats.vIsNa.assign (nbSamples, false);
  iSnpStats.g_init.assign (nbSamples, 0);
  j = 0;
  for (i = 0; i < vIdxSamplesToSkip.size(); ++i)
  {
    if (vIdxSamplesToSkip[i])
      continue;
    AA = atof(tokens[5+3*i].c_str());
    AB = atof(tokens[5+3*i+1].c_str());
    BB = atof(tokens[5+3*i+2].c_str());
    if (AA == 0 && AB == 0 && BB == 0)
      iSnpStats.vIsNa[j] = true;
    else
    {
      iSnpStats.g_init[j] = 0 * AA + 1 * AB + 2 * BB;
      maf += iSnpStats.g_init[j];
    }
    ++j;
  }
  maf /= 2 * (vIdxSamplesToSkip.size()
	      - count (vIdxSamplesToSkip.begin(),
		       vIdxSamplesToSkip.end(),
		       true)
	      - count (iSnpStats.vIsNa.begin(),
		       iSnpStats.vIsNa.end(),
		       true));
  iSnpStats.maf = maf <= 0.5 ? maf : (1 - maf);
}

/** \brief Return the minor allele frequency.
 *  \note The input comes from a line in the IMPUTE format that was splitted.
 */
double
getMaf (
  const vector<string> & tokens,
  const vector<bool> & vIdxSamplesToSkip)
{
  size_t i, nbNAs = 0;
  double AA, AB, BB, maf = 0;
  for (i = 0; i < vIdxSamplesToSkip.size(); ++i)
  {
    if(vIdxSamplesToSkip[i])
      continue;
    AA = atof(tokens[5+3*i].c_str());
    AB = atof(tokens[5+3*i+1].c_str());
    BB = atof(tokens[5+3*i+2].c_str());
    if (AA == 0 && AB == 0 && BB == 0) // missing value
      ++nbNAs;
    maf += 0 * AA + 1 * AB + 2 * BB;
  }
  maf /= 2 * (vIdxSamplesToSkip.size()
	      - count (vIdxSamplesToSkip.begin(),
		       vIdxSamplesToSkip.end(),
		       true)
	      - nbNAs);
  return maf <= 0.5 ? maf : (1 - maf);
}

/** \brief Index the file with the genotype values (IMPUTE format).
 */
void
indexSnps (
  ifstream & genoStream,
  const vector<string> & vSnpsToKeep,
  const string chrToKeep,
  const vector<string> & vSamplesToSkip,
  vector<string> & vSamples,
  vector<bool> & vIdxSamplesToSkip,
  const double minMaf,
  map<string, streampos> & mSnpNameCoord2Pos,
  const int verbose)
{
  string line;
  vector<string> tokens;
  streampos snpPos;
  size_t totNbSamples, nbSamples;
  
  if (verbose > 0)
    cout << "index SNPs in genotype file ..." << endl;
  
  // read header line and record sample names with their column indices
  genoStream.clear ();
  genoStream.seekg (0, ios::beg);
  getline (genoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  totNbSamples = tokens.size() - 5;
  if (verbose > 0)
    printf ("total number of samples: %zu\n", totNbSamples);
  vSamples.assign (totNbSamples, "");
  vIdxSamplesToSkip.assign (totNbSamples, false);
  for(size_t colIdx = 5; colIdx < tokens.size(); ++colIdx)
  {
    vSamples[colIdx-5] = tokens[colIdx];
    if(! vSamplesToSkip.empty() &&
       find(vSamplesToSkip.begin(), vSamplesToSkip.end(), tokens[colIdx])
       != vSamplesToSkip.end())
      vIdxSamplesToSkip[colIdx-5] = true;
  }
  nbSamples = count (vIdxSamplesToSkip.begin(),
		     vIdxSamplesToSkip.end(),
		     false);
  if (verbose > 0 && nbSamples != totNbSamples)
    printf ("samples to discard: %zu\n", totNbSamples - nbSamples);
  
  while (genoStream.good())
  {
    snpPos = genoStream.tellg();
    getline (genoStream, line);
    if (line.empty())
      break;
    if (line.find('\t') != string::npos)
      split (line, '\t', tokens);
    else
      split (line, ' ', tokens);
    
    if (tokens.size() != 5+3*totNbSamples)
    {
      cerr << line << endl;
      cerr << "ERROR: SNP lines in genotype file should have "
	   << 5+3*totNbSamples << " columns" << endl;
      exit (1);
    }
    
    if (! vSnpsToKeep.empty()
	&& find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[1])
	== vSnpsToKeep.end())
      continue;
    if (! chrToKeep.empty() && chrToKeep.compare(tokens[0]) != 0)
      continue;
    if (minMaf > 0 && getMaf (tokens, vIdxSamplesToSkip) < minMaf)
      continue;
    
    mSnpNameCoord2Pos.insert (make_pair(tokens[1], snpPos));
  }
  
  if (verbose > 0)
    cout << "nb of indexed SNPs: " << mSnpNameCoord2Pos.size() << endl;
}

/** \brief Quantile-normalize an input vector to a standard normal.
 *  \note Missing values should be removed beforehand.
*/
void qnorm (vector<double> & vData)
{
  size_t i, n = vData.size();
  size_t * p;
  double q;
  
  p = (size_t*) calloc (n, sizeof(size_t));
  if (p == NULL)
  {
    cerr << "ERROR: can't allocate memory for p" << endl;
    exit (1);
  }
  
  gsl_sort_index (p, &vData[0], 1, n);
  
  for (i=0; i<n; ++i)
  {
    q = ((double)(i + 0.5)) / ((double)(n));
    vData[p[i]] = gsl_cdf_ugaussian_Pinv (q);
  }
  
  free (p);
}

/** \brief Compute the summary statistics of the linear regression.
 *  \note phenotype = mu + genotype * beta + error
 *  \note missing values should have been already filtered out
 */
void
ols (
  const string yName,
  const string xName,
  const vector<double> & g,
  const vector<double> & y,
  double * betahat,
  double * sebetahat,
  double * sigmahat,
  double * pval,
  double * R2,
  int verbose)
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

/** \brief Compute the permutation P-value at the feature level using, as SNP-
 *  level test statistic, the  P-values on beta=0.
 *
 *  T_i,j: test statistic for SNP i and gene j
 *  T_min,j = min_i (T_i,j)
 *  T~_i,j,k: test statistic for SNP i, gene j and permutation k
 *  T~_min,j,k = min_i (T~_i,j,k)
 *  feature-level Pval = (#{k: T~_min,j,k <= T_min,j} + 1) / (#permutations + 1)
 */
void
computePermutationPvaluesAtFeatureLevel (
  FtrStats & iFtrStats,
  const size_t & nbPermutations,
  const bool & needQnorm,
  const double & minBetaPval,
  const bool & calcSpearman,
  const double & minRsPval,
  const bool & trickPerm,
  gsl_rng * rng,
  const int & verbose)
{
  size_t i, p, perm_id, snp_id;
  double betahat, sebetahat, sigmahat, R2;
  gsl_permutation * perm; // for the indices of the phenotypes vector
  vector<double> resPerm;
  
  if (verbose > 0)
  {
    printf ("perform %zu permutations on the phenotypes (minBetaPval=%f) ...\n",
	    nbPermutations, minBetaPval);
    fflush (stdout);
  }
  
  perm = gsl_permutation_calloc (iFtrStats.y_init.size());
  if (perm == 0)
  {
    cerr << "ERROR: can't allocate memory for the permutation" << endl;
    exit (1);
  }
  
  if (! calcSpearman)
    iFtrStats.betaPermPval = 1;
  else
    iFtrStats.rsPermPval = 1;
  
  for(perm_id=0; perm_id<nbPermutations; ++perm_id)
  {
    double minBetaPvalPerm = 1, betaPvalPerm, minRsPvalPerm = 1, rsPvalPerm,
      rsZscore;
    gsl_ran_shuffle (rng, perm->data, perm->size, sizeof(size_t));
    
    for(snp_id=0; snp_id<iFtrStats.vSnpStats.size(); ++snp_id)
    {
      vector<double> yPerm, gPerm;
      for (i=0; i<perm->size; ++i)
      {
	p = gsl_permutation_get (perm, i);
	if (! iFtrStats.vIsNa[p] && ! iFtrStats.vSnpStats[snp_id].vIsNa[i])
	{
	  yPerm.push_back (iFtrStats.y_init[p]);
	  gPerm.push_back (iFtrStats.vSnpStats[snp_id].g_init[i]);
	}
      }
      if (needQnorm)
	qnorm (yPerm);
      if (! calcSpearman)
      {
	ols ("", "", gPerm, yPerm, &betahat, &sebetahat, &sigmahat,
	     &betaPvalPerm, &R2, 0);
	if (betaPvalPerm < minBetaPvalPerm)
	  minBetaPvalPerm = betaPvalPerm;
      }
      else
      {
	gsl_vector_const_view gsl_gPerm =
	  gsl_vector_const_view_array (&gPerm[0], gPerm.size());
	gsl_vector_const_view gsl_yPerm =
	  gsl_vector_const_view_array (&yPerm[0], yPerm.size());
	rsPvalPerm = my_stats_correlation_spearman (gsl_gPerm.vector.data, 1,
						    gsl_yPerm.vector.data, 1,
						    gsl_gPerm.vector.size);
      rsZscore = sqrt((gPerm.size() - 3) / 1.06) * 1/2
	* (log(1 + rsPvalPerm) - log(1 - rsPvalPerm));
      rsPvalPerm = 2 * gsl_cdf_ugaussian_Q (fabs(rsZscore));
      if (rsPvalPerm < minRsPvalPerm)
	minRsPvalPerm = rsPvalPerm;
      }
    }
    
    if (! calcSpearman)
    {
      resPerm.push_back (minBetaPvalPerm);
      if (minBetaPvalPerm <= minBetaPval)
	++iFtrStats.betaPermPval;
      if (trickPerm && iFtrStats.betaPermPval == 11)
	break;
    }
    else
    {
      resPerm.push_back (minRsPvalPerm);
      if (minRsPvalPerm <= minRsPval)
	++iFtrStats.rsPermPval;
      if (trickPerm && iFtrStats.rsPermPval == 11)
	break;
    }
  }
  
  if (!calcSpearman)
  {
    if (resPerm.size() == nbPermutations)
      iFtrStats.betaPermPval /= (resPerm.size() + 1);
    else
      iFtrStats.betaPermPval = gsl_ran_flat (rng, 
					     (11 / ((double) (resPerm.size() + 2))),
					     (11 / ((double) (resPerm.size() + 1))));
  }
  else
  {
    if (resPerm.size() == nbPermutations)
      iFtrStats.rsPermPval /= (resPerm.size() + 1);
    else
      iFtrStats.rsPermPval = gsl_ran_flat (rng, 
					   (11 / ((double) (resPerm.size() + 2))),
					   (11 / ((double) (resPerm.size() + 1))));
  }
  
  gsl_permutation_free (perm);
}

/** \brief Compute the permutation P-value at the feature level using, as SNP-
 *  level test statistic, the  P-values on beta=0.
 *
 *  T_i,j: test statistic for SNP i and gene j
 *  T_min,j = min_i (T_i,j)
 *  T~_i,j,k: test statistic for SNP i, gene j and permutation k
 *  T~_min,j,k = min_i (T~_i,j,k)
 *  feature-level Pval = (#{k: T~_min,j,k <= T_min,j} + 1) / (#permutations + 1)
 */
void
computePermutationPvaluesAtFeatureLevelParallel (
  FtrStats & iFtrStats,
  const size_t nbPermutations,
  const bool needQnorm,
  const double minBetaPval,
  const bool calcSpearman,
  const double minRsPval,
  const int nbThreads,
  vector<gsl_rng *> & vRngs,
  const int verbose)
{
  long int perm_id;
  int t;
  vector<double> resPerm (nbPermutations, 1);
  vector<gsl_permutation *> vPerms; // for the indices of the phenotypes vector
  
  if (verbose > 0)
  {
    printf ("perform %zu permutations on the phenotypes (%d threads) ...\n",
	    nbPermutations, nbThreads);
    fflush (stdout);
  }
  
  for (t=0; t<nbThreads; ++t)
  {
    gsl_permutation * perm = gsl_permutation_calloc (iFtrStats.y_init.size());
    if (perm == NULL)
    {
      cerr << "ERROR: can't allocate memory for the permutation #" << t+1 << endl;
      exit (1);
    }
    vPerms.push_back (perm);
  }
  
#pragma omp parallel shared(resPerm) private(perm_id)
  {
#pragma omp for nowait
    for(perm_id=0; perm_id<(long int)nbPermutations; ++perm_id)
    {
      int tid = omp_get_thread_num();
      double minBetaPvalPerm = 1, betaPvalPerm, betahat, sebetahat, sigmahat, R2,
	minRsPvalPerm = 1, rsPvalPerm, rsZscore;
      gsl_ran_shuffle (vRngs[tid], vPerms[tid]->data, vPerms[tid]->size, sizeof(size_t));
      
      for(size_t snp_id=0; snp_id<iFtrStats.vSnpStats.size(); ++snp_id)
      {
	size_t p;
	vector<double> yPerm, gPerm;
	for (size_t i=0; i<vPerms[tid]->size; ++i)
	{
	  p = gsl_permutation_get (vPerms[tid], i);
	  if (! iFtrStats.vIsNa[p] && ! iFtrStats.vSnpStats[snp_id].vIsNa[i])
	  {
	    yPerm.push_back (iFtrStats.y_init[p]);
	    gPerm.push_back (iFtrStats.vSnpStats[snp_id].g_init[i]);
	  }
	}
	if (needQnorm)
	  qnorm (yPerm);
	if (! calcSpearman)
	{
	  ols ("", "", gPerm, yPerm, &betahat, &sebetahat, &sigmahat,
	       &betaPvalPerm, &R2, 0);
	  if (betaPvalPerm < minBetaPvalPerm)
	    minBetaPvalPerm = betaPvalPerm;
	}
	else
	{
	  gsl_vector_const_view gsl_gPerm =
	    gsl_vector_const_view_array (&gPerm[0], gPerm.size());
	  gsl_vector_const_view gsl_yPerm =
	    gsl_vector_const_view_array (&yPerm[0], yPerm.size());
	  rsPvalPerm = my_stats_correlation_spearman (gsl_gPerm.vector.data, 1,
						      gsl_yPerm.vector.data, 1,
						      gsl_gPerm.vector.size);
	  rsZscore = sqrt((gPerm.size() - 3) / 1.06) * 1/2
	    * (log(1 + rsPvalPerm) - log(1 - rsPvalPerm));
	  rsPvalPerm = 2 * gsl_cdf_ugaussian_Q (fabs(rsZscore));
	  if (rsPvalPerm < minRsPvalPerm)
	    minRsPvalPerm = rsPvalPerm;
	}
      }
      
      if (! calcSpearman)
	resPerm[perm_id] = minBetaPvalPerm;
      else
	resPerm[perm_id] = minRsPvalPerm;
    }
  }
  
  for (t=0; t<nbThreads; ++t)
    gsl_permutation_free (vPerms[t]);
  
  if (! calcSpearman)
  {
    iFtrStats.betaPermPval = 1;
    for (perm_id=0; perm_id<(long int)nbPermutations; ++perm_id)
      if (resPerm[perm_id] <= minBetaPval)
	++iFtrStats.betaPermPval;
    iFtrStats.betaPermPval /= (nbPermutations + 1);
  }
  else
  {
    iFtrStats.rsPermPval = 1;
    for (perm_id=0; perm_id<(long int)nbPermutations; ++perm_id)
      if (resPerm[perm_id] <= minRsPval)
	++iFtrStats.rsPermPval;
    iFtrStats.rsPermPval /= (nbPermutations + 1);
  }
}

void
computeSummaryStatsForOneFeature (
  FtrStats & iFtrStats,
  ifstream & genoStream,
  map<string, streampos> & mSnpNameCoord2Pos,
  const vector<bool> & vIdxSamplesToSkip,
  const size_t & nbPermutations,
  const bool & needQnorm,
  const bool & calcSpearman,
  const int & nbThreads,
  const bool & trickPerm,
  gsl_rng * rng,
  const int & verbose)
{
  size_t snp_id, i;
  vector<double> y, g;
  size_t nbSamples = count (vIdxSamplesToSkip.begin(),
			    vIdxSamplesToSkip.end(),
			    false);
  SnpStats * pt_iSnpStats;
  double minBetaPval = 1, minRsPval = 1;
  
  for (snp_id = 0; snp_id < iFtrStats.vSnpStats.size(); ++snp_id)
  {
    pt_iSnpStats = &(iFtrStats.vSnpStats[snp_id]);
    SnpStats_init (*pt_iSnpStats, genoStream, mSnpNameCoord2Pos,
		   vIdxSamplesToSkip, verbose-1);
    if (verbose > 1)
      printf ("ftr=%s snp=%s\n", iFtrStats.name.c_str(),
	      pt_iSnpStats->name.c_str());
    
    y.clear();
    g.clear();
    for (i = 0; i < nbSamples; ++i)
      if (! iFtrStats.vIsNa[i] && ! pt_iSnpStats->vIsNa[i])
      {
	y.push_back (iFtrStats.y_init[i]);
	g.push_back (pt_iSnpStats->g_init[i]);
      }
    pt_iSnpStats->n = g.size();
    if (needQnorm)
      qnorm (y);
    
    if (! calcSpearman)
    {
      ols (iFtrStats.name, pt_iSnpStats->name, g, y, &(pt_iSnpStats->betahat),
	   &(pt_iSnpStats->sebetahat), &(pt_iSnpStats->sigmahat),
	   &(pt_iSnpStats->betaPval), &(pt_iSnpStats->R2), verbose-1);
      if (pt_iSnpStats->betaPval < minBetaPval)
	minBetaPval = pt_iSnpStats->betaPval;
    }
    else
    {
      gsl_vector_const_view gsl_g = gsl_vector_const_view_array (&g[0],
								 g.size());
      gsl_vector_const_view gsl_y = gsl_vector_const_view_array (&y[0],
								 y.size());
      pt_iSnpStats->rs = my_stats_correlation_spearman (gsl_g.vector.data, 1,
						    gsl_y.vector.data, 1,
						    g.size());
      pt_iSnpStats->rsZscore = sqrt((g.size() - 3) / 1.06) * 1/2
	* (log(1 + pt_iSnpStats->rs) - log(1 - pt_iSnpStats->rs));
      pt_iSnpStats->rsPval = 2 * gsl_cdf_ugaussian_Q (fabs(pt_iSnpStats->rsZscore));
      if (pt_iSnpStats->rsPval < minRsPval)
	minRsPval = pt_iSnpStats->rsPval;
    }
  }
  
  if (nbPermutations > 0)
    computePermutationPvaluesAtFeatureLevel (iFtrStats,
					     nbPermutations,
					     needQnorm,
					     minBetaPval,
					     calcSpearman,
					     minRsPval,
					     trickPerm,
					     rng,
					     verbose);
}

void
computeSummaryStatsForOneFeatureParallel (
  FtrStats & iFtrStats,
  const string & genoFile,
  map<string, streampos> & mSnpNameCoord2Pos,
  const vector<bool> & vIdxSamplesToSkip,
  const size_t & nbPermutations,
  const bool & needQnorm,
  const bool & calcSpearman,
  const int & nbThreads,
  vector<gsl_rng *> & vRngs,
  const int & verbose)
{
  long int snp_id;
  int t;
  size_t nbSamples = count (vIdxSamplesToSkip.begin(),
			    vIdxSamplesToSkip.end(),
			    false);
  
  // need one file handle per thread
  vector<ifstream *> vPtGenoStreams;
  for (t=0; t<nbThreads; ++t)
  {
    ifstream * pt_genoStream = new ifstream;
    pt_genoStream->open(genoFile.c_str());
    if (! pt_genoStream->is_open())
    {
      cerr << "ERROR: can't open file " << genoFile
	   << " for thread " << t << endl;
      exit (1);
    }
    vPtGenoStreams.push_back (pt_genoStream);
  }
  
#pragma omp parallel private(snp_id)
  {
#pragma omp for nowait
    for (snp_id = 0; snp_id < (long int) iFtrStats.vSnpStats.size(); ++snp_id)
    {
      SnpStats * pt_iSnpStats = &(iFtrStats.vSnpStats[snp_id]);
      SnpStats_init (*pt_iSnpStats, *(vPtGenoStreams[omp_get_thread_num()]),
		     mSnpNameCoord2Pos, vIdxSamplesToSkip, verbose-1);
      if (verbose > 1)
	printf ("ftr=%s snp=%s\n", iFtrStats.name.c_str(),
		pt_iSnpStats->name.c_str());
      
      vector<double> y, g;
      for (size_t i = 0; i < nbSamples; ++i)
	if (! iFtrStats.vIsNa[i] && ! pt_iSnpStats->vIsNa[i])
	{
	  y.push_back (iFtrStats.y_init[i]);
	  g.push_back (pt_iSnpStats->g_init[i]);
	}
      pt_iSnpStats->n = g.size();
      if (needQnorm)
	qnorm (y);
      
      if (! calcSpearman)
      {
	ols (iFtrStats.name, pt_iSnpStats->name, g, y, &(pt_iSnpStats->betahat),
	     &(pt_iSnpStats->sebetahat), &(pt_iSnpStats->sigmahat),
	     &(pt_iSnpStats->betaPval), &(pt_iSnpStats->R2), verbose-1);
      }
      else
      {
	gsl_vector_const_view gsl_g = gsl_vector_const_view_array (&g[0],
								   g.size());
	gsl_vector_const_view gsl_y = gsl_vector_const_view_array (&y[0],
								   y.size());
	pt_iSnpStats->rs = my_stats_correlation_spearman (gsl_g.vector.data, 1,
							  gsl_y.vector.data, 1,
							  g.size());
	pt_iSnpStats->rsZscore = sqrt((g.size() - 3) / 1.06) * 1/2
	  * (log(1 + pt_iSnpStats->rs) - log(1 - pt_iSnpStats->rs));
	pt_iSnpStats->rsPval = 2 * gsl_cdf_ugaussian_Q (fabs(pt_iSnpStats->rsZscore));
      }
    }
  }
  
  for (t=0; t<nbThreads; ++t)
  {
    vPtGenoStreams[t]->close ();
    delete vPtGenoStreams[t];
  }
  vPtGenoStreams.clear();
  
  if (nbPermutations > 0)
  {
    double minBetaPval = 1, minRsPval = 1;
    if (! calcSpearman)
    {
      for (snp_id = 0; snp_id < (long int) iFtrStats.vSnpStats.size(); ++snp_id)
	if (iFtrStats.vSnpStats[snp_id].betaPval < minBetaPval)
	  minBetaPval = iFtrStats.vSnpStats[snp_id].betaPval;
    }
    else
    {
      for (snp_id = 0; snp_id < (long int) iFtrStats.vSnpStats.size(); ++snp_id)
	if (iFtrStats.vSnpStats[snp_id].rsPval < minRsPval)
	  minRsPval = iFtrStats.vSnpStats[snp_id].rsPval;
    }
    computePermutationPvaluesAtFeatureLevelParallel (iFtrStats,
						     nbPermutations,
						     needQnorm,
						     minBetaPval,
						     calcSpearman,
						     minRsPval,
						     nbThreads,
						     vRngs,
						     verbose);
  }
}

void
computeAndWriteSummaryStatsFtrPerFtr (
  const string & genoFile,
  const string & phenoFile,
  const string & outFile,
  const string & ftrCoordsFile,
  const string & linksFile,
  const string & anchor,
  const size_t & lenCis,
  const vector<string> & vFtrsToKeep,
  const vector<string> & vSnpsToKeep,
  const vector<string> & vSamplesToSkip,
  const string & chrToKeep,
  const bool & needQnorm,
  const double & minMaf,
  const size_t & nbPermutations,
  const bool & calcSpearman,
  const int & nbThreads,
  const vector<size_t> & vPermIndices,
  const bool & trickPerm,
  const int & seed,
  const int & verbose)
{
  FtrStats iFtrStats;
  string linePheno, lineLinks;
  size_t nbFtrs = 0, nbAnalyzedPairs = 0, nbAnalyzedFtrs = 0;
  vector<string> tokens, vSamples, vCisSnps;
  vector<bool> vIdxSamplesToSkip;
  map<string, streampos> mSnpNameCoord2Pos;
  ifstream phenoStream, genoStream, ftrCoordsStream, linksStream;
  ofstream outStream;
  gsl_permutation * perm = 0;
  
  // open input files
  phenoStream.open(phenoFile.c_str());
  if (! phenoStream.is_open())
  {
    cerr << "ERROR: can't open file " << phenoFile << endl;
    exit (1);
  }
  genoStream.open(genoFile.c_str());
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  if (! linksFile.empty())
  {
    linksStream.open(linksFile.c_str());
    if (! linksStream.is_open())
    {
      cerr << "ERROR: can't open file " << linksFile << endl;
      exit (1);
    }
  }
  ftrCoordsStream.open(ftrCoordsFile.c_str());
  if (! ftrCoordsStream.is_open())
  {
    cerr << "ERROR: can't open file " << ftrCoordsFile << endl;
    exit (1);
  }
  
  // open output file and write the header line
  outStream.open(outFile.c_str());
  if (! outStream.is_open())
  {
    cerr << "ERROR: can't open file " << outFile << endl;
    exit (1);
  }
  outStream << "ftr chrF start end snp chrS coord maf n";
  if (! calcSpearman)
  {
    outStream << " betahat sebetahat sigmahat betaPval R2";
    if (nbPermutations > 0)
      outStream << " betaPermPvalFtr";
  }
  else
  {
    outStream << " rs rsZscore rsPval";
    if (nbPermutations > 0)
      outStream << " rsPermPvalFtr";
  }
  outStream << endl;
  
  // index the genotype file
  indexSnps (genoStream, vSnpsToKeep, chrToKeep, vSamplesToSkip, vSamples,
	     vIdxSamplesToSkip, minMaf, mSnpNameCoord2Pos, verbose);
  
  if (verbose > 0)
  {
    cout << "look for association between each pair feature-SNP";
    if (nbPermutations > 0)
    {
      cout << nbPermutations << " permutation"
	   << (nbPermutations > 1 ? "s" : "");
      if (trickPerm)
	cout << ", with trick";
      else
	cout << ", without trick";
      if (nbThreads <= 1)
	cout << ", " << nbThreads << " thread"
	     << (nbThreads > 1 ? "s" : "");
      cout << ")";
    }
    cout << " ..." << endl;
    fflush (stdout);
  }
  
  // read header line of the phenotype file
  getline (phenoStream, linePheno);
  if (linePheno.find('\t') != string::npos)
    split (linePheno, '\t', tokens);
  else
    split (linePheno, ' ', tokens);
  if (tokens.size() != vSamples.size())
  {
    cerr << "ERROR: different number of samples between genotype"
	 << " and phenotype files" << endl;
    exit (1);
  }
  
  // check the samples are in same order with genotype file
  for(size_t colIdx = 0; colIdx < tokens.size(); ++colIdx)
    if (tokens[colIdx].compare(vSamples[colIdx]) != 0)
    {
      cerr << "ERROR: samples are not in same order between genotype"
	   << " and phenotype files" << endl;
      exit (1);
    }
  
  // initialize RNG (only if permutations)
  gsl_rng * rng = NULL;
  vector<gsl_rng *> vRngs;
  if (nbPermutations > 0)
  {
    gsl_rng_env_setup();
    if (nbThreads <= 1)
    {
      rng = gsl_rng_alloc (gsl_rng_default);
      if (rng == NULL)
      {
	cerr << "ERROR: can't allocate memory for the RNG" << endl;
	exit (1);
      }
      gsl_rng_set (rng, seed);
    }
    else
    {
      for (int t=0; t<nbThreads; ++t)
      {
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default);
	if (rng == NULL)
	{
	  cerr << "ERROR: can't allocate memory for the RNG #" << t+1 << endl;
	  exit (1);
	}
	gsl_rng_set (rng, seed * (t+1));
	vRngs.push_back (rng);
      }
    }
  }
  
  // initialize the permutation struct if required
  if (! vPermIndices.empty())
  {
    if (vPermIndices.size() != vSamples.size())
    {
      cerr << "ERROR: wrong number of permutation indices ("
	   << vPermIndices.size() << " versus " << vSamples.size()
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
  
  // for each feature
  while (phenoStream.good())
  {
    // initialize it
    FtrStats_init (iFtrStats, phenoStream, ftrCoordsStream, chrToKeep,
		   vIdxSamplesToSkip, vFtrsToKeep, perm, verbose-1);
    if (iFtrStats.name.empty())
      continue;
    ++nbFtrs;
    
    // retrieve its SNPs in cis
    if (! linksFile.empty())
      FtrStats_getCisSnps (iFtrStats, linksStream, mSnpNameCoord2Pos);
    else
      FtrStats_getCisSnps (iFtrStats, genoStream, mSnpNameCoord2Pos,
			   anchor, lenCis);
    if (iFtrStats.vSnpStats.empty())
      continue;
    
    if (verbose > 1)
    {
      printf ("analyzing feature %s (%s, %zu cis SNPs) ...\n",
	      iFtrStats.name.c_str(), iFtrStats.chr.c_str(),
	      iFtrStats.vSnpStats.size());
      fflush (stdout);
    }
    ++nbAnalyzedFtrs;
    
    // loop over SNPs in cis
    if (nbThreads <= 1)
      computeSummaryStatsForOneFeature (iFtrStats, genoStream,
					mSnpNameCoord2Pos, vIdxSamplesToSkip,
					nbPermutations, needQnorm,
					calcSpearman, nbThreads, trickPerm,
					rng, verbose-1);
    else
      computeSummaryStatsForOneFeatureParallel (iFtrStats, genoFile,
						mSnpNameCoord2Pos,
						vIdxSamplesToSkip,
						nbPermutations, needQnorm,
						calcSpearman, nbThreads,
						vRngs, verbose-1);
    
    // write the results (one line per SNP)
    if (iFtrStats.vSnpStats.size() > 0)
    {
      FtrStats_write (iFtrStats, outStream, nbPermutations, calcSpearman);
      nbAnalyzedPairs += iFtrStats.vSnpStats.size();
    }
  }
  
  if (! vPermIndices.empty())
    gsl_permutation_free (perm);
  gsl_rng_free (rng);
  for (size_t t=0; t<vRngs.size(); ++t)
    gsl_rng_free (vRngs[t]);
  phenoStream.close();
  genoStream.close();
  outStream.close();
  linksStream.close();
  ftrCoordsStream.close();
  
  if (verbose > 0)
  {
    cout << "nb of features: " << nbFtrs << endl;
    cout << "nb of analyzed features: " << nbAnalyzedFtrs << endl;
    cout << "nb of analyzed feature-SNP pairs: " << nbAnalyzedPairs << endl;
    cout << "results written in " << outFile << endl;
  }
}

#ifdef GET_SUMSTATS_MAIN

/** \brief Display the help on stdout.
*/
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " assesses associations between genotypes and phenotypes" << endl
       << "(one genetic variant per univariate phenotype) and returns the summary stats" << endl
       << "betahat, se(betahat) and sigmahat, as well as the P-value"
       << " for H0:\"beta=0\"." << endl
       << "Spearman correlation coefficient and permutation P-values at the SNP and" << endl
       << "gene level can also be computed." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS]..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (default=1)" << endl
       << "  -g, --geno\tfile with genotypes in IMPUTE format (delimiter=<space/tab>)" << endl
       << "\t\ta header line with sample names is required" << endl
       << "\t\tsamples in columns should be in same order as in phenotype file" << endl
       << "  -p, --pheno\tfile with phenotypes" << endl
       << "\t\trow 1 for sample names, column 1 for feature names" << endl
       << "\t\tdelimiter=<space/tab>" << endl
       << "  -o, --output\tfile that will contain the summary stats" << endl
       << "      --fcoord\tBED file with the features coordinates" << endl
       << "\t\tfeatures should be in same order than in phenotype file" << endl
       << "  -l, --links\tfile with links between genes and SNPs" << endl
       << "\t\tcustom format: feature<space/tab>SNP|coord" << endl
       << "\t\tfeatures should be in same order than in phenotype file" << endl
       << "\t\tuseful to focus on genetic variants in cis (use windowBed)" << endl
       << "  -a, --anchor\tfor the cis region if links is not specified" << endl
       << "\t\tdefault=TSS+TES, can also be only TSS" << endl
       << "      --cis\tlength of half of the cis region (in bp)" << endl
       << "\t\tapart from the anchor(s), default=100000" << endl
       << "  -c, --chr\tname of the chromosome to analyze (eg. 'chr21')" << endl
       << "  -f, --ftr\tfile with a list of features to analyze" << endl
       << "\t\tone feature name per line" << endl
       << "  -s, --snp\tfile with a list of SNPs to analyze" << endl
       << "\t\tone SNP identifier per line" << endl
       << "  -d, --discard\tfile with a list of samples to discard" << endl
       << "\t\tone individual per line, should match header of phenotype file" << endl
       << "  -q, --qnorm\tquantile-normalize phenotypes to a standard normal" << endl
       << "\t\tvery useful when using linear regression (versus Spearman coef)" << endl
       << "  -m, --maf\tthreshold for the minor allele frequency (default=0.0)" << endl
       << "\t\twhatever the option, the MAF will still be computed and saved" << endl
       << "  -P, --perm\tnumber of phenotype permutations at each feature" << endl
       << "\t\tdefault=0, recommended=10000" << endl
       << "      --trick\tstop after the tenth permutation for which the test statistic" << endl
       << "\t\tis better than or equal to the true value, and sample from a" << endl
       << "\t\tuniform between 11/(nbPermSoFar+2) and 11/(nbPermSoFar+1)" << endl
       << "\t\tnot used when nbThreads > 1" << endl
       << "  -S, --sp\tcompute the Spearman rank correlation coefficient (and Z score)" << endl
       << "\t\tinstead of performing linear regressions" << endl
       << "      --seed\tseed for the random number generator" << endl
       << "\t\tby default, the RNG is initialized via microseconds from epoch" << endl
       << "  -t, --thread\tnumber of threads (default=1)" << endl
       << "\t\tused for SNPs in cis of the same feature, and for permutations" << endl
       << "      --permf\tfile with a permutation of the samples" << endl
       << "\t\tshould start at 0" << endl
       << "\t\tshould have as many lines as the number of samples in --pheno" << endl
       << endl
       << "Examples:" << endl
       << "  " << argv[0] << " -g <genotypes> -p <phenotypes> -o <output> \\" << endl
       << "  -l <link> --fcoord <bed>" << endl
       << endl
       << "Remarks:" << endl
       << "  Samples with missing phenotypes (NA) and/or genotypes (0 0 0) are skipped." << endl
       << "For non-variable genotypes, a zero effect size is returned, along with an" << endl
       << "infinite std error and a P-value of 1." << endl;
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
  string & genoFile,
  string & phenoFile,
  string & outFile,
  string & ftrCoordsFile,
  string & linksFile,
  string & chrToKeep,
  string & ftrsFile,
  string & snpsFile,
  string & samplesFile,
  double & minMaf,
  size_t & nbPermutations,
  bool & needQnorm,
  bool & calcSpearman,
  int & nbThreads,
  string & permFile,
  string & anchor,
  size_t & lenCis,
  bool & trickPerm,
  int & seed,
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
	{"output", required_argument, 0, 'o'},
	{"fcoord", required_argument, 0, 0},
	{"links", required_argument, 0, 'l'},
	{"chr", required_argument, 0, 'c'},
	{"ftr", required_argument, 0, 'f'},
	{"snp", required_argument, 0, 's'},
	{"discard", required_argument, 0, 'd'},
	{"qnorm", no_argument, 0, 'q'},
	{"maf", required_argument, 0, 'm'},
	{"perm", required_argument, 0, 'P'},
	{"sp", no_argument, 0, 'S'},
	{"thread", required_argument, 0, 't'},
	{"permf", required_argument, 0, 0},
	{"anchor", required_argument, 0, 'a'},
	{"cis", required_argument, 0, 0},
	{"trick", no_argument, 0, 0},
	{"seed", required_argument, 0, 0},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:o:l:c:f:s:d:qm:P:St:a:",
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
      if (strcmp(long_options[option_index].name, "permf") == 0)
      {
	permFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "cis") == 0)
      {
	lenCis = atol(optarg);
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
      phenoFile = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'l':
      linksFile = optarg;
      break;
    case 'c':
      chrToKeep = optarg;
      break;
    case 'f':
      ftrsFile = optarg;
      break;
    case 's':
      snpsFile = optarg;
      break;
    case 'd':
      samplesFile = optarg;
      break;
    case 'q':
      needQnorm = true;
      break;
    case 'm':
      minMaf = atof(optarg);
      break;
    case 'P':
      nbPermutations = atol(optarg);
      break;
    case 'S':
      calcSpearman = true;
      break;
    case 't':
      nbThreads = atoi(optarg);
      break;
    case 'a':
      anchor = optarg;
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
  if (phenoFile.empty())
  {
    fprintf (stderr, "ERROR: missing file with phenotypes (-p).\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (phenoFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", phenoFile.c_str());
    help (argv);
    exit (1);
  }
  if (outFile.empty())
  {
    fprintf (stderr, "ERROR: missing output file (-o).\n\n");
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
  if (linksFile.empty() && anchor.empty())
  {
    fprintf (stderr, "ERROR: missing links file and anchor (-l and -a).\n\n");
    help (argv);
    exit (1);
  }
  if (! linksFile.empty() && ! doesFileExist (linksFile))
  {
    fprintf (stderr, "ERROR: can't find file '%s'.\n\n", linksFile.c_str());
    help (argv);
    exit (1);
  }
  if (minMaf < 0 || minMaf > 1)
  {
    fprintf (stderr, "ERROR: min MAF should be between 0 and 1 (-m).\n\n");
    help (argv);
    exit (1);
  }
  if (nbThreads > 1)
  {
    nbThreads = min(nbThreads, omp_get_max_threads());
    omp_set_num_threads (nbThreads);
  }
  if (nbPermutations > 0)
    if (seed < 0)
      seed = (int) getSeed();
}

int main (int argc, char ** argv)
{
  string genoFile, phenoFile, outFile, ftrCoordsFile, linksFile,
    ftrsFile, snpsFile, samplesFile, chrToKeep = "", permFile,
    anchor = "TSS+TES";
  double minMaf = 0.0;
  size_t nbPermutations = 0, lenCis = 100000;
  bool needQnorm = false, calcSpearman = false, trickPerm = false;
  int verbose = 1, nbThreads = 1, seed = -1;
  
  parse_args (argc, argv, genoFile, phenoFile, outFile, ftrCoordsFile,
	      linksFile, chrToKeep, ftrsFile, snpsFile, samplesFile,
	      minMaf, nbPermutations, needQnorm, calcSpearman, nbThreads,
	      permFile, anchor, lenCis, trickPerm, seed, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl;
  }
  
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsFile, verbose);
  vector<string> vSamplesToSkip = loadOneColumnFile (samplesFile, verbose);
  vector<size_t> vPermIndices = loadOneColumnFileAsNumbers (permFile, verbose);
  
  computeAndWriteSummaryStatsFtrPerFtr (genoFile,
					phenoFile,
					outFile,
					ftrCoordsFile,
					linksFile,
					anchor,
					lenCis,
					vFtrsToKeep,
					vSnpsToKeep,
					vSamplesToSkip,
					chrToKeep,
					needQnorm,
					minMaf,
					nbPermutations,
					calcSpearman,
					nbThreads,
					vPermIndices,
					trickPerm,
					seed,
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

#endif
