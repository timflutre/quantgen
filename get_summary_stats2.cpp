/** \file get_summary_stats2.cpp
 *
 *  `get_summary_stats2' computes summary statistics of association between genotypes and phenotypes.
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
 *  g++ -Wall -O3 -DGET_SUMSTATS_MAIN get_summary_stats2.cpp utils.cpp -lgsl -lgslcblas -o get_summary_stats2
 *  help2man -o get_summary_stats2.man ./get_summary_stats2
 *  groff -mandoc get_summary_stats2.man > get_summary_stats2.ps
*/

#include <math.h>
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
#include <iomanip>
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
    double inv_xtx[2][2];
    inv_xtx[0][0] = gtg;
    inv_xtx[0][1] = - n * gm;
    inv_xtx[1][0] = - n * gm;
    inv_xtx[1][1] = n;
    double xty[2];
    xty[0] = n * ym;
    xty[1] = gty;
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

struct Snp
{
  string name; // eg. rs2205177
  string chr; // eg. chr21
  size_t coord; // 0-based coordinate
  vector<double> vGenos; // genotypes of samples
  vector<bool> vIsNa; // missing values
  double maf; // minor allele frequency
};

// assume both features are on the same chromosome
bool Snp_compByCoord (
  const vector<Snp>::const_iterator & it_iSnp1,
  const vector<Snp>::const_iterator & it_iSnp2)
{
  if (it_iSnp1->coord < it_iSnp2->coord)
    return true;
  return false;
}

void Snp_init (
  Snp & iSnp,
  const vector<string> & tokens,
  const vector<bool> & vIdxSamplesToSkip)
{
  if (tokens.size() != vIdxSamplesToSkip.size() * 3 + 5)
  {
    cerr << "ERROR: a line in genotype file has a wrong number of samples"
	 << endl;
    exit (1);
  }
  iSnp.chr = tokens[0];
  iSnp.name = tokens[1];
  iSnp.coord = atol (tokens[2].c_str());
  double maf = 0, AA, AB, BB;
  iSnp.vIsNa.assign (vIdxSamplesToSkip.size(), false);
  iSnp.vGenos.assign (vIdxSamplesToSkip.size(), 0);
  size_t j = 0;
  for (size_t i = 0; i < vIdxSamplesToSkip.size(); ++i)
  {
    if (vIdxSamplesToSkip[i])
      continue;
    AA = atof (tokens[5+3*i].c_str());
    AB = atof (tokens[5+3*i+1].c_str());
    BB = atof (tokens[5+3*i+2].c_str());
    if (AA == 0 && AB == 0 && BB == 0)
      iSnp.vIsNa[j] = true;
    else
    {
      iSnp.vGenos[j] = 0 * AA + 1 * AB + 2 * BB;
      maf += iSnp.vGenos[j];
    }
    ++j;
  }
  maf /= 2 * (vIdxSamplesToSkip.size()
	      - count (vIdxSamplesToSkip.begin(),
		       vIdxSamplesToSkip.end(),
		       true)
	      - count (iSnp.vIsNa.begin(),
		       iSnp.vIsNa.end(),
		       true));
  iSnp.maf = maf <= 0.5 ? maf : (1 - maf);
}

struct StatsFtr
{
  string ftrName;
  vector<string> vSnpNames;
  vector<double> vMafs;
  vector<size_t> n; // sample size
  vector<double> betahat; // MLE of beta
  vector<double> sebetahat; // standard error
  vector<double> sigmahat; // MLE of sigma
  vector<double> betaPval; // P-value of the test where H0:"beta=0" and H1:"beta!=0"
  vector<double> R2; // coef of determination (proportion of variance explained)
  vector<double> rs; // Spearman rank correlation coefficient
  vector<double> rsZscore; // using Fisher transformation
  vector<double> rsPval; // using rsZscore
  double betaPermPval; // permutation P-value based on beta P-value
  double rsPermPval; // permutation P-value based on Spearman-derived P-value
};

void StatsFtr_reset (
  StatsFtr & iStatsFtr)
{
  iStatsFtr.ftrName.clear();
  iStatsFtr.vSnpNames.clear();
  iStatsFtr.vMafs.clear();
  iStatsFtr.n.clear();
  iStatsFtr.betahat.clear();
  iStatsFtr.sebetahat.clear();
  iStatsFtr.sigmahat.clear();
  iStatsFtr.betaPval.clear();
  iStatsFtr.R2.clear();
  iStatsFtr.rs.clear();
  iStatsFtr.rsZscore.clear();
  iStatsFtr.rsPval.clear();
}

void StatsFtr_init (
  StatsFtr & iStatsFtr,
  const string & ftrName)
{
  StatsFtr_reset (iStatsFtr);
  iStatsFtr.ftrName = ftrName;
}

void StatsFtr_addSnp (
  StatsFtr & iStatsFtr,
  const Snp & iSnp)
{
  iStatsFtr.vSnpNames.push_back (iSnp.name);
  iStatsFtr.vMafs.push_back (iSnp.maf);
  iStatsFtr.n.push_back (0);
  iStatsFtr.betahat.push_back (0);
  iStatsFtr.sebetahat.push_back (0);
  iStatsFtr.sigmahat.push_back (0);
  iStatsFtr.betaPval.push_back (0);
  iStatsFtr.R2.push_back (0);
  iStatsFtr.rs.push_back (0);
  iStatsFtr.rsZscore.push_back (0);
  iStatsFtr.rsPval.push_back (0);
}

void StatsFtr_write (
  const StatsFtr & iStatsFtr,
  ostream & outStream,
  const bool & calcSpearman,
  const size_t & nbPermutations)
{
  for (size_t i = 0; i < iStatsFtr.vSnpNames.size(); i++)
  {
    outStream << iStatsFtr.ftrName
	      << " " << iStatsFtr.vSnpNames[i]
	      << " " << iStatsFtr.vMafs[i]
	      << " " << iStatsFtr.n[i];
    if (! calcSpearman)
    {
      outStream << " " << iStatsFtr.betahat[i]
		<< " " << iStatsFtr.sebetahat[i]
		<< " " << iStatsFtr.sigmahat[i]
		<< " " << iStatsFtr.betaPval[i]
		<< " " << iStatsFtr.R2[i];
      if (nbPermutations > 0)
	outStream << " " << iStatsFtr.betaPermPval;
    }
    else
    {
      outStream << " " << iStatsFtr.rs[i]
		<< " " << iStatsFtr.rsZscore[i]
		<< " " << iStatsFtr.rsPval[i];
      if (nbPermutations > 0)
	outStream << " " << iStatsFtr.rsPermPval;
    }
    outStream << endl;
  }
}

struct FtrSingleS
{
  string name; // eg. ENSG00000182816
  string chr; // eg. chr21
  size_t start; // 0-based coordinate
  size_t end; // idem
  vector<double> vPhenos; // phenotypes of samples
  vector<bool> vIsNa; // missing values
};

// http://stackoverflow.com/a/2025952/597069
struct FtrSingleS_findByName
{
  string name;
  FtrSingleS_findByName (const string & name) : name(name) {}
  bool operator () (const FtrSingleS & iFtr) const
  {
    return (iFtr.name.compare(name) == 0);
  }
};

// assume both features are on the same chromosome
bool FtrSingleS_compByCoord (
  const vector<FtrSingleS>::iterator & it_iFtr1,
  const vector<FtrSingleS>::iterator & it_iFtr2)
{
  if (it_iFtr1->start < it_iFtr2->start)
    return true;
  else if (it_iFtr1->start == it_iFtr2->start && 
    it_iFtr1->end < it_iFtr2->end )
    return true;
  return false;
}

void FtrSingleS_init (
  FtrSingleS & iFtr,
  const vector<string> & tokens,
  const vector<bool> & vIdxSamplesToSkip,
  gsl_permutation * perm)
{
  if (tokens.size() != vIdxSamplesToSkip.size() + 1)
  {
    cerr << "ERROR: a line in phenotype file has a wrong number of samples"
	 << endl;
    exit (1);
  }
  iFtr.name = tokens[0];
  iFtr.vIsNa.assign (vIdxSamplesToSkip.size(), false);
  iFtr.vPhenos.assign (vIdxSamplesToSkip.size(), 0);
  size_t j = 0;
  string tok;
  for (size_t colIdx = 1; colIdx < tokens.size(); ++colIdx)
  {
    if(perm == NULL && vIdxSamplesToSkip[colIdx-1])
      continue;
    if (perm != NULL &&
	vIdxSamplesToSkip[gsl_permutation_get (perm, colIdx-1) + 1])
      continue;
    if (perm == NULL)
      tok = tokens[colIdx];
    else
      tok = tokens[gsl_permutation_get (perm, colIdx-1) + 1];
    if (tok.compare("NA") == 0)
      iFtr.vIsNa[j] = true;
    else
      iFtr.vPhenos[j] = atof (tok.c_str());
    ++j;
  }
}

void getAndCheckSamples (
  const string & genoFile,
  const string & phenoFile,
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

void loadFeatures (
  const string & phenoFile,
  const string & ftrCoordsFile,
  const vector<string> & vFtrsToKeep,
  const string & chrToKeep,
  const size_t & nbSamples,
  const vector<bool> & vIdxSamplesToSkip,
  gsl_permutation * perm,
  vector<FtrSingleS> & vFeatures,
  map<string, vector< vector<FtrSingleS>::iterator > > & mChr2VecItFtrs,
  const int & verbose)
{
  string line, tok;
  vector<string> tokens;
  
  // parse the phenotype file
  ifstream phenoStream;
  openFile (phenoFile, phenoStream);
  getline (phenoStream, line); // skip header
  while (phenoStream.good())
  {
    getline (phenoStream, line);
    if (line.empty())
      break;
    split (line, " \t", tokens);
    if (! vFtrsToKeep.empty() &&
	find (vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0])
	== vFtrsToKeep.end())
      continue;
    FtrSingleS iFtr;
    FtrSingleS_init (iFtr, tokens, vIdxSamplesToSkip, perm);
    vFeatures.push_back (iFtr);
  }
  phenoStream.close();
  
  // retrieve their coordinates from BED file
  ifstream ftrCoordsStream;
  openFile (ftrCoordsFile, ftrCoordsStream);
  while (ftrCoordsStream.good())
  {
    getline (ftrCoordsStream, line);
    if (line.empty())
      break;
    split (line, " \t", tokens);
    vector<FtrSingleS>::iterator it = find_if (vFeatures.begin(), vFeatures.end(), FtrSingleS_findByName (tokens[3]));
    if (it == vFeatures.end())
      continue;
    if (! chrToKeep.empty() && chrToKeep.compare(tokens[0]) != 0)
    {
      vFeatures.erase (it);
      continue;
    }
    
    it->chr = tokens[0];
    it->start = atol (tokens[1].c_str());
    it->end = atol (tokens[2].c_str());
    
    if (mChr2VecItFtrs.find(it->chr) == mChr2VecItFtrs.end())
      mChr2VecItFtrs.insert (make_pair (it->chr, vector< vector<FtrSingleS>::iterator > ()));
    mChr2VecItFtrs[it->chr].push_back (it);
  }
  ftrCoordsStream.close();
  
  // sort the features per chr
  for (map<string, vector< vector<FtrSingleS>::iterator > >::iterator it = mChr2VecItFtrs.begin();
       it != mChr2VecItFtrs.end(); ++it)
    sort (it->second.begin(), it->second.end(), FtrSingleS_compByCoord);
  
  if (vFeatures.size() == 0)
  {
    cerr << "ERROR: no feature to analyze" << endl;
    exit (1);
  }
  
  if (verbose > 0)
    cout << "nb of features: " << vFeatures.size() << endl
	 << "nb of chromosomes: " << mChr2VecItFtrs.size() << endl << flush;
}

void loadSnps (
  const string & genoFile,
  const vector<string> & vSnpsToKeep,
  const string & chrToKeep,
  const vector<bool> & vIdxSamplesToSkip,
  const double & minMaf,
  vector<Snp> & vSnps,
  map<string, vector< vector<Snp>::const_iterator > > & mChr2VecItSnps,
  const int & verbose)
{
  // parse the genotype file
  ifstream genoStream;
  openFile (genoFile, genoStream);
  string line;
  vector<string> tokens;
  getline (genoStream, line); // skip header
  while (genoStream.good())
  {
    getline (genoStream, line);
    if (line.empty())
      break;
    split (line, " \t", tokens);
    if (! vSnpsToKeep.empty() &&
	find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[1])
	== vSnpsToKeep.end())
      continue;
    if (! chrToKeep.empty() && chrToKeep.compare(tokens[0]) != 0)
      continue;
    Snp iSnp;
    Snp_init (iSnp, tokens, vIdxSamplesToSkip);
    if (iSnp.maf < minMaf)
      continue;
    vSnps.push_back (iSnp);
  }
  genoStream.close();
  if (verbose > 0)
    cout << "nb of SNPs: " << vSnps.size() << endl << flush;
  
  // organize SNPs per chr
  for (vector<Snp>::const_iterator it = vSnps.begin(); it != vSnps.end(); ++it)
  {
    if (mChr2VecItSnps.find(it->chr) == mChr2VecItSnps.end())
      mChr2VecItSnps.insert (make_pair (it->chr, vector< vector<Snp>::const_iterator > ()));
    mChr2VecItSnps[it->chr].push_back (it);
  }
  
  // sort the SNPs per chr
  for (map<string, vector< vector<Snp>::const_iterator > >::iterator it = mChr2VecItSnps.begin();
       it != mChr2VecItSnps.end(); ++it)
    sort (it->second.begin(), it->second.end(), Snp_compByCoord);
}

int isSnpInCis (
  const size_t & snpCoord,
  const size_t & ftrStart,
  const size_t & ftrEnd,
  const string & anchor,
  const size_t & lenCis)
{
  int res = -1;
  if (anchor.compare("FSS+FES") == 0)
  {
    if (((ftrStart >= lenCis &&
	 snpCoord >= ftrStart - lenCis) ||
	(ftrStart < lenCis &&
	 snpCoord >= 0)) &&
	snpCoord <= ftrEnd + lenCis)
      res = 0;
    else if (snpCoord > ftrEnd + lenCis)
      res = 1;
  }
  else if (anchor.compare("FSS") == 0)
  {
    if (((ftrStart >= lenCis &&
	 snpCoord >= ftrStart - lenCis) ||
	(ftrStart < lenCis &&
	 snpCoord >= 0)) &&
	snpCoord <= ftrStart + lenCis)
      res = 0;
    else if (snpCoord > ftrStart + lenCis)
      res = 1;
  }
  return res;
}

void initAllRngs(
  const size_t & nbPermutations,
  const int & trick,
  const size_t & seed,
  gsl_rng * & rngPerm,
  gsl_rng * & rngTrick)
{
  if (nbPermutations > 0)
  {
    gsl_rng_env_setup();
    rngPerm = gsl_rng_alloc (gsl_rng_default);
    if (rngPerm == NULL)
    {
      cerr << "ERROR: can't allocate memory for the RNG" << endl;
      exit (1);
    }
    gsl_rng_set (rngPerm, seed);
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
  }
}

void freeAllRngs(
  const size_t & nbPermutations,
  const int & trick,
  gsl_rng * & rngPerm,
  gsl_rng * & rngTrick)
{
  if (nbPermutations > 0)
  {
    gsl_rng_free (rngPerm);
    if (trick != 0)
      gsl_rng_free (rngTrick);
  }
}

void testOnePairFtrSnp (
  const FtrSingleS & iFtr,
  const Snp & iSnp,
  const bool & needQnorm,
  const bool & calcSpearman,
  StatsFtr & iStatsFtr)
{
  vector<double> y, g;
  size_t idx = iStatsFtr.vSnpNames.size() - 1;
  
  for (size_t i = 0; i < iFtr.vPhenos.size(); ++i)
    if (! iFtr.vIsNa[i] && ! iSnp.vIsNa[i])
    {
      y.push_back (iFtr.vPhenos[i]);
      g.push_back (iSnp.vGenos[i]);
    }
  iStatsFtr.n[idx] = g.size();
  if (needQnorm)
    qnorm (y);
  
  if (! calcSpearman)
    ols (g, y, iStatsFtr.betahat[idx],
	 iStatsFtr.sebetahat[idx], iStatsFtr.sigmahat[idx],
	 iStatsFtr.betaPval[idx], iStatsFtr.R2[idx]);
  else
  {
    gsl_vector_const_view gsl_g = gsl_vector_const_view_array (&g[0],
							       g.size());
    gsl_vector_const_view gsl_y = gsl_vector_const_view_array (&y[0],
							       y.size());
    iStatsFtr.rs[idx] = my_stats_correlation_spearman (gsl_g.vector.data, 1,
							gsl_y.vector.data, 1,
							g.size());
    iStatsFtr.rsZscore[idx] = sqrt((g.size() - 3) / 1.06) * 1/2
      * (log(1 + iStatsFtr.rs[idx]) - log(1 - iStatsFtr.rs[idx]));
    iStatsFtr.rsPval[idx] = 2 * gsl_cdf_ugaussian_Q (fabs(iStatsFtr.rsZscore[idx]));
  }
}

void getPermPvalAtFtrLevel (
  const FtrSingleS & iFtr,
  const vector< vector<Snp>::const_iterator > & vItSnps,
  const vector<size_t> & vIdxCisSnps,
  const bool & needQnorm,
  const bool & calcSpearman,
  const size_t & nbPermutations,
  const int & trick,
  gsl_rng * rngPerm,
  gsl_rng * rngTrick,
  StatsFtr & iStatsFtr,
  const int & verbose)
{
  double minBetaPval = 1, minRsPval = 1,
    betahat, sebetahat, sigmahat, R2;
  vector<double> resPerm;
  
  gsl_permutation * perm = NULL;
  perm = gsl_permutation_calloc (iFtr.vPhenos.size());
  if (perm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the permutation" << endl;
    exit (1);
  }
  
  if (! calcSpearman)
  {
    iStatsFtr.betaPermPval = 1;
    minBetaPval = *min_element (iStatsFtr.betaPval.begin(),
				iStatsFtr.betaPval.end());
  }
  else
  {
    iStatsFtr.rsPermPval = 1;
    minRsPval = *min_element (iStatsFtr.rsPval.begin(),
			      iStatsFtr.rsPval.end());
  }
  
  bool shuffleOnly = false;
  Snp iSnp;
  for(size_t permId = 0; permId < nbPermutations; ++permId)
  {
    double minBetaPvalPerm = 1, betaPvalPerm, minRsPvalPerm = 1, rsPvalPerm,
      rsZscore;
    gsl_ran_shuffle (rngPerm, perm->data, perm->size, sizeof(size_t));
    if (shuffleOnly)
      continue;
    // if (verbose > 0 && (permId+1) % 500 == 0)
    //   cout << setfill('0') << setw((int)floor(log10(nbPermutations))+1)
    // 	   << (permId+1) << "/" << nbPermutations << "\r" << flush;
    
    for(size_t snpId = 0; snpId < vIdxCisSnps.size(); ++snpId)
    {
      iSnp = *(vItSnps[vIdxCisSnps[snpId]]);
      vector<double> yPerm, gPerm;
      for (size_t i = 0; i < perm->size; ++i)
      {
	size_t p = gsl_permutation_get (perm, i);
	if (! iFtr.vIsNa[p] && ! iSnp.vIsNa[i])
	{
	  yPerm.push_back (iFtr.vPhenos[p]);
	  gPerm.push_back (iSnp.vGenos[i]);
	}
      }
      if (needQnorm)
	qnorm (yPerm);
      if (! calcSpearman)
      {
	ols (gPerm, yPerm, betahat, sebetahat, sigmahat, betaPvalPerm, R2);
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
	++iStatsFtr.betaPermPval;
      if (trick != 0 && iStatsFtr.betaPermPval == 11)
      {
	if (trick == 1)
	  break;
	else if (trick == 2)
	  shuffleOnly = true;
      }
    }
    else
    {
      resPerm.push_back (minRsPvalPerm);
      if (minRsPvalPerm <= minRsPval)
	++iStatsFtr.rsPermPval;
      if (trick != 0 && iStatsFtr.rsPermPval == 11)
      {
	if (trick == 1)
	  break;
	else if (trick == 2)
	  shuffleOnly = true;
      }
    }
  }
  
  if (!calcSpearman)
  {
    if (resPerm.size() == nbPermutations)
      iStatsFtr.betaPermPval /= (nbPermutations + 1);
    else
      iStatsFtr.betaPermPval = gsl_ran_flat (rngTrick, 
					     (11 / ((double) (resPerm.size() + 2))),
					     (11 / ((double) (resPerm.size() + 1))));
  }
  else
  {
    if (resPerm.size() == nbPermutations)
      iStatsFtr.rsPermPval /= (nbPermutations + 1);
    else
      iStatsFtr.rsPermPval = gsl_ran_flat (rngTrick, 
					   (11 / ((double) (resPerm.size() + 2))),
					   (11 / ((double) (resPerm.size() + 1))));
  }
  
  gsl_permutation_free (perm);
}

void testAllPairsFtrSnpInCis (
  vector<FtrSingleS> & vFeatures,
  const map<string, vector< vector<FtrSingleS>::iterator > > & mChr2VecItFtrs,
  const vector<Snp> & vSnps,
  map<string, vector< vector<Snp>::const_iterator > > & mChr2VecItSnps,
  const string & anchor,
  const size_t & lenCis,
  const bool & needQnorm,
  const bool & calcSpearman,
  const size_t & nbPermutations,
  const int & trick,
  const size_t & seed,
  vector<StatsFtr> & vStatsFtrs,
  const int & verbose)
{
  gsl_rng * rngPerm = NULL;
  gsl_rng * rngTrick = NULL;
  initAllRngs (nbPermutations, trick, seed, rngPerm, rngTrick);
  
  if (verbose > 0)
  {
    cout << "look for association between each pair feature-SNP ..." << endl
	 << "anchor=" << anchor << " lenCis=" << lenCis << endl << flush;
    if (nbPermutations > 0)
    {
      cout << "permutation"<< (nbPermutations > 1 ? "s=" : "=")
	   << nbPermutations
	   << ", seed=" << seed
	   << ", trick=" << trick
	   << endl << flush;
    }
  }
  
  // for each chr
  vector<FtrSingleS>::iterator itF;
  vector<Snp>::const_iterator itS;
  vector<size_t> vCounters = getCounters (vFeatures.size(), 5);
  size_t countFtrs = 0, nbAnalyzedPairs = 0;
  for (map<string, vector< vector<FtrSingleS>::iterator > >::const_iterator itC = mChr2VecItFtrs.begin();
       itC != mChr2VecItFtrs.end(); ++itC)
  {
    // for each ftr on this chr
    for (size_t i = 0; i < itC->second.size(); ++i)
    {
      itF = itC->second[i];
#ifdef DEBUG
      cout << itF->name << " " << itF->chr << " " << itF->start << " " << itF->end << endl << flush;
#endif
      StatsFtr iStatsFtr;
      StatsFtr_init (iStatsFtr, itF->name);
      ++countFtrs;
      if (verbose > 0)
	printCounter (countFtrs, vCounters);
      
      // for each SNP in cis on same chr
      vector<size_t> vIdxCisSnps;
      for (size_t j = 0; j < mChr2VecItSnps[itC->first].size(); ++j)
      {
	itS = mChr2VecItSnps[itC->first][j];
#ifdef DEBUG
	cout << itS->name << " " << itS->chr << " " << itS->coord << endl << flush;
#endif
	int inCis = isSnpInCis (itS->coord, itF->start, itF->end, anchor, lenCis);
	if (inCis == 1)
	  break;
	else if (inCis == -1)
	  continue;
	vIdxCisSnps.push_back (j);
	++nbAnalyzedPairs;
	StatsFtr_addSnp (iStatsFtr, *itS);
	
	testOnePairFtrSnp (*itF, *itS, needQnorm, calcSpearman, iStatsFtr);
      }
      
      // perform permutations
      if (nbPermutations > 0)
	getPermPvalAtFtrLevel (*itF, mChr2VecItSnps[itC->first], vIdxCisSnps, needQnorm, calcSpearman, nbPermutations, trick, rngPerm, rngTrick, iStatsFtr, verbose);
      
      vStatsFtrs.push_back (iStatsFtr);
    }
  }
  if (verbose > 0)
    cout << "nb of analyzed feature-SNP pairs: " << nbAnalyzedPairs << endl;
  
  freeAllRngs (nbPermutations, trick, rngPerm, rngTrick);
}

void testAllPairsFtrSnpInCisAndTrans (
  vector<FtrSingleS> & vFeatures,
  const vector<Snp> & vSnps,
  ofstream & outStream,
  const size_t & nbPermutations,
  const int & trick,
  const size_t & seed,
  const int & verbose)
{
  // TODO
}

void writeResults (
  const vector<StatsFtr> & vStatsFtrs,
  const string & outFile,
  const bool & calcSpearman,
  const size_t & nbPermutations,
  const int & verbose)
{
  ofstream outStream;
  openFile (outFile, outStream);
  
  outStream << "ftr snp maf n";
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
  
  for (vector<StatsFtr>::const_iterator it = vStatsFtrs.begin();
       it != vStatsFtrs.end(); ++it)
    StatsFtr_write (*it, outStream, calcSpearman, nbPermutations);
  
  if (verbose > 0)
    cout << "results written in " << outFile << endl;
}

void
computeAndWriteSummaryStatsFtrPerFtr (
  const string & genoFile,
  const string & phenoFile,
  const string & outFile,
  const string & ftrCoordsFile,
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
  const vector<size_t> & vPermIndices,
  const int & trick,
  const size_t & seed,
  const int & verbose)
{
  size_t nbSamples = 0;
  vector<bool> vIdxSamplesToSkip;
  gsl_permutation * perm = NULL;
  getAndCheckSamples (genoFile, phenoFile, vSamplesToSkip, vPermIndices, nbSamples, vIdxSamplesToSkip, perm, verbose);
  
  vector<FtrSingleS> vFeatures;
  map<string, vector< vector<FtrSingleS>::iterator > > mChr2VecItFtrs;
  loadFeatures (phenoFile, ftrCoordsFile, vFtrsToKeep, chrToKeep, nbSamples, vIdxSamplesToSkip, perm, vFeatures, mChr2VecItFtrs, verbose);
  if (! vPermIndices.empty())
    gsl_permutation_free (perm);
  
  vector<Snp> vSnps;
  map<string, vector< vector<Snp>::const_iterator > > mChr2VecItSnps;
  loadSnps (genoFile, vSnpsToKeep, chrToKeep, vIdxSamplesToSkip, minMaf, vSnps, mChr2VecItSnps, verbose);
  
  vector<StatsFtr> vStatsFtrs;
  testAllPairsFtrSnpInCis (vFeatures, mChr2VecItFtrs, vSnps, mChr2VecItSnps, anchor, lenCis, needQnorm, calcSpearman, nbPermutations, trick, seed, vStatsFtrs, verbose);
  
  writeResults (vStatsFtrs, outFile, calcSpearman, nbPermutations, verbose);
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
       << "\t\tthe first word on row 1 can be 'Id', followed by sample names" << endl
       << "  -o, --output\tfile that will contain the summary stats" << endl
       << "      --fcoord\tBED file with the features coordinates" << endl
       << "\t\tfeatures should be in same order than in phenotype file" << endl
       << "  -a, --anchor\tfeature boundary(ies) for the cis region" << endl
       << "\t\tdefault=FSS+FES, can also be only FSS" << endl
       << "\t\tif empty, all SNPs (cis and trans) are tested" << endl
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
       << "      --trick\tapply trick to speed-up permutations" << endl
       << "\t\tstop after the tenth permutation for which the test statistic" << endl
       << "\t\tis better than or equal to the true value, and sample from" << endl
       << "\t\ta uniform between 11/(nbPermsSoFar+2) and 11/(nbPermsSoFar+1)" << endl
       << "\t\tif '1', the permutations really stops" << endl
       << "\t\tif '2', all permutations are done but the test statistics are not computed" << endl
       << "\t\tallow to compare different test statistics on the same permutations" << endl
       << "  -S, --sp\tcompute the Spearman rank correlation coefficient (and Z score)" << endl
       << "\t\tinstead of performing linear regressions" << endl
       << "      --seed\tseed for the random number generator" << endl
       << "\t\tby default, the RNG is initialized via microseconds from epoch" << endl
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
  string & chrToKeep,
  string & ftrsFile,
  string & snpsFile,
  string & samplesFile,
  double & minMaf,
  size_t & nbPermutations,
  bool & needQnorm,
  bool & calcSpearman,
  string & permFile,
  string & anchor,
  size_t & lenCis,
  int & trick,
  size_t & seed,
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
	{"chr", required_argument, 0, 'c'},
	{"ftr", required_argument, 0, 'f'},
	{"snp", required_argument, 0, 's'},
	{"discard", required_argument, 0, 'd'},
	{"qnorm", no_argument, 0, 'q'},
	{"maf", required_argument, 0, 'm'},
	{"perm", required_argument, 0, 'P'},
	{"sp", no_argument, 0, 'S'},
	{"permf", required_argument, 0, 0},
	{"anchor", required_argument, 0, 'a'},
	{"cis", required_argument, 0, 0},
	{"trick", required_argument, 0, 0},
	{"seed", required_argument, 0, 0},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:o:c:f:s:d:qm:P:Sa:",
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
	lenCis = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "trick") == 0)
      {
	trick = atoi (optarg);
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
      genoFile = optarg;
      break;
    case 'p':
      phenoFile = optarg;
      break;
    case 'o':
      outFile = optarg;
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
      minMaf = atof (optarg);
      break;
    case 'P':
      nbPermutations = atol (optarg);
      break;
    case 'S':
      calcSpearman = true;
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
  ofstream outStream;
  openFile (outFile, outStream);
  outStream.close();
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
  if (anchor.empty())
  {
    fprintf (stderr, "ERROR: SNPs in trans not yet implemented, see --anchor.\n\n");
    help (argv);
    exit (1);
  }
  if (minMaf < 0 || minMaf > 1)
  {
    fprintf (stderr, "ERROR: min MAF should be between 0 and 1 (-m).\n\n");
    help (argv);
    exit (1);
  }
  if (trick != 0 && trick != 1 && trick != 2)
  {
    fprintf (stderr, "ERROR: unrecognized trick %d (--trick).\n\n", trick);
    help (argv);
    exit (1);
  }
  if (nbPermutations > 0 && seed == string::npos)
    seed = getSeed();
}

int main (int argc, char ** argv)
{
  string genoFile, phenoFile, outFile, ftrCoordsFile, ftrsFile,
    snpsFile, samplesFile, chrToKeep = "", permFile,
    anchor = "FSS+FES";
  double minMaf = 0.0;
  size_t nbPermutations = 0, seed = string::npos, lenCis = 100000;
  bool needQnorm = false, calcSpearman = false;
  int verbose = 1, trick = 0;
  
  parse_args (argc, argv, genoFile, phenoFile, outFile, ftrCoordsFile,
	      chrToKeep, ftrsFile, snpsFile, samplesFile, minMaf,
	      nbPermutations, needQnorm, calcSpearman, 
	      permFile, anchor, lenCis, trick, seed, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl << flush;
  }
  
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsFile, verbose);
  vector<string> vSamplesToSkip = loadOneColumnFile (samplesFile, verbose);
  vector<size_t> vPermIndices = loadOneColumnFileAsNumbers (permFile, verbose);
  
  computeAndWriteSummaryStatsFtrPerFtr (genoFile,
					phenoFile,
					outFile,
					ftrCoordsFile,
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
					vPermIndices,
					trick,
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
