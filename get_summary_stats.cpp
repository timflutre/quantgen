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

struct SnpStats
{
  string name; // eg. rs2205177
  string chr; // eg. chr21
  size_t coord; // 1-based coordinate
  double maf; // minor allele frequency
  size_t n; // sample size
  double betahat; // MLE of beta
  double sebetahat; // standard error
  double sigmahat; // MLE of sigma
  double pval; // P-value of the test where H0:"beta=0" and H1:"beta!=0"
  double R2; // coef of determination (proportion of variance explained)
  double rs; // Spearman rank correlation coefficient
  double rsZscore; // using Fisher transformation
  double betaPermPval; // permutation P-value based on beta P-value
  double rsPermPval; // permutation P-value based on Spearman coefs
};

struct FtrStats
{
  string name; // eg. ENSG00000182816
  string chr; // eg. chr21
  size_t start; // 1-based coordinate
  size_t end; // idem
  vector <SnpStats> vSnpStats;
  double betaPermPval; // permutation P-value based on beta P-value
  double rsPermPval; // permutation P-value based on Spearman coef
};

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
       << "  -p, --pheno\tfile with phenotypes (row 1 for sample names,"
       << " column 1" << endl
       << "\t\tfor feature names, delimiter=<space/tab>)" << endl
       << "  -o, --output\toutput file for the summary stats" << endl
       << "      --fcoord\tBED file with the features coordinates" << endl
       << "\t\tfeatures should be in same order than in phenotype file" << endl
       << "  -l, --links\tfile with links between genes and SNPs" << endl
       << "\t\tcustom format: feature<space/tab>SNP|coord" << endl
       << "\t\tfeatures should be in same order than in phenotype file" << endl
       << "\t\tuseful to focus on genetic variants in cis (use windowBed)" << endl
       << "  -c, --chr\tname of the chromosome to analyze (eg. 'chr21')" << endl
       << "  -f, --ftr\tgzipped file with a list of features to analyze" << endl
       << "\t\t(one feature name per line)" << endl
       << "  -s, --snp\tgzipped file with a list of SNPs to analyze" << endl
       << "\t\t(one SNP coordinate per line)" << endl
       << "  -d, --discard\tgzipped file with a list of samples to discard" << endl
       << "\t\t(one individual per line, should match header of phenotype file)" << endl
       << "  -m, --maf\tthreshold for the minor allele frequency (default=0.0)" << endl
       << "\t\t(whatever the option, the MAF will still be computed and saved)" << endl
       << "  -P, --perm\tnumber of phenotype permutations at each feature" << endl
       << "\t\tdefault=0, recommended=10000 (but stop after 100 if P-value > 0.1)" << endl
       << "  -S, --sp\tcompute the Spearman rank correlation coefficient (and Z score)" << endl
       << "\t\tinstead of performing linear regressions" << endl
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
  bool & calcSpearman,
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
	{"maf", required_argument, 0, 'm'},
	{"perm", required_argument, 0, 'P'},
	{"sp", no_argument, 0, 'S'},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:o:l:c:f:s:d:m:P:S",
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
    case 'm':
      minMaf = atof(optarg);
      break;
    case 'P':
      nbPermutations = atol(optarg);
      break;
    case 'S':
      calcSpearman = true;
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
  if (minMaf < 0 || minMaf > 1)
  {
    fprintf (stderr, "ERROR: min MAF should be between 0 and 1 (-m).\n\n");
    help (argv);
    exit (1);
  }
}

void FtrStats_reset (FtrStats * pt_iFtrStats)
{
  pt_iFtrStats->name.clear();
  pt_iFtrStats->chr.clear();
  pt_iFtrStats->start = string::npos;
  pt_iFtrStats->end = string::npos;
  pt_iFtrStats->vSnpStats.clear();
  pt_iFtrStats->betaPermPval = numeric_limits<double>::quiet_NaN();
  pt_iFtrStats->rsPermPval = numeric_limits<double>::quiet_NaN();
}

void
FtrStats_init (
  FtrStats & iFtrStats,
  ifstream & phenoStream,
  ifstream & ftrCoordsStream,
  const string chrToKeep,
  const vector<bool> & vIdxSamplesToSkip,
  const vector<string> vFtrsToKeep,
  vector<bool> & vIsPhenoNa,
  vector<double> & y_init,
  const int verbose)
{
  string linePheno, lineFtrCoords, lineLinks, tok;
  vector<string> tokensPheno, tokensFtrCoords, tokensLinks;
  size_t nbSamples = count (vIdxSamplesToSkip.begin(),
			    vIdxSamplesToSkip.end(),
			    false);
  
  FtrStats_reset (&iFtrStats);
  
  while (true)
  {
    getline (phenoStream, linePheno);
    getline (ftrCoordsStream, lineFtrCoords);
    
    if (linePheno.empty() && lineFtrCoords.empty())
      break;
    
    // both files should have the same number of lines
    if ((linePheno.empty() && ! lineFtrCoords.empty())
	|| (! linePheno.empty() && lineFtrCoords.empty()))
    {
      cerr << "ERROR: phenotype file and features coordinates file" << endl
	   << "should have same number of features" << endl;
      exit (1);
    }
    
    if (linePheno.find('\t') != string::npos)
      split (linePheno, '\t', tokensPheno);
    else
      split (linePheno, ' ', tokensPheno);
    if (lineFtrCoords.find('\t') != string::npos)
      split (lineFtrCoords, '\t', tokensFtrCoords);
    else
      split (lineFtrCoords, ' ', tokensFtrCoords);
    
    // features should be sorted similarly in both files
    if (tokensPheno[0].compare(tokensFtrCoords[3]) != 0)
    {
      cerr << "ERROR: features in phenotype and features coordinates files" << endl
	   << "should be sorted (" << tokensPheno[0]
	   << " versus " << tokensFtrCoords[3] << ")." << endl
	   << "Use the following commands:" << endl
	   << "cat ftr_coords.bed | sort -k4,4 > ftr_coords_sort.bed" << endl
	   << "cat pheno.txt \\" << endl
	   << "| (read -r; printf \"%s\\n\" \"$REPLY\"; sort -k1,1) \\"
	   << endl
	   << "> pheno_sort.txt" << endl;
      exit (1);
    }
    
    // skip it if requested
    if (! vFtrsToKeep.empty()
	&& find(vFtrsToKeep.begin(), vFtrsToKeep.end(), tokensPheno[0])
	== vFtrsToKeep.end())
      break;
    if (! chrToKeep.empty() && chrToKeep.compare(tokensFtrCoords[0]) != 0)
      break;
    
    if (tokensPheno.size()-1 != vIdxSamplesToSkip.size())
    {
      cerr << "ERROR: different number of samples for feature "
	   << tokensPheno[0] << " (" << tokensPheno.size()-1
	   << " vs " << vIdxSamplesToSkip.size() << ")" << endl;
      exit (1);
    }
    
    iFtrStats.name = tokensPheno[0];
    
    // retrieve its mapping information
    iFtrStats.chr = tokensFtrCoords[0];
    iFtrStats.start = atol(tokensFtrCoords[1].c_str()) + 1;
    iFtrStats.end = atol(tokensFtrCoords[2].c_str());
    
    // retrieve its values
    vIsPhenoNa.assign (nbSamples, false);
    y_init.assign (nbSamples, 0);
    size_t j = 0;
    for (size_t colIdx = 1; colIdx < tokensPheno.size(); ++colIdx)
    {
      if(vIdxSamplesToSkip[colIdx-1])
	continue;
      tok = tokensPheno[colIdx];
      if (tok.compare("NA") == 0)
	vIsPhenoNa[j] = true;
      else
	y_init[j] = atof(tok.c_str());
      ++j;
    }
    
    break;
  }
}

void
FtrStats_getCisSnps (
  const FtrStats & iFtrStats,
  ifstream & linksStream,
  map<string, long int> & mSnpNameCoord2Pos,
  vector<string> & vCisSnps)
{
  string line;
  long int linksPos;
  vector<string> tokens;
  bool hasFtrBeenSeen = false;
  
  vCisSnps.clear();
  
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
    if (mSnpNameCoord2Pos.find(tokens[1]) == mSnpNameCoord2Pos.end())
      continue; // SNP to skip (see indexSnps)
    vCisSnps.push_back (tokens[1]);
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
		<< " " << iSnpStats.pval
		<< " " << iSnpStats.R2;
      if (nbPermutations > 0)
	outStream << " " << iSnpStats.betaPermPval
		  << " " << iFtrStats.betaPermPval;
    }
    else
    {
      outStream << " " << iSnpStats.rs
		<< " " << iSnpStats.rsZscore;
      if (nbPermutations > 0)
	outStream << " " << iSnpStats.rsPermPval
		  << " " << iFtrStats.rsPermPval;
    }
    outStream << endl;
  }
}

void
SnpStats_init (
  SnpStats & iSnpStats,
  ifstream & genoStream,
  const string & snpNameCoord,
  map<string, long int> & mSnpNameCoord2Pos,
  const vector<bool> & vIdxSamplesToSkip,
  vector<bool> & vIsGenoNa,
  vector<double> & g_init,
  const int verbose)
{
  size_t i, j;
  double AA, AB, BB, maf;
  string line;
  vector<string> tokens;
  size_t nbSamples = count (vIdxSamplesToSkip.begin(),
			    vIdxSamplesToSkip.end(),
			    false);
  
  genoStream.seekg (mSnpNameCoord2Pos[snpNameCoord]);
  getline (genoStream, line);
  if (line.find('\t') != string::npos)
    split (line, '\t', tokens);
  else
    split (line, ' ', tokens);
  
  // retrieve its mapping information
  iSnpStats.chr = tokens[0];
  iSnpStats.name = tokens[1];
  iSnpStats.coord = atol(tokens[2].c_str());
  
  // retrieve its values
  maf = 0;
  vIsGenoNa.assign (nbSamples, false);
  g_init.assign (nbSamples, 0);
  j = 0;
  for (i = 0; i < vIdxSamplesToSkip.size(); ++i)
  {
    if (vIdxSamplesToSkip[i])
      continue;
    AA = atof(tokens[5+3*i].c_str());
    AB = atof(tokens[5+3*i+1].c_str());
    BB = atof(tokens[5+3*i+2].c_str());
    if (AA == 0 && AB == 0 && BB == 0)
      vIsGenoNa[j] = true;
    else
    {
      g_init[j] = 0 * AA + 1 * AB + 2 * BB;
      maf += g_init[j];
    }
    ++j;
  }
  maf /= 2 * (vIdxSamplesToSkip.size()
	      - count (vIdxSamplesToSkip.begin(),
		       vIdxSamplesToSkip.end(),
		       true)
	      - count (vIsGenoNa.begin(),
		       vIsGenoNa.end(),
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
  const string genoFile,
  const vector<string> & vSnpsToKeep,
  const string chrToKeep,
  const vector<string> & vSamplesToSkip,
  vector<string> & vSamples,
  vector<bool> & vIdxSamplesToSkip,
  const double minMaf,
  map<string, long int> & mSnpNameCoord2Pos,
  int verbose)
{
  ifstream genoStream;
  string line;
  vector<string> tokens;
  long int snpPos;
  size_t totNbSamples, nbSamples;
  
  genoStream.open(genoFile.c_str());
  if (! genoStream.is_open())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  if (verbose > 0)
    cout << "index SNPs in genotype file ..." << endl;
  
  // read header line and record sample names with their column indices
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
  
  // for each SNP
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
    
    // check the format
    if (tokens.size() != 5+3*totNbSamples)
    {
      cerr << line << endl;
      cerr << "ERROR: SNP lines in genotype file should have "
	   << 5+3*totNbSamples << " columns" << endl;
      exit (1);
    }
    
    // skip it if requested
    if (! vSnpsToKeep.empty()
	&& find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[2])
	== vSnpsToKeep.end())
      continue;
    if (! chrToKeep.empty() && chrToKeep.compare(tokens[0]) != 0)
      continue;
    if (minMaf > 0 && getMaf (tokens, vIdxSamplesToSkip) < minMaf)
      continue;
    
    stringstream ss;
    ss << tokens[1] << "|" << tokens[2];  // SNP_name|SNP_coord
    mSnpNameCoord2Pos.insert (make_pair(ss.str(), snpPos));
  }
  
  genoStream.close();
  
  if (verbose > 0)
    cout << "nb of indexed SNPs: " << mSnpNameCoord2Pos.size() << endl;
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
 *  @param yAllSnps vector of vectors of phenotypes for all SNPs
 *  @param gAllSnps vector of vectors of genotypes for all SNPs
 *  @param minBetaPval min P-value for beta=0 over all SNPs
 *
 *  T_i,j: test statistic for SNP i and gene j
 *  T_min,j = min_i (T_i,j)
 *  T~_i,j,k: test statistic for SNP i, gene j and permutation k
 *  T~_min,j,k = min_i (T~_i,j,k)
 *  feature-level Pval = (#{k: T~_min,j,k <= T_min,j} + 1) / (#permutations + 1)
 *  SNP-level Pval = (#{k: T~_min,j,k <= T_i,j} + 1) / (#permutations + 1)
 */
void
computePermutationPvaluesAtFeatureLevel (
  FtrStats & iFtrStats,
  const size_t nbPermutations,
  const vector< vector <double> > yAllSnps,
  const vector< vector <double> > gAllSnps,
  const double minBetaPval,
  const int verbose)
{
  size_t perm_id, snp_id, countNbPerms = 0, seed = 1859;
  double minBetaPvalPerm, betaPvalPerm, betahat, sebetahat, sigmahat, R2;
  vector<double> yPerm, g;
  
  if (verbose > 0)
  {
    printf ("perform %zu permutations on the phenotypes ...\n",
	    nbPermutations);
    fflush (stdout);
  }
  
  // initialize the rng for the permutations (TODO: try multi-threading?)
  gsl_rng_env_setup();
  gsl_rng * r = gsl_rng_alloc (gsl_rng_default);
  if (r == 0)
  {
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit (1);
  }
  gsl_rng_set (r, seed);
  
  iFtrStats.betaPermPval = 1;
  
  for(perm_id=0; perm_id<nbPermutations; ++perm_id)
  {
    ++countNbPerms;
    yPerm = yAllSnps[0];
    gsl_ran_shuffle (r, &yPerm[0], yPerm.size(), sizeof(double));
    minBetaPvalPerm = 1;
    
    for(snp_id=0; snp_id<gAllSnps.size(); ++snp_id)
    {
      g = gAllSnps[snp_id];
      ols ("", "", g, yPerm, &betahat, &sebetahat, &sigmahat,
	   &betaPvalPerm, &R2, 0);
      if (betaPvalPerm < minBetaPvalPerm)
	minBetaPvalPerm = betaPvalPerm;
      if (perm_id == 0)
	iFtrStats.vSnpStats[snp_id].betaPermPval = 1;
      if (betaPvalPerm <= iFtrStats.vSnpStats[snp_id].betaPermPval)
	++(iFtrStats.vSnpStats[snp_id].betaPermPval);
    }
    
    if (minBetaPvalPerm <= minBetaPval)
      ++(iFtrStats.betaPermPval);
    
    // after 100 permutations, see if it's worth doing more of them
    if (countNbPerms == 100)
      if (iFtrStats.betaPermPval / (countNbPerms + 1) > 0.1)
    	break;
  }
  
  // compute the SNP-level P-values
  for(snp_id=0; snp_id<gAllSnps.size(); ++snp_id)
    iFtrStats.vSnpStats[snp_id].betaPermPval /= (countNbPerms + 1);
  
  // compute the feature-level P-value
  iFtrStats.betaPermPval /= (countNbPerms + 1);
  
  gsl_rng_free (r);
}

void
computeSummaryStatsForOneFeature (
  FtrStats & iFtrStats,
  ifstream & genoStream,
  const vector<string> & vCisSnps,
  map<string, long int> & mSnpNameCoord2Pos,
  const vector<bool> vIdxSamplesToSkip,
  const vector<bool> vIsPhenoNa,
  const vector<double> y_init,
  const size_t nbPermutations,
  const bool calcSpearman,
  const int verbose)
{
  size_t snp_id, i;
  string snpNameCoord;
  vector<string> tokens;
  vector<bool> vIsGenoNa;
  vector<double> g_init, y, g;
  vector< vector<double> > yAllSnps, gAllSnps;
  size_t nbSamples = count (vIdxSamplesToSkip.begin(),
			    vIdxSamplesToSkip.end(),
			    false);
  double minBetaPval = 1; // min P-value for beta=0 over all SNPs (no permutation)
  
  // for each SNP in cis
  for (snp_id = 0; snp_id < vCisSnps.size(); ++snp_id)
  {
    SnpStats iSnpStats;
    snpNameCoord = vCisSnps[snp_id];
    SnpStats_init (iSnpStats, genoStream, snpNameCoord, mSnpNameCoord2Pos,
		   vIdxSamplesToSkip, vIsGenoNa, g_init, verbose-1);
    
    if (iSnpStats.name.empty())
      continue;
    if (verbose > 0)
      printf ("SNP %s\n", iSnpStats.name.c_str());
    
    // match phenotypes and genotypes missing values
    y.clear();
    g.clear();
    for (i = 0; i < nbSamples; ++i)
      if (! vIsPhenoNa[i] && ! vIsGenoNa[i])
      {
	y.push_back (y_init[i]);
	g.push_back (g_init[i]);
      }
    iSnpStats.n = y.size();
    yAllSnps.push_back (y);
    gAllSnps.push_back (g);
    
    // perform the ordinary-least-square regression
    if (! calcSpearman)
      ols (iFtrStats.name, iSnpStats.name, g, y, &iSnpStats.betahat,
	   &iSnpStats.sebetahat, &iSnpStats.sigmahat,
	   &iSnpStats.pval, &iSnpStats.R2, verbose-1);
    if (iSnpStats.pval < minBetaPval)
      minBetaPval = iSnpStats.pval;
    
    // or calculate the Spearman coef
    else
    {
      gsl_vector_const_view gsl_g = gsl_vector_const_view_array (&g[0],
								 g.size());
      gsl_vector_const_view gsl_y = gsl_vector_const_view_array (&y[0],
								 y.size());
      iSnpStats.rs = my_stats_correlation_spearman (gsl_g.vector.data, 1,
						    gsl_y.vector.data, 1,
						    g.size());
      iSnpStats.rsZscore = sqrt((g.size() - 3) / 1.06) * 1/2
	* (log(1 + iSnpStats.rs) - log(1 - iSnpStats.rs));
    }
    
    iFtrStats.vSnpStats.push_back (iSnpStats);
  }
  
  // compute permutation P-value at the feature level
  if (nbPermutations > 0)
    computePermutationPvaluesAtFeatureLevel (iFtrStats,
					     nbPermutations,
					     yAllSnps,
					     gAllSnps,
					     minBetaPval,
					     verbose);
}

void
computeAndWriteSummaryStatsFtrPerFtr (
  const string genoFile,
  const string phenoFile,
  const string outFile,
  const string ftrCoordsFile,
  const string linksFile,
  const vector<string> vFtrsToKeep,
  const vector<string> vSnpsToKeep,
  const vector<string> vSamplesToSkip,
  const string chrToKeep,
  const double minMaf,
  const size_t nbPermutations,
  const bool calcSpearman,
  const int verbose)
{
  FtrStats iFtrStats;
  string linePheno, lineLinks;
  size_t nbFtrs = 0, nbAnalyzedPairs = 0, nbAnalyzedFtrs = 0;
  vector<size_t> vCounters;
  vector<string> tokens, vSamples, vCisSnps;
  vector<double> y_init;
  vector<bool> vIdxSamplesToSkip, vIsPhenoNa;
  map<string, long int> mSnpNameCoord2Pos;
  ifstream phenoStream, genoStream, ftrCoordsStream, linksStream;
  ofstream outStream;
  
  // open input files
  phenoStream.open(phenoFile.c_str());
  if (! phenoStream.good())
  {
    cerr << "ERROR: can't open file " << phenoFile << endl;
    exit (1);
  }
  genoStream.open(genoFile.c_str());
  if (! genoStream.good())
  {
    cerr << "ERROR: can't open file " << genoFile << endl;
    exit (1);
  }
  linksStream.open(linksFile.c_str());
  if (! linksStream.good())
  {
    cerr << "ERROR: can't open file " << linksFile << endl;
    exit (1);
  } 
  ftrCoordsStream.open(ftrCoordsFile.c_str());
  if (! ftrCoordsStream.good())
  {
    cerr << "ERROR: can't open file " << ftrCoordsFile << endl;
    exit (1);
  }
  
  // open output file and write the header line
  outStream.open(outFile.c_str());
  if (! outStream.good())
  {
    cerr << "ERROR: can't open file " << outFile << endl;
    exit (1);
  }
  outStream << "ftr chr start end snp chr coord maf n";
  if (! calcSpearman)
  {
    outStream << " betahat sebetahat sigmahat pvalBeta R2";
    if (nbPermutations > 0)
      outStream << " betaPermPvalSnp betaPermPvalFtr";
  }
  else
  {
    outStream << " rs rsZscore";
    if (nbPermutations > 0)
      outStream << " rsPermPvalSnp rsPermPvalFtr";
  }
  outStream << endl;
  
  // index the genotype file
  indexSnps (genoFile, vSnpsToKeep, chrToKeep, vSamplesToSkip, vSamples,
	     vIdxSamplesToSkip, minMaf, mSnpNameCoord2Pos, verbose);
  
  if (verbose > 0)
    printf ("look for association between each pair feature-SNP ...\n");
  
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
  
  // for each feature
  while (phenoStream.good())
  {
    // initialize it
    FtrStats_init (iFtrStats, phenoStream, ftrCoordsStream, chrToKeep,
		   vIdxSamplesToSkip, vFtrsToKeep, vIsPhenoNa, y_init,
		   verbose-1);
    if (iFtrStats.name.empty())
      continue;
    ++nbFtrs;
    
    // retrieve its SNPs in cis
    FtrStats_getCisSnps (iFtrStats, linksStream, mSnpNameCoord2Pos, vCisSnps);
    if (vCisSnps.empty())
      continue;
    
    if (verbose > 1)
    {
      printf ("analyzing feature %s (%s, %zu cis SNPs) ...\n",
	      iFtrStats.name.c_str(), iFtrStats.chr.c_str(),
	      vCisSnps.size());
      fflush (stdout);
    }
    ++nbAnalyzedFtrs;
    
    // loop over SNPs in cis
    computeSummaryStatsForOneFeature (iFtrStats, genoStream, vCisSnps,
				      mSnpNameCoord2Pos, vIdxSamplesToSkip,
				      vIsPhenoNa, y_init, nbPermutations,
				      calcSpearman, verbose-2);
    
    // write the results (one line per SNP)
    if (iFtrStats.vSnpStats.size() > 0)
    {
      FtrStats_write (iFtrStats, outStream, nbPermutations, calcSpearman);
      nbAnalyzedPairs += iFtrStats.vSnpStats.size();
    }
  }
  
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

int main (int argc, char ** argv)
{
  string genoFile, phenoFile, outFile, ftrCoordsFile, linksFile,
    ftrsFile, snpsFile, samplesFile, chrToKeep = "";
  double minMaf = 0.0;
  size_t nbPermutations = 0;
  bool calcSpearman = false;
  int verbose = 1;
  parse_args (argc, argv, genoFile, phenoFile, outFile, ftrCoordsFile,
	      linksFile, chrToKeep, ftrsFile, snpsFile, samplesFile,
	      minMaf, nbPermutations, calcSpearman, verbose);
  
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
  
  computeAndWriteSummaryStatsFtrPerFtr (genoFile,
					phenoFile,
					outFile,
					ftrCoordsFile,
					linksFile,
					vFtrsToKeep,
					vSnpsToKeep,
					vSamplesToSkip,
					chrToKeep,
					minMaf,
					nbPermutations,
					calcSpearman,
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
