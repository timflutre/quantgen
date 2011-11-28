/** \file impute2bimbam.cpp
 *
 *  `impute2bimbam' converts genotype data format from IMPUTE into BIMBAM.
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
 *  gcc -Wall impute2bimbam.cpp -lstdc++ -lgzstream -lz -o impute2bimbam
 *  help2man -o impute2bimbam.man ./impute2bimbam
 *  groff -mandoc impute2bimbam.man > impute2bimbam.ps
*/

#include <cmath>
#include <ctime>
#include <getopt.h>

#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <fstream>
using namespace std;

#include "utils.cpp"

//-----------------------------------------------------------------------------

void help (char ** argv);
void version (char ** argv);
void parse_args (int argc, char ** argv,
		 string * pt_input,
		 string * pt_output,
		 int * pt_verbose);
void convertImputeFileToBimbamFile (string input,
				    string output,
				    int verbose);
int main (int argc, char ** argv);

//-----------------------------------------------------------------------------

/** \brief Display the usage on stdout.
*/
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " converts genotype data format from IMPUTE into BIMBAM." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS]..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (default=1)" << endl
       << "  -i, --input\tpath to input directory and generic file name in the IMPUTE format" << endl
       << "\t\teg. '~/data/chrXX_chunkAll.impute2'" << endl
       << "  -o, --output\tgeneric name for the output files" << endl
       << "\t\tdefault is 'chrXX.bimbam'" << endl
       << endl
       << "Examples:" << endl
       << "  " << argv[0] << " -i ~/data/chrXX_chunkAll.impute2" << endl;
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
		 string * pt_input,
		 string * pt_output,
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
	{"output", required_argument, 0, 'o'},
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:i:o:",
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
      *pt_input = optarg;
      break;
    case 'o':
      *pt_output = optarg;
      break;
    case '?':
      break;
    default:
      abort ();
    }
  }
  if ((*pt_input).empty())
  {
    fprintf (stderr, "ERROR: missing input (-i).\n");
    help (argv);
    exit (1);
  }
  if ((*pt_output).empty())
  {
    *pt_output = "chrXX.bimbam";
  }
}

void convertImputeFilesToBimbamFiles (string input,
				      string output,
				      int verbose)
{
  string line;
  ifstream inStream;
  vector<string> tokens;
  ofstream outStream;
  size_t nbSamples = 0;
  
  for (int chrNb = 1; chrNb <= 22; ++chrNb)
  {
    if (verbose > 0)
    {
      cout << "convert genotypes on chr" << chrNb << "..." << endl;
      fflush (stdout);
    }
    
    string inFile = copyString (input);
    replaceAll (inFile, "XX", toString(chrNb));
    inStream.open (inFile.c_str());
    if (! inStream.is_open())
    {
      cerr << "ERROR: can't open file " << inFile << endl;
      exit (1);
    }
    
    string outFile = copyString (output);
    replaceAll (outFile, "XX", toString(chrNb));
    outStream.open (outFile.c_str());
    if (! outStream.is_open())
    {
      cerr << "ERROR: can't open file " << outFile << endl;
      exit (1);
    }
    
    while (inStream.good())
    {
      getline (inStream, line);
      if (line.empty())
      {
	break;
      }
      if (line.find('\t') != string::npos)
	split (line, '\t', tokens);
      else
	split (line, ' ', tokens);
      outStream << tokens[1]          // SNP id
		<< " " << tokens[3]   // allele A (minor allele for BimBam)
		<< " " << tokens[4];  // allele B (major allele for BimBam)
      nbSamples = (size_t) floor ((tokens.size() - 5) / 3);
      for (size_t i = 0; i < nbSamples; ++i)
      {
	outStream << " " << 2 * atof(tokens[5+3*i].c_str())
	  + 1 * atof(tokens[5+3*i+1].c_str())
	  + 0 * atof(tokens[5+3*i+2].c_str());
      }
      outStream << endl;
    }
    
    inStream.close();
    outStream.close();
  }
}

int main (int argc, char ** argv)
{
  string input, output;
  int verbose = 1;
  parse_args (argc, argv, &input, &output, &verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")" << endl;
  }
  
  convertImputeFilesToBimbamFiles (input, output, verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime)
	 << ": elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << ")" << endl;
  }
  
  return 0;
}
