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
       << "  -o, --output\tgeneric prefix for the output files" << endl
       << "\t\tdefault is 'chrXX' leading to 'chrXX.bimbam' and 'chrXX_snpAnnot.txt'" << endl
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
    fprintf (stderr, "ERROR: missing input (-i).\n\n");
    help (argv);
    exit (1);
  }
  if ((*pt_output).empty())
  {
    *pt_output = "chrXX";
  }
}

void convertImputeFilesToBimbamFiles (string input,
				      string output,
				      int verbose)
{
  string line;
  ifstream inStream;
  vector<string> tokens;
  ofstream outStream1, outStream2;
  size_t nbSamples = 0;
  stringstream ss;
  
  for (int chrNb = 1; chrNb <= 22; ++chrNb)
  {
    string inFile = copyString (input);
    replaceAll (inFile, "XX", toString(chrNb));
    if (! doesFileExist (inFile))
    {
      continue;
    }
    if (verbose > 0)
    {
      cout << "convert genotypes from file '" << inFile << "' ..." << endl;
      fflush (stdout);
    }
    
    inStream.open (inFile.c_str());
    if (! inStream.is_open())
    {
      cerr << "ERROR: can't open file " << inFile << endl;
      exit (1);
    }
    
    ss.clear();
    ss.str(string());  // http://stackoverflow.com/a/834631/597069
    ss << output << ".bimbam";
    string outFile1 = ss.str();
    replaceAll (outFile1, "XX", toString(chrNb));
    outStream1.open (outFile1.c_str());
    if (! outStream1.is_open())
    {
      cerr << "ERROR: can't open file " << outFile1 << endl;
      exit (1);
    }
    ss.clear();
    ss.str(string());
    ss << output << "_snpAnnot.txt";
    string outFile2 = ss.str();
    replaceAll (outFile2, "XX", toString(chrNb));
    outStream2.open (outFile2.c_str());
    if (! outStream2.is_open())
    {
      cerr << "ERROR: can't open file " << outFile2 << endl;
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
      outStream1 << tokens[1]          // SNP id
		<< " " << tokens[3]    // allele A (minor allele for BimBam)
		<< " " << tokens[4];   // allele B (major allele for BimBam)
      nbSamples = (size_t) floor ((tokens.size() - 5) / 3);
      for (size_t i = 0; i < nbSamples; ++i)
      {
	outStream1 << " " << 2 * atof(tokens[5+3*i].c_str())
	  + 1 * atof(tokens[5+3*i+1].c_str())
	  + 0 * atof(tokens[5+3*i+2].c_str());
      }
      outStream1 << endl;
      outStream2 << tokens[1]          // SNP id
		 << " " << tokens[2]   // SNP coordinate
		 << " " << chrNb
		 << endl;
    }
    
    inStream.close();
    outStream1.close();
    outStream2.close();
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
