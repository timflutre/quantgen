/** \file extract_bed_from_names.cpp
 *
 *  `extract_bed_from_names' extracts a list of BED records based on their names
 *  Copyright (C) 2011-2013 Timothee Flutre
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
 *  g++ -Wall -g utils_io.cpp extract_bed_from_names.cpp -lgsl -lgslcblas -lz -o extract_bed_from_names
 */

#include <cmath>
#include <ctime>
#include <getopt.h>

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "utils_io.hpp"
using namespace utils;

#ifndef VERSION
#define VERSION "1.0"
#endif

/** \brief Display the help on stdout.
 */
void help(char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " extracts a list of BED records based on their names." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "      --names\tfile with one record name per line" << endl
       << "      --in\tinput BED file" << endl
       << "      --out\toutput BED file (gzipped)" << endl
    ;
}
/** \brief Display version and license information on stdout.
 */
void version(char ** argv)
{
  cout << argv[0] << " " << VERSION << endl
       << endl
       << "Copyright (C) 2011-2013 Timothee Flutre." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by Timothee Flutre." << endl;
}

/** \brief Parse the command-line arguments and check the values of the 
 *  compulsory ones.
 */
void parseCmdLine(
  int argc,
  char ** argv,
  string & namesFile,
  string & inBedFile,
  string & outBedFile,
  int & verbose)
{
  int c = 0;
  while(true)
  {
    static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"verbose", required_argument, 0, 'v'},
      {"names", required_argument, 0, 0},
      {"in", required_argument, 0, 0},
      {"out", required_argument, 0, 0},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "hVv:",
                    long_options, &option_index);
    if(c == -1)
      break;
    switch(c)
    {
    case 0:
      if(long_options[option_index].flag != 0)
        break;
      if(strcmp(long_options[option_index].name, "names") == 0)
      {
        namesFile = optarg;
        break;
      }
      if(strcmp(long_options[option_index].name, "in") == 0)
      {
        inBedFile = optarg;
        break;
      }
      if(strcmp(long_options[option_index].name, "out") == 0)
      {
        outBedFile = optarg;
        break;
      }
    case 'h':
      help(argv);
      exit(0);
    case 'V':
      version(argv);
      exit(0);
    case 'v':
      verbose = atoi(optarg);
      break;
    case '?':
      printf("\n"); help(argv);
      abort();
    default:
      printf("\n"); help(argv);
      abort();
    }
  }
  if(namesFile.empty())
  {
    getCmdLine(argc, argv);
    fprintf(stderr, "ERROR: missing compulsory option --names\n\n");
    help(argv);
    exit(1);
  }
  if(! doesFileExist(namesFile))
  {
    getCmdLine(argc, argv);
    fprintf(stderr, "ERROR: can't find '%s'\n\n", namesFile.c_str());
    help(argv);
    exit(1);
  }
  if(inBedFile.empty())
  {
    getCmdLine(argc, argv);
    fprintf(stderr, "ERROR: missing compulsory option --in\n\n");
    help(argv);
    exit(1);
  }
  if(! doesFileExist(inBedFile))
  {
    getCmdLine(argc, argv);
    fprintf(stderr, "ERROR: can't find '%s'\n\n", inBedFile.c_str());
    help(argv);
    exit(1);
  }
  if(outBedFile.empty())
  {
    getCmdLine(argc, argv);
    fprintf(stderr, "ERROR: missing compulsory option --out\n\n");
    help(argv);
    exit(1);
  }
}

void loadNames(
  const string & namesFile,
  const int & verbose,
  vector<string> & names)
{
  if(verbose > 0)
    cout << "load names from file " << namesFile << " ..." << endl;
  
  gzFile stream;
  vector<string> tokens;
  string line;
  openFile(namesFile, stream, "rb");
  while(getline(stream, line)){
    split(line, " \t", tokens);
    if(find(names.begin(), names.end(), tokens[0]) == names.end())
      names.push_back(tokens[0]);
  }
  if(! gzeof(stream)){
    cerr << "ERROR: can't read successfully file "
	 << namesFile << " up to the end" << endl;
    exit(1);
  }
  closeFile(namesFile, stream);
  
  if (verbose > 0)
    cout << "nb of names: " << names.size() << endl;
}

void extractBedRecords(
  const vector<string> & names,
  const string & inBedFile,
  const string & outBedFile,
  const int & verbose)
{
  if(verbose > 0)
    cout << "extract records from file " << inBedFile << " ..." << endl;
  
  gzFile inStream, outStream;
  vector<string> tokens;
  string line;
  stringstream txt;
  size_t nb_lines_out = 0;
  openFile(inBedFile, inStream, "rb");
  openFile(outBedFile, outStream, "wb");
  while(getline(inStream, line)){
    split(line, " \t", tokens);
    if(find(names.begin(), names.end(), tokens[3]) != names.end()){
      txt.str("");
      txt << tokens[0];
      for(size_t i = 1; i < tokens.size(); ++i)
	txt << "\t" << tokens[i];
      txt << endl;
      ++nb_lines_out;
      gzwriteLine(outStream, txt.str(), outBedFile, nb_lines_out);
    }
  }
  if(! gzeof(inStream)){
    cerr << "ERROR: can't read successfully file "
	 << inBedFile << " up to the end" << endl;
    exit(1);
  }
  closeFile(inBedFile, inStream);
  closeFile(outBedFile, outStream);
}

void run(
  const string & namesFile,
  const string & inBedFile,
  const string & outBedFile,
  const int & verbose)
{
  vector<string> names;
  loadNames(namesFile, verbose, names);
  
  extractBedRecords(names, inBedFile, outBedFile, verbose);
}

int main(int argc, char ** argv)
{
  string namesFile, inBedFile, outBedFile;
  int verbose = 1;
  
  parseCmdLine(argc, argv, namesFile, inBedFile, outBedFile, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << basename(argv[0])
         << " " << getDateTime (startRawTime) << endl
         << "version " << VERSION << " compiled " << __DATE__
         << " " << __TIME__ << endl
         << "cmd-line: " << getCmdLine (argc, argv) << endl
         << "cwd: " << getCurrentDirectory() << endl;
    cout << flush;
  }
  
  run(namesFile, inBedFile, outBedFile, verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << basename(argv[0])
         << " " << getDateTime (endRawTime) << endl
         << "elapsed -> " << getElapsedTime(startRawTime, endRawTime) << endl
         << "max.mem -> " << getMaxMemUsedByProcess2Str () << endl;
  }
  
  return EXIT_SUCCESS;
}
