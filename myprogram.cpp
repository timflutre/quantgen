/** \file myprogram.cpp
 *
 *  `myprogram' does this and that.
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
 *  g++ -Wall -g utils_io.cpp myprogram.cpp -lgsl -lgslcblas -lz -o myprogram
 */

#include <cmath>
#include <ctime>
#include <getopt.h>

#include <iostream>
#include <string>
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
       << " does this and that." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "      --in\tinput" << endl
       << endl
       << "Examples:" << endl
       << "  " << argv[0] << " --in <input>" << endl
       << endl
       << "Remarks:" << endl
       << "  This is my typical template file for C++." << endl
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
void
parseCmdLine(
  int argc,
  char ** argv,
  string & input,
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
      {"in", required_argument, 0, 0},
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
      if(strcmp(long_options[option_index].name, "in") == 0)
      {
        input = optarg;
        break;
      }
    case 'h':
      help (argv);
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
  if(input.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --in" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! doesFileExist(input)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << input << endl << endl;
    help(argv);
    exit(1);
  }
}

int main(int argc, char ** argv)
{
  string input;
  int verbose = 1;
  
  parseCmdLine(argc, argv, input, verbose);
  
  time_t startRawTime, endRawTime;
  if(verbose > 0)
  {
    time(&startRawTime);
    cout << "START " << basename(argv[0])
         << " " << getDateTime(startRawTime) << endl
         << "version " << VERSION << " compiled " << __DATE__
         << " " << __TIME__ << endl
         << "cmd-line: " << getCmdLine(argc, argv) << endl
         << "cwd: " << getCurrentDirectory() << endl;
    cout << flush;
  }
  
  // ... specific code ...
  
  if(verbose > 0)
  {
    time(&endRawTime);
    cout << "END " << basename(argv[0])
         << " " << getDateTime(endRawTime) << endl
         << "elapsed -> " << getElapsedTime(startRawTime, endRawTime) << endl
         << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
  }
  
  return EXIT_SUCCESS;
}
