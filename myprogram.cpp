/** \file myprogram.cpp
 *
 *  `myprogram' does this and that.
 * choose between:
 * Author: Timothée Flutre
 * Not copyrighted -- provided to the public domain
 * or:
 *  Copyright (C) 2011-2014 Timothée Flutre
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
 *  Compile with: g++ -Wall -g utils_io.cpp myprogram.cpp -lgsl -lgslcblas -lz -o myprogram
 *  "-lgsl -lgslcblas" are just provided as example
 */

#include <cmath>
#include <ctime>
#include <cstring>
#include <getopt.h>
#include <libgen.h>

#include <iostream>
#include <string>
using namespace std;

#include "utils_io.hpp"
using namespace utils;

#ifndef VERSION
#define VERSION "1.0"
#endif

/** \brief Display the help on stdout.
 *  \note The format complies with help2man (http://www.gnu.org/s/help2man)
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
       << "  -i, --input\tpath to the input file" << endl
       << endl
       << "Examples:" << endl
       << "  " << argv[0] << " -i <input>" << endl
       << endl
       << "Report bugs to <>." << endl
    ;
}

/** \brief Display version and license information on stdout.
 */
void version(char ** argv)
{
  cout << argv[0] << " " << VERSION << endl
       << endl
       << "Copyright (C) 2011-2014 Timothée Flutre." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by Timothée Flutre." << endl
    ;
  // or choose "Not copyrighted -- provided to the public domain"
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
      {"input", required_argument, 0, 0},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "hVv:i:",
                    long_options, &option_index);
    if(c == -1)
      break;
    switch(c)
    {
    case 0:
      if(long_options[option_index].flag != 0)
        break;
      if(strcmp(long_options[option_index].name, "input") == 0)
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
    case 'i':
      input = optarg;
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
	 << "ERROR: missing compulsory option --input" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! doesFileExist(input)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find file " << input << endl << endl;
    help(argv);
    exit(1);
  }
}

void run(const string & input, const int & verbose)
{
  
  // specific code ...
  
}

int main(int argc, char ** argv)
{
  string input;
  int verbose = 1;
  
  parseCmdLine(argc, argv, input, verbose);
  
  time_t startRawTime, endRawTime;
  if(verbose > 0){
    time(&startRawTime);
    cout << "START " << basename(argv[0])
         << " " << getDateTime(startRawTime) << endl
         << "version " << VERSION << " compiled " << __DATE__
         << " " << __TIME__ << endl
         << "cmd-line: " << getCmdLine(argc, argv) << endl
         << "cwd: " << getCurrentDirectory() << endl;
    cout << flush;
  }
  
  run(input, verbose);
  
  if(verbose > 0){
    time(&endRawTime);
    cout << "END " << basename(argv[0])
         << " " << getDateTime(endRawTime) << endl
         << "elapsed -> " << getElapsedTime(startRawTime, endRawTime) << endl
         << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
  }
  
  return EXIT_SUCCESS;
}
