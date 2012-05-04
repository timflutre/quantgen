/** \file utils.cpp
 *
 *  `utils.cpp' gathers functions useful for any programs.
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
 */

#include <ctime>
#include <cmath>
#include <cstring>
#include <sys/stat.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;

// http://stackoverflow.com/questions/1644868/c-define-macro-for-debug-printing/1644898#1644898
#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif
#define debug_print(fmt, ...)						\
  do { if (DEBUG_TEST) fprintf(stderr, fmt, __VA_ARGS__); } while (0)

/** \brief Split a string.
 *  \note http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
 */
vector<string> & split (const string & s, char delim, vector<string> & tokens)
{
  tokens.clear();
  stringstream ss(s);
  string item;
  while(getline(ss, item, delim)) {
    tokens.push_back(item);
  }
  return tokens;
}

/** \brief Split a string.
 *  \note http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
 */
vector<string> split (const string & s, char delim)
{
  vector<string> tokens;
  return split (s, delim, tokens);
}

/** \brief Return a string with the elapsed time.
 *  \note Output contains only s, or only m and s, or only h, m and s,
 *  or d, h, m and s.
 */
string elapsedTime (time_t startRawTime, time_t endRawTime)
{
  char str[128];
  time_t elapsedSec = static_cast<time_t>(difftime (endRawTime, startRawTime));
  struct tm * ptm;
  ptm = gmtime (&elapsedSec);
  if (ptm->tm_mday == 1 & ptm->tm_hour == 0 & ptm->tm_min == 0)
    sprintf (str, "%is", ptm->tm_sec);
  else if (ptm->tm_mday == 1 & ptm->tm_hour == 0)
    sprintf (str, "%im %is", ptm->tm_min, ptm->tm_sec);
  else if (ptm->tm_mday == 1)
    sprintf (str, "%ih %im %is", ptm->tm_hour, ptm->tm_min, ptm->tm_sec);
  else
    sprintf (str, "%id %ih %im %is", ptm->tm_mday, ptm->tm_hour,
	     ptm->tm_min, ptm->tm_sec);
  return string(str);
}

/** \brief Return a string with the given date-time, without end-of-line.
 */
string time2string (time_t inTime)
{
  char * ptr = ctime (&inTime);
  char buffer[126];
  strcpy (buffer, ptr);
  buffer[strlen(buffer)-1] = 0;
  return string(buffer);
}

/** \brief Load a one-column file.
 */
vector<string> loadOneColumnFile (string inFile, int verbose)
{
  vector<string> vItems;
  
  if (inFile.empty())
    return vItems;
  
  string line;
  ifstream stream;
  vector<string> tokens;
  size_t line_id = 0;
  
  stream.open(inFile.c_str());
  if (! stream.good())
  {
    cerr << "ERROR: can't open file " << inFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout <<"load file " << inFile << "..." << endl;
  }
  
  while (stream.good())
  {
    getline (stream, line);
    if (line.empty())
    {
      break;
    }
    line_id++;
    split (line, '\t', tokens);
    if (tokens.size() != 1)
    {
      cerr << "ERROR: file " << inFile << " should have only one column"
	   << " at line " << line_id << endl;
      exit (1);
    }
    if (find(vItems.begin(), vItems.end(), tokens[0]) == vItems.end())
    {
      vItems.push_back (tokens[0]);
    }
  }
  
  stream.close();
  
  if (verbose > 0)
  {
    cout << "items loaded: " << vItems.size() << endl;
  }
  
  return vItems;
}

/** \brief Load a one-column file into a vector of size_t.
 */
vector<size_t> loadOneColumnFileAsNumbers (string inFile, int verbose)
{
  vector<size_t> vItems;
  
  if (inFile.empty())
    return vItems;
  
  string line;
  ifstream stream;
  vector<string> tokens;
  size_t line_id = 0;
  
  stream.open(inFile.c_str());
  if (! stream.good())
  {
    cerr << "ERROR: can't open file " << inFile << endl;
    exit (1);
  }
  if (verbose > 0)
  {
    cout <<"load file " << inFile << "..." << endl;
  }
  
  while (stream.good())
  {
    getline (stream, line);
    if (line.empty())
    {
      break;
    }
    line_id++;
    split (line, '\t', tokens);
    if (tokens.size() != 1)
    {
      cerr << "ERROR: file " << inFile << " should have only one column"
	   << " at line " << line_id << endl;
      exit (1);
    }
    size_t idx = strtoul (tokens[0].c_str(), NULL, 0);
    if (find(vItems.begin(), vItems.end(), idx) == vItems.end())
    {
      vItems.push_back (idx);
    }
  }
  
  stream.close();
  
  if (verbose > 0)
  {
    cout << "items loaded: " << vItems.size() << endl;
  }
  
  return vItems;
}

/** \brief Used by scandir.
 */
static int dummy_selector (const struct dirent * dir_entry)
{
  return 1;
}

/** \brief Return a vector with the iterations corresponding to nbSteps.
 *  \note Useful with verbose to print at which iteration a loop is.
 */
vector<size_t> getCounters (size_t nbIterations, size_t nbSteps = 5)
{
  vector<size_t> vCounters;
  size_t step = (size_t) floor (nbIterations / nbSteps);
  for (size_t i = 1; i < nbSteps; ++i)
    vCounters.push_back (i * step);
  vCounters.push_back (nbIterations);
  return vCounters;
}

/** \brief Print the nb of iterations already complete in percentage of
 *  the total loop size.
 */
void printCounter (size_t currentIter, vector<size_t> vCounters)
{
  size_t i = 0;
  while (i < vCounters.size())
  {
    if (currentIter == vCounters[i])
    {
      printf ("%.0f%%\n", (float) 100 * currentIter / vCounters[vCounters.size()-1]);
      fflush (stdout);
      break;
    }
    ++i;
  } 
}

/** \brief Convert int, float, etc into a string.
 *  \note http://notfaq.wordpress.com/2006/08/30/c-convert-int-to-string/
 */
template <class T>
inline string toString (const T & t)
{
  stringstream ss;
  ss << t;
  return ss.str();
}

/** \brief Copy a string into another.
 */
string copyString (const string input)
{
  string output;
  for (string::const_iterator it = input.begin();
       it != input.end();
       ++it)
  {
    output += *it;
  }
  return output;
}

/** \brief Replace part of a string with another string.
 *  \note http://stackoverflow.com/a/3418285/597069
 */
void replaceAll (string & str, const string & from, const string & to)
{
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != string::npos)
  {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();  // in case 'to' contains 'from', eg. replacing 'x' with 'yx'
  }
}

/** \brief Round the given value.
 *  \note http://stackoverflow.com/a/485549/597069
 */
double round (double x)
{
  return (x > 0.0) ? floor(x + 0.5) : ceil(x - 0.5);
}

/** Return true if file exists.
 */
bool doesFileExist (const string filename)
{
  bool fexists = false;
  struct stat buffer;
  fexists = ( stat(filename.c_str(), &buffer) == 0);
  return fexists;
}
