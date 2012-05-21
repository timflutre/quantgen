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
#include <dirent.h>
#include <cerrno>

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

/** \brief Return a string with the elapsed time in d, h, m and s.
 *  \note http://stackoverflow.com/a/2419597/597069
 */
string elapsedTime (time_t startRawTime, time_t endRawTime)
{
  char str[128];
  double elapsed = difftime (endRawTime, startRawTime); // in sec
  sprintf (str, "%01.0fd %01.0fh %01.0fm %01.0fs", floor(elapsed/(24*60*60)),
	   floor(elapsed/(60*60)), floor(fmod(elapsed,60*60)/60.0),
	   fmod(elapsed,60));
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
    cout <<"load file " << inFile << " ..." << endl;
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
    cout <<"load file " << inFile << " ..." << endl;
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

/** \brief List a given directory.
 */
vector<string>
scanInputDirectory (
  const string & inDir,
  const int & verbose)
{
  vector<string> vInFiles;
  struct dirent ** inFiles = NULL;
  int nbInFiles;
  if (verbose > 0)
  {
    cout << "scan directory " << inDir << " ..." << endl;
  }
  nbInFiles = scandir(inDir.c_str(), &inFiles, dummy_selector, alphasort);
  if (nbInFiles == -1)
  {
    cerr << "ERROR: can't scan " << inDir << endl;
    exit (1);
  }
  else if (nbInFiles == 0)
  {
    cerr << "ERROR: " << inDir << " contains no file" << endl;
    exit (1);
  }
  else
  {
    for (int s = 0; s < nbInFiles; ++s)
    {
      if (string(inFiles[s]->d_name) == "." ||
	  string(inFiles[s]->d_name) == "..")
      {
	free (inFiles[s]);
	continue;
      }
      char path[1024];
      int nbChar;
      if (inDir[inDir.size()-1] != '/')
	nbChar = sprintf (path, "%s/%s", inDir.c_str(), inFiles[s]->d_name);
      else
	nbChar = sprintf (path, "%s%s", inDir.c_str(), inFiles[s]->d_name);
      if (nbChar < 0)
      {
	cerr << "ERROR: variable 'path' is not big enough" << endl;
      }
      vInFiles.push_back (string(path));
      free (inFiles[s]);
    }
    if (verbose > 0)
      cout << "nb of files: " << vInFiles.size() << endl;
  }
  free (inFiles);
  return vInFiles;
}

/** \brief Return true if the given path is a directory.
 *  \note http://stackoverflow.com/a/1149769/597069
 */
bool isDirectory(char path[]) {
  bool res = false;
  if (path[strlen(path)] == '.') // exception for \. and \..
    res = true;
  else
  {
    for(int i = strlen(path) - 1; i >= 0; i--)
    {
      if (path[i] == '.')
      {
	res = false; // if we first encounter a . then it's a file
	break;
      }
      else if (path[i] == '\\' || path[i] == '/')
      {
	res = true; // if we first encounter a \ it's a dir
	break;
      }
    }
  }
  return res;
}

/** \brief Remove a directory even if it is not empty.
 *  \note http://stackoverflow.com/a/1149769/597069
 *  \note Don't do anything if the supplied path is empty
 *  or if the directory doesn't exist.
 */
int removeDir(string path) {
  if (path.empty())
    return 0;
  
  if (path[path.size()] == '.')
    return 0;
  
  if (path[path.length()-1] != '/')
    path += "/";
  
  // create a pointer to a directory
  DIR *pdir = NULL;
  pdir = opendir (path.c_str());
  if (pdir == NULL)
  {
    if (errno == 2) // No such file or directory
      return 0;
    else
    {
      cerr << "ERROR: opendir returned NULL for path " << path << endl;
      fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
      return errno;
    }
  }
  
  struct dirent *pent = NULL;
  char file[1024];
  int counter = 1; // use this to skip the first TWO which cause an infinite loop (and eventually, stack overflow)
  while (true)
  {
    pent = readdir (pdir); // while there is still something in the directory
    if (pent == NULL)
    {
      if (errno != 0) // if pent has not been initialised correctly
      {
	cerr << "ERROR: readdir returned NULL for path " << path << endl;
	fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
	return errno; // we couldn't do it
      }
      else // if the directory is empty
	break;
    }
    if (counter > 2)
    {
      for (int i = 0; i < 256; i++)
	file[i] = '\0';
      strcat(file, path.c_str());
      // otherwise, it was initialised correctly, so let's delete the file~
      strcat(file, pent->d_name); // concatenate the strings to get the complete path
      if (isDirectory(file) == true)
	removeDir(file);
      else // it's a file, we can use remove
	remove(file);
    }
    counter++;
  }
  
  // finally, let's clean up
  closedir (pdir); // close the directory
  if (rmdir(path.c_str()) != 0)
  {
    if (errno != 0)
    {
      cerr << "ERROR: rmdir returned an error" << endl;
      fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
      return errno;
    }
  }
  
  return 0;
}
