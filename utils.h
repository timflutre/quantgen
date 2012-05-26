/** \file utils.h
 *
 *  `utils.cpp' gathers functions useful for any programs.
 *  Copyright (C) 2012  T. Flutre
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

#ifndef UTILS_H
#define UTILS_H

#include <ctime>

#include <vector>
#include <string>
using namespace std;

vector<string> & split (const string & s, char delim, vector<string> & tokens);

vector<string> split (const string & s, char delim);

string elapsedTime (time_t startRawTime, time_t endRawTime);

string time2string (time_t inTime);

vector<string> loadOneColumnFile (string inFile, int verbose);

vector<size_t> loadOneColumnFileAsNumbers (string inFile, int verbose);

vector<size_t> getCounters (size_t nbIterations, size_t nbSteps);

void printCounter (size_t currentIter, vector<size_t> vCounters);

template <class T> inline string toString (const T & t);

string copyString (const string input);

void replaceAll (string & str, const string & from, const string & to);

double round (double x);

bool doesFileExist (const string filename);

vector<string> scanInputDirectory (const string & inDir, const int & verbose);

bool isDirectory(const char path[]);

int removeDir(string path);

size_t getSeed (void);

#endif
