#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>

const  int NAMEL = 200;
const  int CPU0 = 1;                  // start CPU
const  int MAXLINEWIDTH = 20000;      // Maximum character number in a line to be gotten with getline()
// const  double MAXVALUE = 1.0e+15;     // for extracting the minimum value
// const  double MINVALUE = -1.0e+15;    // for extracting the maximum value
const double ACCEPTEDDIFF = 1.0e-15;  // the minimum meaningful difference between D and Mock element
const double EPS = 1.0e-15;           // precision control
const double UPPERBOUND = 50.0;       // upper limit cut for exponential power
const double CONVERGETHRESHOLD = 0.01; // V/M is [0.98, 1.02]
// const int CONVERGECHECK = 20; // every 50 steps check convergence
const double OUTLINERTHRESHOLD = 2.5;

// default parameter
const int nSteps_ = 2000;
const int nChains_ = 20;
const std::string scheme_ = "KL";
const std::string target_ = "PATTERN (H mtrix)";
const int nAlphas_ = 1;
const bool Normalization_ = false;
const bool InitialScaling_ = false;
const bool ConvergeTest_ = true;
const std::string clustScheme_ = "Binary";
const std::string RunType_ = "whole"; // could be "simulation" or "evaluation"
const long SEED = 12458698;
const double THETA = 0.0;

typedef std::vector< std::vector< double > > array2d;
typedef std::vector< double > array1d;
typedef std::map<int, int> IIMap;
typedef std::map<std::string, std::string> SSMap;
typedef std::map<std::string, double> SDMap;

#endif

