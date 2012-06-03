#ifndef _MUSCLE_H_
#define _MUSCLE_H_

#include <ctype.h>
#include <stdio.h>


// insane value for uninitialized variables
const unsigned uInsane = 8888888;
const int iInsane = 8888888;
const char cInsane = (char) 0xcd;		// int 3 instruction, used e.g. for unint. memory
const double dInsane = -9e29;
const double fInsane = (double) -9e29;
const double PLUS_INFINITY = 1e+37f;
const char INVALID_STATE = '*';
const double g_dSUEFF = (double) 0.1;

// extern double g_dNAN;
extern unsigned long g_tStart;

#endif

