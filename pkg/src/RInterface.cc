#include "R.h"           // R functions
#include "Rmath.h"       // Rmath
#include "JobHandler.h"  // JobHandler

//-----------------------------------------------------------------------------
//- This file defines TWO functions.
//
//  The first one is a C++ function named 'CppWrapper' that directly accesses
//  the C++ JobHandler class. It's so named because it corresponds to the
// 'main' function of the original stand-alone program that was provided.
//
//  The second one is a C function named 'CWrapper', that INdirectly accesses
//  the C++ JobHandler class, via a call to 'CppWrapper'.

//-----------------------------------------------------------------------------
//- CppWrapper
//
//  C++ function which uses the JobHandler class defined in JobHandler.h and
//  JobHandler.cpp.
//  Since it is a C++ function, it cannot be called directly from R,
//  but instead will be called through an intermediary C function,
//  CppWrapper, defined below.
//
//  'result' is a raw pointer to the result of matrix multiplication.
//  'numRows' and 'numCols' is the number of rows and columns of the first
//  matrix; 'numCols2' is the number of columns of the second matrix.
//  The number of rows of the second matrix must be the number of
//  columns of the first matrix; this will be checked at the level
//  of the R code.
//  'mat1' and 'mat2' are raw pointers to the matrix elements of the two
//   matrices.
//
void CppWrapper(double *matrixX,
                int    *numRows,
                int    *numCols,
                char   **scheme,
                int    *nsteps,
                int    *repeats,
                int    *rankrange,
                int    *numRankRange,
                char   **cltarget,
                char   **clscheme,
                char   **reffile,
                char   **scaling,
                char   **normalizing,
                double *alphas,
                int    *nalphas,
                char   **runtype,
                int    *cstepsize,
                double *idealization,
                double *matrixH,
                double *matrixW,
                double *converged,
                double *totalSteps,
                double *error,
                double *consecutiveError,
                double *Wsparseness,
                double *Hsparseness)
{
//Rprintf("CppWrapper(): entered function...\n"); R_FlushConsole(); R_ProcessEvents();
    // Instantiate JobHandler object.
    JobHandler myJob(matrixX,
                     numRows,
                     numCols,
                     scheme,
                     nsteps,
                     repeats,
                     rankrange,
                     numRankRange,
                     cltarget,
                     clscheme,
                     reffile,
                     scaling,
                     normalizing,
                     alphas,
                     nalphas,
                     runtype,
                     cstepsize,
                     idealization,
                     matrixH,
                     matrixW,
                     converged,
                     totalSteps,
                     error,
                     consecutiveError,
                     Wsparseness,
                     Hsparseness);

    // Perform NMF.
//Rprintf("CppWrapper(): invoking JobHandler::Run()...\n");
//R_FlushConsole();
//R_ProcessEvents();
	myJob.Run();
//Rprintf("CppWrapper(): exiting function...\n");
//R_FlushConsole();
//R_ProcessEvents();

}

//-----------------------------------------------------------------------------
//- CWrapper
//
//  C function that in turn invokes the above C++ function 'CppWrapper'.
//  R can access C code but can't access C++ code directly.
//  This C function provides a C interface to the C++ code that R can access.
//  See: http://www.parashift.com/c++-faq-lite/mixing-c-and-cpp.html
//  In this C function, you must NOT include class declarations,
//  instantiate any C++ objects, or do any oher obvious C++
//  things.
//  The EXTERN statement tells the C++ compiler that the
//  enclosed function 'CWrapper' is a C function.
//  Although apparently we can insert C++ style comments and
//  we can even declare variables in the middle of the function,
//  which I thought you can't do in regular C.
//
extern "C" {
    void CWrapper(double *matrixX,
                  int    *numRows,
                  int    *numCols,
                  char   **scheme,
                  int    *nsteps,
                  int    *repeats,
                  int    *rankrange,
                  int    *numRankRange,
                  char   **cltarget,
                  char   **clscheme,
                  char   **reffile,
                  char   **scaling,
                  char   **normalizing,
                  double *alphas,
                  int    *nalphas,
                  char   **runtype,
                  int    *cstepsize,
                  double *idealization,
                  double *matrixH,
                  double *matrixW,
                  double *converged,
                  double *totalSteps,
                  double *error,
                  double *consecutiveError,
                  double *Wsparseness,
                  double *Hsparseness)
    {
//Rprintf("CWrapper(): entered function...\n"); R_FlushConsole(); R_ProcessEvents();
        CppWrapper(matrixX,
                   numRows,
                   numCols,
                   scheme,
                   nsteps,
                   repeats,
                   rankrange,
                   numRankRange,
                   cltarget,
                   clscheme,
                   reffile,
                   scaling,
                   normalizing,
                   alphas,
                   nalphas,
                   runtype,
                   cstepsize,
                   idealization,
                   matrixH,
                   matrixW,
                   converged,
                   totalSteps,
                   error,
                   consecutiveError,
                   Wsparseness,
                   Hsparseness);
//Rprintf("CWrapper(): exiting function...\n");
//R_FlushConsole();
//R_ProcessEvents();
    }
}
