#ifndef _CLUSTSETNMF_H_
#define _CLUSTSETNMF_H_

class Clust;

#include "clust.h"
#include "DataMatrix.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

// JMM (7/11/2016): Include R header file AFTER system headers. 
// Reference: First few paragraphs in in Section 6 The R API: entry points for C code
// in the online documentation for "Writing R Extensions".
// Also see email from Prof. B. Ripley sent 7/10/2016
// at 3:58 AM entitled "CRAN packages failing to install with modern C++"
// Note: Inclusion of R headers used to be done in the GNMF header files,
// e.g. in 
#define R_NO_REMAP
#include "R.h"           // R functions

class ClustSetNMF : public ClustSet
{
public:
	/*
	ClustSetNMF(const array2d &dataPoints) : dataSet(dataPoints)
	{
	}
	*/

	ClustSetNMF(const DataMatrix &consensus) : dataSet(consensus)
	{
	}

	// An object of class ClustSetNMF is instantiated in Evaluation::IndHC() and in Evaluation::HierarchCustering()
	// CLANG is complaining that an Abstract Class cannot be instantiated.
	// To avoid this problem, do NOT declare GetLeafCount(), GetLeafName(unsigned uNodeIndex),
	// GetLeafId(unsigned uNodeIndex), or ComputeDist() as virtual.
	// No other classes appear to be derived from this class, so this should not cause a problem.
    unsigned GetLeafCount()
	{
		return static_cast<unsigned>( dataSet.nRows ); // (unsigned)()
	}

	const char *GetLeafName(unsigned uNodeIndex)
	{
		Rprintf("ClustSetDF::GetLeafCount (virtual) was called\n");
		R_FlushConsole();
		R_ProcessEvents();
		return 0;
	}

	unsigned GetLeafId(unsigned uNodeIndex)
	{
		return uNodeIndex;
	}
/*
	virtual void JoinNodes(const Clust &C, unsigned uLeftNodeIndex,
	  unsigned uRightNodeIndex, unsigned uJoinedNodeIndex,
	  double *ptrdLeftLength, double *ptrdRightLength)
	{
		Rprintf("ClustSetDF::JoinNodes (virtual), should never be called\n");
		R_FlushConsole();
		R_ProcessEvents();
		return;
	}
*/
	double ComputeDist(const Clust &C, unsigned uNodeIndex1,
	  unsigned uNodeIndex2)
	{
		return (1.0 - dataSet.data[uNodeIndex1][uNodeIndex2]);
	}

private:
	const DataMatrix &dataSet;
};

#endif

