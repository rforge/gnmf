#ifndef _CLUSTSETNMF_H_
#define _CLUSTSETNMF_H_

class Clust;

#include "clust.h"
#include "DataMatrix.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
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

	virtual unsigned GetLeafCount()
	{
		return static_cast<unsigned>( dataSet.nRows ); // (unsigned)()
	}

	virtual const char *GetLeafName(unsigned uNodeIndex)
	{
Rprintf("ClustSetDF::GetLeafCount\n");
R_FlushConsole();
R_ProcessEvents();
return 0;
//		std::ostringstream out;
//		out << dataSet.rowNames[uNodeIndex];
//		return out.str().c_str();
	}

	virtual unsigned GetLeafId(unsigned uNodeIndex)
	{
		return uNodeIndex;
	}

	virtual void JoinNodes(const Clust &C, unsigned uLeftNodeIndex,
	  unsigned uRightNodeIndex, unsigned uJoinedNodeIndex,
	  double *ptrdLeftLength, double *ptrdRightLength)
	{
Rprintf("ClustSetDF::JoinNodes, should never be called\n");
R_FlushConsole();
R_ProcessEvents();
return;
//		std::cout << "ClustSetDF::JoinNodes, should never be called" << std::endl;
//		exit (-1);
	}

	virtual double ComputeDist(const Clust &C, unsigned uNodeIndex1,
	  unsigned uNodeIndex2)
	{
		return (1.0 - dataSet.data[uNodeIndex1][uNodeIndex2]);
	}

private:
	const DataMatrix &dataSet;
};

#endif

