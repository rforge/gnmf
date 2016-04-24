#pragma once

// JMM (12/7/2014): Replace call to C function RAND()
// with call to internal R function UNIF_RAND()
//
// Reference: http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Random-numbers

#include "DataMatrix.h"
#include "CMatrix.h"
#include "ClustDataSet.h"
#include <vector>
#include <math.h>
//
// JMM (2/8/2015): Addressing handling of R headers.
// Enclose in EXTERN "C" statement, do not include RMATH.h header file --
// R.h is apparently all that is necessary.
#define R_NO_REMAP
extern "C" {
    #include <R.h>
}

class FactorMatrix :
	public DataMatrix
{
public:
	array1d modes;
	CMatrix *connectivity;
	CMatrix *realConnet;

	FactorMatrix(void) 
	{
		connectivity = 0;
		realConnet = 0;
	}
	FactorMatrix(const int nRows, const int nColumns, const long _seed) 
		: DataMatrix(nRows, nColumns), seed(_seed)
	{
		connectivity = 0;
		realConnet = 0;
		Initializing();
	}

	FactorMatrix(const int nRows, const int nColumns) 
		: DataMatrix(nRows, nColumns)
	{
		connectivity = 0;
		realConnet = 0;
	}


	FactorMatrix(const DataMatrix & matrix )
		: DataMatrix(matrix)
	{}
	
	// Declare destructor virtual to avoid compiler ambiguity.
	virtual ~FactorMatrix(void) 
	{
		if(connectivity != 0) delete connectivity;
		if(realConnet != 0) delete realConnet;
	}

	virtual void Normalizing( const array1d & modes ) = 0;
	virtual void GetModes() = 0;
	void GetConnectivity();
	void GetConnectivityReal( std::string );
	double GetSparseness();
	int  Multipling(
		const array2d &, 
		array2d &
		) const;

	int  Multipling(
		const array2d &, 
		const array2d &,
		array2d &
		) const;

	virtual int  ParaProducting(
		const array2d &,
		array2d &
		) const = 0;
	virtual int  ParaDividing(
		const array2d &,
		array2d &
		) const = 0;
	virtual int  ParaProducting(
		const array2d &,
		array2d &,
		const array2d &ref
		) const = 0;
	virtual int  ParaDividing( 
		const array2d &,
		array2d &,
		const array2d &ref
		) const = 0;

protected:
	long seed;
	double sparseness;
	array1d sparseValues;

	void Initializing();
	void GetConnectivityVDP();
	void GetConnectivityPearsons();
};
