#pragma once
#include "FactorMatrix.h"

class WMatrix :
	public FactorMatrix
{
	// friend std::istream &operator >> (std::istream &in, WMatrix &myMatrix);
	// friend std::ostream &operator << (std::ostream &out, const WMatrix &myMatrix);

public:
	WMatrix(void) {}

	WMatrix(const int nRows, const int nColumns, const long seed) 
		: FactorMatrix(nRows, nColumns, seed)
	{}

	WMatrix(const int nRows, const int nColumns) 
		: FactorMatrix(nRows, nColumns)
	{}

	WMatrix( const DataMatrix & _matrix) 
		: FactorMatrix( _matrix )
	{}

	// Declare destructor virtual to avoid compiler ambiguity.
	virtual ~WMatrix(void) {
		if(connectivity != 0) delete connectivity;
		if(realConnet != 0) delete realConnet;
	}

	void Normalizing( const array1d & );
	void GetModes();
	void GetModes(const double bdn, const std::vector<double> &v);
	int  ParaProducting(
		const array2d &,
		array2d &
		) const;
	int  ParaDividing(
		const array2d &,
		array2d &
		) const;
	int  ParaProducting(
		const array2d &,
		array2d &,
		const array2d &ref
		) const;
	int  ParaDividing( 
		const array2d &,
		array2d &,
		const array2d &ref
		) const;
};
