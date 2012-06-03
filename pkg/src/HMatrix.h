#pragma once
#include "FactorMatrix.h"

class HMatrix :
	public FactorMatrix
{
	// friend std::istream &operator >> (std::istream &in, HMatrix &myMatrix);
	// friend std::ostream &operator << (std::ostream &out, const HMatrix &myMatrix);
public:
	HMatrix(void) {}

	HMatrix(int nRows, int nColumns, long seed) 
		: FactorMatrix(nRows, nColumns, seed)
	{}

	HMatrix(const int nRows, const int nColumns) 
		: FactorMatrix(nRows, nColumns)
	{}


	HMatrix( const DataMatrix & _matrix )
		: FactorMatrix( _matrix )
	{}


	~HMatrix(void)
	{
	}

	// void Initializing( long baseSeed );
	void Normalizing( const array1d & modes );
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
private:
};
