#pragma once
#include "DataMatrix.h"
#include <limits>

class VMatrix :
	public DataMatrix
{
	// friend std::istream &operator >> (std::istream &in, VMatrix &myMatrix);
	// friend std::ostream &operator << (std::ostream &out, const VMatrix &myMatrix);

public:
	VMatrix(void) {}
	VMatrix(int _nRows, int _nColumns) 
		: DataMatrix(_nRows, _nColumns) 
	{}

	~VMatrix(void) {}

	double EDist(array2d &) const;
	double KL_Div(const array2d &) const;
	double BDist(const array2d &, const double) const;
	double NBDist(const array2d &, const double) const;
	double Renyis(const array2d &, const double alpha) const;
	double InverseLink_GammaKL_Div(const array2d &) const;
	double Pareto_Div(const array2d &, const double alpha) const;
	double GammaKL_Div(const array2d &) const;
	double GammaJD_Div(const array2d &) const;
	double Gamma_Div(const array2d &) const;
	double Div2(const array2d &) const;
	double ODPDivergence(const array2d &, const array1d &) const;
	double ChiSquare(
		const array2d &, 
		const array2d &
		);
	void Idealize( double bottom = 0.0 );
	void Mock2Weight(const array2d &);
	int  Dividing(
		const array2d &,
		array2d &
		) const;

private:
	double minValue;

};
