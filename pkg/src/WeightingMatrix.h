#pragma once
#include "DataMatrix.h"

class WeightingMatrix :
	public DataMatrix
{
public:

	WeightingMatrix(void)
	{
	}

	WeightingMatrix(int nRows, int nColumns) 
		: DataMatrix(nRows, nColumns)
	{}

	~WeightingMatrix(void)
	{
	}

	int	TakeLog();
	int TakeExponential();
	void Deviation(array2d &);
	int ComputeData(
		const array2d &tobeScaled, 
		const array2d &scaler
		);

};
