#include "WeightingMatrix.h"

using namespace std;

int WeightingMatrix::TakeLog()
{
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] = log(data[i][j]);
		}
	}

	return 1;
}

int WeightingMatrix::TakeExponential()
{
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] = exp(data[i][j]);
		}
	}

	return 1;
}


void WeightingMatrix::Deviation(array2d &vec2)
{
	vec2.resize( nRows, array1d(nColumns) );
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			vec2[i][j] = fabs(data[i][j] - 1.0);
		}
	}
}

int WeightingMatrix::ComputeData(
							const array2d &tobeScaled, 
							const array2d &scaler
							)
{
	for (unsigned i=0; i<(unsigned)nRows; i++)
	{
		for (unsigned j=0; j<(unsigned)nColumns; j++)
		{
			// try
			{
				data[i][j] = tobeScaled[i][j] / scaler[i][j];
			}
			/*
			catch( std::runtime_error &e )
			{
				std::cerr << "DataMatrix::ComputeData(const array2d &, const array2d &, array2d &): " << e.what() << endl;
				return -1;
			}
			*/
		}
	}

	return 1;
}

