#include "WMatrix.h"

using namespace std;

int WMatrix::ParaProducting( 
		const array2d &vec2,
		array2d &product
	) const
{
	if(nRows != (int)vec2.size())
	{
//		cout << "MatrixParaProduct has different rows for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( nColumns, array1d(vec2[0].size()) );

	for (int i=0; i<nColumns; i++)
	{
		for (int j=0; j<(int)vec2[0].size(); j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nRows; k++)
			{
				product[i][j] += data[k][i] * vec2[k][j];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}

	}

	return 1;
}

int WMatrix::ParaDividing( 
		const array2d &vec2,
		array2d &product
	) const
{
	if(nRows != (int)vec2.size())
	{
//		cout << "MatrixParaProduct has different rows for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( nColumns, array1d(vec2[0].size()) );

	for (int i=0; i<nColumns; i++)
	{
		for (int j=0; j<(int)vec2[0].size(); j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nRows; k++)
			{
				product[i][j] += data[k][i] / vec2[k][j];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}

	}

	return 1;
}


int WMatrix::ParaProducting(
		const array2d &vec2,
		array2d &product,
		const array2d &ref
	) const
{
	if(nRows != (int)vec2.size())
	{
//		cout << "MatrixParaProduct has different rows for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( nColumns, array1d(vec2.size()) );

	for (int i=0; i<nColumns; i++)
	{
		for (int j=0; j<(int)vec2[0].size(); j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nRows; k++)
			{
				if(fabs(ref[k][j]) < OUTLINERTHRESHOLD) product[i][j] += data[k][i] * vec2[k][j];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}

	}

	return 1;
}

int WMatrix::ParaDividing(
		const array2d &vec2,
		array2d &product,
		const array2d &ref
	) const
{
	if(nRows != (int)vec2.size())
	{
//		cout << "MatrixParaProduct has different rows for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( nColumns, array1d(vec2[0].size()) );

	for (int i=0; i<nColumns; i++)
	{
		for (int j=0; j<(int)vec2[0].size(); j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nRows; k++)
			{
				if(fabs(ref[k][j]) < OUTLINERTHRESHOLD)  product[i][j] += data[k][i] / vec2[k][j];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}

	}

	return 1;
}


void WMatrix::Normalizing( const array1d & modes )
{
	for (int i=0; i<nColumns; i++)
	{
		for (int j=0; j<nRows; j++)
		{
			data[j][i] *= modes[i];
		}
	}
}

void WMatrix::GetModes()
{
	modes.resize(nColumns);
	for (int j=0; j<nColumns; j++)
	{
		modes[j] = 0.0;
		for (int i=0; i<nRows; i++)
		{
			modes[j] += data[i][j];
		}
	}
}

void WMatrix::GetModes(const double bdn, const vector<double> &v)
{
	modes.resize(nColumns);
	for (int j=0; j<nColumns; j++)
	{
		modes[j] = 0.0;
		for (int i=0; i<nRows; i++)
		{
			modes[j] += data[i][j] * (bdn - v[i]) / bdn;
		}
	}
}
