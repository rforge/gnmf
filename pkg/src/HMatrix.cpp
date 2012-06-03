#include "HMatrix.h"

using namespace std;

void HMatrix::Normalizing( const array1d & modes )
{
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] /= modes[i];
		}
	}
}

void HMatrix::GetModes()
{
	modes.resize(nRows);
	for (int i=0; i<nRows; i++)
	{
		modes[i] = 0.0;
		for (int j=0; j<nColumns; j++)
		{
			modes[i] += data[i][j];
		}
	}
}

void HMatrix::GetModes(const double bdn, const vector<double> &v)
{
	modes.resize(nRows);
	for (int i=0; i<nRows; i++)
	{
		modes[i] = 0.0;
		for (int j=0; j<nColumns; j++)
		{
			modes[i] += data[i][j] * (bdn - v[j]) / bdn;
		}
	}
}


int HMatrix::ParaProducting( 
		const array2d &vec2,
		array2d &product
	) const
{
	if(nColumns != (int)vec2[0].size())
	{
		cout << "MatrixParaProduct has different columns for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( vec2.size(), array1d(nRows) );

	for (int i=0; i<(int)vec2.size(); i++)
	{
		for (int j=0; j<nRows; j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nColumns; k++)
			{
				product[i][j] += data[j][k] * vec2[i][k];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}

	}

	return 1;
}

int HMatrix::ParaDividing( 
		const array2d &vec2,
		array2d &product
	) const
{
	if(nColumns != (int)vec2[0].size())
	{
		cout << "MatrixParaProduct has different columns for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( vec2.size(), array1d(nRows) );

	for (int i=0; i<(int)vec2.size(); i++)
	{
		for (int j=0; j<nRows; j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nColumns; k++)
			{
				product[i][j] += data[j][k] / vec2[i][k];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}
	}

	return 1;
}



int HMatrix::ParaProducting( 
		const array2d &vec2,
		array2d &product,
		const array2d &ref
	) const
{
	if(nColumns != (int)vec2[0].size())
	{
		cout << "MatrixParaProduct has different columns for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( vec2.size(), array1d(nRows) );

	for (int i=0; i<(int)vec2.size(); i++)
	{
		for (int j=0; j<nRows; j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nColumns; k++)
			{
				if(fabs(ref[i][k]) < OUTLINERTHRESHOLD) product[i][j] += data[j][k] * vec2[i][k];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}
	}

	return 1;
}

int HMatrix::ParaDividing( 
		const array2d &vec2,
		array2d &product,
		const array2d &ref
	) const
{
	if(nColumns != (int)vec2[0].size())
	{
		cout << "MatrixParaProduct has different columns for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( vec2.size(), array1d(nRows) );

	for (int i=0; i<(int)vec2.size(); i++)
	{
		for (int j=0; j<nRows; j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nColumns; k++)
			{
				if(fabs(ref[i][k]) < OUTLINERTHRESHOLD) product[i][j] += data[j][k] / vec2[i][k];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}
				
	}

	return 1;
}
