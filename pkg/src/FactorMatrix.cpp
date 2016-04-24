#include "FactorMatrix.h"

using namespace std;

void FactorMatrix::Initializing()
{
    // JMM (12/7/2014): Replace call to system random number generator
    // with call to internal R function
	// Reference: http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Random-numbers
	GetRNGstate();

	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++) 
		{
            // JMM (12/7/2014): Replace call to system random number generator
            // with call to internal R function
			data[i][j] = unif_rand();
		}
	}

    PutRNGstate();              
}

int FactorMatrix::Multipling( 
		const array2d &vec2,
		array2d &product
	) const
{
	if(nColumns != (int)vec2.size())
	{
//		cout << "MatrixProduct has different common dimensions for vec1 and vec2!" << endl;
		return -1;
	}

	product.resize( nRows, array1d(vec2[0].size()) );

	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<(int)vec2[0].size(); j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<nColumns; k++)
			{
				product[i][j] += data[i][k] * vec2[k][j];
			}

			// no zero element is allowed
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}
	}

	return 1;
}

int FactorMatrix::Multipling( 
		const array2d &vec1,
		const array2d &vec2,
		array2d &product
	) const
{
	if(nColumns != (int)vec2.size() || vec2[0].size() != vec1.size())
	{
//		cout << "MatrixProduct has different common dimensions for vec1 and vec2!" << endl;
		return -1;
	}

	array2d tmp(nRows, array1d(vec1.size()));
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<(int)vec1.size(); j++)
		{
			tmp[i][j] = 0.0;
			for (int k=0; k<nColumns; k++)
			{
				tmp[i][j] += data[i][k] * vec2[k][j];
			}
		}
	}

	product.resize( nRows, array1d(vec1[0].size()) );
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<(int)vec1[0].size(); j++)
		{
			product[i][j] = 0.0;
			for (int k=0; k<(int)vec1.size(); k++)
			{
				product[i][j] += tmp[i][k] * vec1[k][j];
			}
			if(product[i][j] == 0.0) product[i][j] = EPS;
		}
	}

	return 1;
}

double FactorMatrix::GetSparseness()
{
	double absValue = 0.0;
	double sqrtValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			absValue += data[i][j]; // because all elements are non-negative
			sqrtValue += data[i][j] * data[i][j];
		}
	}

	sqrtValue = sqrt( sqrtValue );
	if(sqrtValue == 0.0) sqrtValue = EPS; // in case all element is close to zero, this still give a valid value
	double sqrtN = sqrt( static_cast<double>(nRows * nColumns) );
	
	sparseness = (sqrtN - absValue / sqrtValue) / (sqrtN - 1);

	return sparseness;
}

void FactorMatrix::GetConnectivity()
{
	connectivity = new CMatrix(nColumns);

	vector<int> index(nColumns);
	for (int i=0; i<nColumns; i++)
	{
		double max = data[0][i];
		index[i] = 0;
		for (int j=0; j<nRows; j++)
		{
			if(data[j][i] > max)
			{
				max = data[j][i];
				index[i] = j;
			}
		}
	}

	// build connectivity matrix and calculate consensus matrix
	for (int i=0; i<nColumns; i++)
	{
		connectivity->data[i][i] = 1;
		for (int j=i+1; j<nColumns; j++) 
		{
			if(index[i] == index[j]) 
			{
				connectivity->data[i][j] = 1;
			}
			else
			{
				connectivity->data[i][j] = 0;
			}
			connectivity->data[j][i] = connectivity->data[i][j];
		}
	}
}

void FactorMatrix::GetConnectivityReal( string scheme )
{
	if(scheme.find("VDP", 0) != scheme.npos)
	{
		GetConnectivityVDP();	
	}
	else
	{
		GetConnectivityPearsons();
	}
}


void FactorMatrix::GetConnectivityVDP()
{
	realConnet = new CMatrix(nColumns);

	// build connectivity matrix by applying to actual data instead of binary values
	vector<double> cModes(nColumns);
	for (int i=0; i<nColumns; i++)
	{
		double mode = 0.0;
		for (int k=0; k<nRows; k++)
		{
			mode += data[k][i] * data[k][i];
		}
		// cModes.push_back( sqrt(mode) );
		cModes[i] = sqrt(mode);
		// cout << "nMode i " << i << " " << cModes[i] <<  " from mode " << mode << endl;
	}
	for (int i=0; i<nColumns; i++)
	{
		realConnet->data[i][i] = 1;
		for (int j=i+1; j<nColumns; j++) 
		{
			realConnet->data[i][j] = 0.0;
			for (int k=0; k<nRows; k++)
			{
				realConnet->data[i][j] += data[k][i] * data[k][j];
			}
			realConnet->data[i][j] /= (cModes[i] * cModes[j]);
			realConnet->data[j][i] = realConnet->data[i][j];
		}
	}
}

void FactorMatrix::GetConnectivityPearsons()
{
	realConnet = new CMatrix(nColumns);

	// build connectivity matrix by applying to actual data instead of binary values
	vector<double> sum(nColumns);
	vector<double> var(nColumns);
	for (int i=0; i<nColumns; i++)
	{
		sum[i] = 0.0;
		var[i] = 0.0;
		for (int k=0; k<nRows; k++)
		{
			sum[i] += data[k][i];
			var[i] += data[k][i] * data[k][i];
		}
		var[i] = sqrt(nRows * var[i] - sum[i] * sum[i]);
	}

	for (int i=0; i<nColumns; i++)
	{
		realConnet->data[i][i] = 1;

		// check var
		if(var[i] == 0.0) continue;

		for (int j=i+1; j<nColumns; j++) 
		{
			realConnet->data[i][j] = 0.0;
			if(var[j] == 0.0) continue;

			for (int k=0; k<nRows; k++)
			{
				realConnet->data[i][j] += data[k][i] * data[k][j];
			}

			realConnet->data[i][j] = (nRows * realConnet->data[i][j] - sum[i] * sum[j] ) / (var[i] * var[j]);
			realConnet->data[j][i] = realConnet->data[i][j];
		}
	}
}
