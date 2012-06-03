#include "VMatrix.h"

using namespace std;

void VMatrix::Idealize( double bottom )
{
	if(bottom > 0.0)
	{
		minValue = numeric_limits<double>::max();
		// for (unsigned i=0; i<nRows; i++)
		for (unsigned i=0; i<data.size(); i++)
		{
			// for (unsigned j=0; j<nColumns; j++)
			if((int)data[i].size() < nColumns)
			{
				cout << "row " << i << " has less columns " << data[i].size() << endl;
			}
			for (unsigned j=0; j<data[i].size(); j++)
			{
				if(data[i][j] > 0.0 && data[i][j] < minValue)
				{
					minValue = data[i][j];
				}
			}
		}

		minValue /= 10.0;
	}
	else
	{
		minValue = bottom;
	}
	
	if(HD)
	{
		for (unsigned i=0; i<hdIndex.size(); i++)
		{
			for (unsigned j=0; j<hdIndex[i].size(); j++)
			{
				if(hdValue[i][j] < minValue) hdValue[i][j] = minValue;
			}
		}
	}
	else
	{
		for (unsigned i=0; i<data.size(); i++)
		{
			for (unsigned j=0; j<data[i].size(); j++)
			{
				if(data[i][j] < minValue) data[i][j] = minValue;
			}
		}
	}
}

void VMatrix::Mock2Weight( const array2d &tobeScaled )
{
	for (unsigned i=0; i<data.size(); i++)
	{
		for (unsigned j=0; j<data[i].size(); j++)
		{
			if(tobeScaled[i][j] == 0.0)
			{
				data[i][j] = 0.0;
			}
			else
			{
				data[i][j] = tobeScaled[i][j] / data[i][j];
			}
		}
	}
}

int VMatrix::Dividing(
		const array2d &scaler,
		array2d &result
	) const
{
	result.resize( nRows, array1d(nColumns) );

	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			result[i][j] = data[i][j] / scaler[i][j];
		}
	}

	return 1;
}


double VMatrix::EDist(array2d &vec) const
{
	double ed = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			double diff = vec[i][j] - data[i][j];
			ed += diff * diff;
		}
	}

	return sqrt(ed);
}

double VMatrix::Div2(const array2d &vec2) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			if(data[i][j] == 0.0)
			{
				klValue += vec2[i][j];
			}
			else
			{
				klValue += vec2[i][j] * log( vec2[i][j] / data[i][j] ) - vec2[i][j] + data[i][j];
			}
		}
	}
	return klValue;
}

double VMatrix::KL_Div(const array2d &vec2) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			if(data[i][j] == 0.0)
			{
				klValue += vec2[i][j];
			}
			else
			{
				klValue += data[i][j] * log( data[i][j] / vec2[i][j] ) - data[i][j] + vec2[i][j];
			}
		}
	}
	return klValue;
}

double VMatrix::BDist(const array2d &vec2, const double bdn) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			if(data[i][j] == 0.0)
			{
				klValue += (bdn - data[i][j]) * log( (bdn - data[i][j]) / (bdn - vec2[i][j]) );
			}
			else
			{
				klValue += data[i][j] * log( data[i][j] / vec2[i][j] ) +  (bdn - data[i][j]) * log( (bdn - data[i][j]) / (bdn - vec2[i][j]) );
			}
		}
	}
	return klValue;
}

double VMatrix::NBDist(const array2d &vec2, const double bdn) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			if(data[i][j] == 0.0)
			{
				klValue -=  (bdn + data[i][j]) * log( (bdn + data[i][j]) / (bdn + vec2[i][j]) );
			}
			else
			{
				klValue += data[i][j] * log( data[i][j] / vec2[i][j] ) -  (bdn + data[i][j]) * log( (bdn + data[i][j]) / (bdn + vec2[i][j]) );
			}
		}
	}
	return klValue;
}

double VMatrix::GammaKL_Div(const array2d &vec2) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			if(data[i][j] == 0.0)
			{
				klValue += - 1.0; // ???
			}
			else
			{
				klValue += log( data[i][j] / vec2[i][j] ) + vec2[i][j] / data[i][j] - 1.0;
			}
		}
	}
	return klValue;
}

double VMatrix::Pareto_Div(const array2d &vec2, const double alpha) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			klValue += log( data[i][j] / (alpha * vec2[i][j]) ) + vec2[i][j] * log( data[i][j] / alpha );
		}
	}
	return klValue;
}

double VMatrix::InverseLink_GammaKL_Div(const array2d &vec2) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			double tmp =  data[i][j] * vec2[i][j];
			klValue +=  tmp - log( tmp ) - 1.0;
		}
	}
	return klValue;
}

double VMatrix::Gamma_Div(const array2d &vec2) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			double diff = data[i][j] / vec2[i][j];
			klValue += diff - log(diff) - 1.0;
		}
	}
	return klValue;
}


double VMatrix::GammaJD_Div(const array2d &vec2) const
{
	double klValue = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			double diff = vec2[i][j] - data[i][j];
			klValue += diff * diff / (data[i][j] * vec2[i][j]);
		}
	}
	return klValue;
}

double VMatrix::Renyis(const array2d &vec2, const double alpha) const
{
	if(alpha == 1.0 || alpha == 0.0) return KL_Div(vec2);

	double klValue = 0.0;
	double alpha_ = 1.0 - alpha;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			klValue += pow(data[i][j], alpha) * pow(vec2[i][j], alpha_) - alpha * data[i][j] - alpha_ * vec2[i][j];
			/*
			if(vec2[i][j] == 0.0) 
			{
				klValue -= alpha * data[i][j];
			}
			else if(data[i][j] == 0.0)
			{
				klValue -= alpha_ * vec2[i][j];
			}
			else
			{
				// double tmp = alpha * log(data[i][j]) + alpha_ * log(vec2[i][j]);
				// klValue += exp(tmp) - alpha * data[i][j] - alpha_ * vec2[i][j];		
				klValue += pow(data[i][j] / vec2[i][j], alpha) * vec2[i][j] - alpha * data[i][j] - alpha_ * vec2[i][j];
			}
			*/
		}
	}
	return klValue / (alpha * (alpha  - 1.0));
}

double VMatrix::ODPDivergence( const array2d &vec2, const array1d &phi) const
{
	double d = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			if(data[i][j] == 0.0)
			{
				d += vec2[i][j] / phi[i];
			}
			else
			{
				d += (data[i][j] * log(data[i][j] / vec2[i][j]) - data[i][j] + vec2[i][j]) / phi[i];
			}
		}
	}

	return d;
}

double VMatrix::ChiSquare(
			const array2d &mock,
			const array2d &uncert
			)
{
	double chiSquare = 0.0;
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			double diffEle = (data[i][j] - mock[i][j]) / uncert[i][j];
			diffEle = diffEle * diffEle;
			chiSquare += diffEle;
		}
	}

	return sqrt(chiSquare);
}
