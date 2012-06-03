#include "Pareto.h"

using namespace std;

void Pareto::UpdatePatternMatrix()
{
	for (int i=0; i< amplitude->nColumns; i++)
	{
		for (int j=0; j<myData.nColumns; j++)
		{
			double upper = 0.0;
			double lower = 0.0;
			for (int k=0; k<amplitude->nRows; k++)
			{
				upper += amplitude->data[k][i] / mock->data[k][j];
				lower += amplitude->data[k][i] * log(myData.data[k][j] / alpha);
			}
			if(upper != 0.0 && lower != 0.0)
			{
				pattern->data[i][j] *= upper / lower;
			}
		}
	}
}

void Pareto::UpdateAmplitudeMatrix()
{
	for (int i=0; i<myData.nRows; i++)
	{
		for (int j=0; j<pattern->nRows; j++)
		{
			double upper = 0.0;
			double lower = 0.0;
			for (int k=0; k<pattern->nColumns; k++)
			{
				upper += pattern->data[j][k] / mock->data[i][k];
				lower += pattern->data[j][k] * log(myData.data[i][k] / alpha);
			}
			if(upper != 0.0 && lower != 0.0)
			{
				amplitude->data[i][j] *= upper / lower;
			}
		}
	}
}

int Pareto::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0 )
	{
		converged = false;
		return 0;
	}
	else 
	{
		errors = myData.Pareto_Div( mock->data, alpha );
		consecutiveError = fabs(errors - myData.Pareto_Div( prevMock, alpha ));
		totalSteps = round;
		if( consecutiveError < CONVERGETHRESHOLD)
		{
			converged = true;
			return 1;
		}
		else
		{
			converged = false;
			return 0;
		}
	}
}
