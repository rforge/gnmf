#include "BD.h"

// JMM (7/11/2016): R header included in BD.h.

using namespace std;

void BD::UpdatePatternMatrix()
{
	CalculateWeights();

	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double upper = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				upper += mock->data[k][j] * amplitude->data[k][i];
			}

			double lower = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				lower += (1 - myData.data[k][j] / bdn) * amplitude->data[k][i];
			}
			if(upper < 0.0)
			{
//				cout << "weightedComb " << upper << endl;
                ;
			}
			if(lower != 0.0)
			{
				pattern->data[i][j] *= upper / lower;
			}
		}
	}
}

void BD::UpdateAmplitudeMatrix()
{
	CalculateWeights();

	for (int j=0; j<rank; j++) 
	{
		for (int i=0; i<myControl.nObservations; i++)
		{
			double upper = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				upper += mock->data[i][k] * pattern->data[j][k];
			}

			double lower = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				lower += (1 - myData.data[i][k] / bdn) * pattern->data[j][k];
			}

			if(upper < 0.0)
			{
//				cout << "weightedComb " << upper << endl;
                ;
			}
			if(lower != 0.0)
			{
				amplitude->data[i][j] *= upper / lower;
			}
		}
	}
}

int BD::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0)
	{
		converged = false;
		return 0;
	}
	else 
	{
		errors = myData.BDist( mock->data, bdn );
		consecutiveError = fabs(errors - myData.BDist( prevMock, bdn ));
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
