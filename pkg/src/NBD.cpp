#include "NBD.h"

using namespace std;

void NBD::UpdatePatternMatrix()
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
				lower += (1 + myData.data[k][j] / bdn) * amplitude->data[k][i];
			}
			if(upper < 0.0)
			{
				cout << "weightedComb " << upper << endl;
			}
			if(lower != 0.0)
			{
				pattern->data[i][j] *= upper / lower;
			}
		}
	}
}

void NBD::UpdateAmplitudeMatrix()
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
				lower += (1 + myData.data[i][k] / bdn) * pattern->data[j][k];
			}

			if(upper < 0.0)
			{
				cout << "weightedComb " << upper << endl;
			}
			if(lower != 0.0)
			{
				amplitude->data[i][j] *= upper / lower;
			}
		}
	}
}


int NBD::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0 )
	{
		converged = false;
		return 0;
	}
	else 
	{
		errors = myData.NBDist( mock->data, bdn );
		consecutiveError = fabs(errors - myData.NBDist( prevMock, bdn ));
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

