#include "Renyi.h"

using namespace std;

void Renyi::UpdatePatternMatrix()
{
	CalculateWeights();
	mock->TakePower( alpha );
	amplitude->GetModes();

	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				weightedComb += mock->data[k][j] * amplitude->data[k][i];
			}
			if(amplitude->modes[i] != 0.0)
			{
				pattern->data[i][j] *= Utility::PowerUp(weightedComb/amplitude->modes[i], alphaR);
			}
		}
	}
}

void Renyi::UpdateAmplitudeMatrix()
{
	CalculateWeights();
	mock->TakePower( alpha );
	pattern->GetModes();

	for (int j=0; j<rank; j++) 
	{
		for (int i=0; i<myControl.nObservations; i++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				weightedComb += mock->data[i][k] * pattern->data[j][k];
			}
			
			if(pattern->modes[j] != 0.0)
			{
				amplitude->data[i][j] *= Utility::PowerUp(weightedComb/pattern->modes[j], alphaR);
			}
		}
	}
}

int Renyi::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0)
	{
		converged = false;
		return 0;
	}
	else 
	{
		errors = myData.Renyis( mock->data, alpha );
//		double klError = myData.KL_Div( mock->data );
		consecutiveError = fabs(errors - myData.Renyis( prevMock, alpha  ));
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


