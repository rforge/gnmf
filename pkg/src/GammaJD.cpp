#include "GammaJD.h"

using namespace std;

void GammaJD::UpdatePatternMatrix()
{
	mock->TakePower(2.0);
	CalculateWeights();

	array2d modes;
	amplitude->ParaDividing( myData.data, modes );
	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				weightedComb += mock->data[k][j] * amplitude->data[k][i];
			}

			pattern->data[i][j] *= sqrt( weightedComb / modes[i][j] );
		}
	}
}

void GammaJD::UpdateAmplitudeMatrix()
{
	mock->TakePower(2.0);
	CalculateWeights();

	array2d modes;
	pattern->ParaDividing(myData.data, modes);

	for (int j=0; j<rank; j++) 
	{
		for (int i=0; i<myControl.nObservations; i++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				weightedComb += mock->data[i][k] * pattern->data[j][k];
			}
			
			amplitude->data[i][j] *= sqrt( weightedComb / modes[i][j] );
		}
	}
}

int GammaJD::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0)
	{
		converged = false;
		return 0;
	}
	else 
	{
		errors = myData.GammaJD_Div( mock->data );
		consecutiveError = fabs(errors - myData.GammaJD_Div( prevMock ));
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
