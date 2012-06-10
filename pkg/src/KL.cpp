#include "KL.h"

using namespace std;
void KL::UpdatePatternMatrix()
{
	CalculateWeights();
	amplitude->GetModes();

	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				// weightedComb += weightMatrix->data[k][j] * amplitude->data[k][i];
				weightedComb += mock->data[k][j] * amplitude->data[k][i];
			}
			if(weightedComb < 0.0)
			{
//				cout << "weightedComb " << weightedComb << endl;
                ;
			}
			if(amplitude->modes[i] != 0.0)
			{
				// updateFactorP->data[i][j] = weightedComb / amplitude->modes[i];
				pattern->data[i][j] *= weightedComb / amplitude->modes[i];
			}
		}
	}
}

void KL::UpdateAmplitudeMatrix()
{
	CalculateWeights();
	pattern->GetModes();

	for (int j=0; j<rank; j++) 
	{
		for (int i=0; i<myControl.nObservations; i++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				// weightedComb += weightMatrix->data[i][k] * pattern->data[j][k];
				weightedComb += mock->data[i][k] * pattern->data[j][k];
			}
			if(weightedComb < 0.0)
			{
//				cout << "weightedComb " << weightedComb << endl;
                ;
			}
			if(pattern->modes[j] != 0.0)
			{
				amplitude->data[i][j] *= weightedComb / pattern->modes[j];
			}
		}
	}
}

