#include "ODP.h"

using namespace std;


void ODP::SaveIntermediaFiles()
{
	string phiFileName = wMatrixName + "_phi";
	ofstream aOUT(phiFileName.c_str(), ios::out);
	for (int i=0; i<(int)phi.size(); i++)
	{
		aOUT << i << " " << phi[i] << endl;
	}
	aOUT.close();

}

void ODP::UpdatePatternMatrix()
{
	phi.resize( myControl.nObservations );
	for (int i=0; i<myControl.nObservations; i++)
	{
		phi[i] = 0.0;
		for (int j=0; j<myControl.nParameters; j++)
		{
			double tmp =  myData.data[i][j] - mock->data[i][j]; 
			tmp = tmp * (tmp / mock->data[i][j]);
			phi[i] += tmp;
		}
		phi[i] /= (myControl.nParameters - rank);
	}

	CalculateWeights();

	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double upper = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				upper += mock->data[k][j] * amplitude->data[k][i] / phi[k] ;
			}

			double lower = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				lower += amplitude->data[k][i] / phi[k];
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

void ODP::UpdateAmplitudeMatrix()
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

int ODP::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0 )
	{
		converged = false;
		return 0;
	}
	else 
	{
		errors = myData.ODPDivergence( mock->data, phi );
		consecutiveError = fabs(errors - myData.ODPDivergence( prevMock, phi ));
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
