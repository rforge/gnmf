#include "Div2.h"

using namespace std;

/*
void Div2::UpdatePatternMatrix()
{
	CalculateWeights();
	amplitude->GetModes();

	HMatrix h1(*pattern);
	HMatrix h2(*pattern);

	WeightingMatrix w1(rank, myControl.nParameters);
	UpdatePatternMatrix( w1 );
	h1.Scaling( w1.data );
	GenerateMock(*amplitude, h1);
	double div1 = myData.KL_Div( mock->data );

	// div2 P matrix
	weightMatrix->TakeLog();

	WeightingMatrix w2(rank, myControl.nParameters);
	UpdatePatternMatrix(w2);
	w2.TakeExponential();
	h2.Scaling(w2.data);
	GenerateMock(*amplitude, h2);
	double div2 = myData.KL_Div(mock->data);

	if(div1 > div2) 
	{
		updateFactorP->data = w2.data;
	}
	else
	{
		updateFactorP->data = w1.data;
	}
}

void Div2::UpdatePatternMatrix(WeightingMatrix &w)
{
	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				weightedComb += weightMatrix->data[k][j] * amplitude->data[k][i];
			}
			if(amplitude->modes[i] == 0.0)
			{
				w.data[i][j] = 1.0;
			}
			else
			{
				w.data[i][j] = weightedComb / amplitude->modes[i];
			}
		}
	}
}

void Div2::UpdateAmplitudeMatrix()
{
	CalculateWeights();
	amplitude->GetModes();

	WMatrix w1(*amplitude);
	WMatrix w2(*amplitude);

	WeightingMatrix p1(myControl.nObservations, rank);
	UpdateAmplitudeMatrix( p1 );
	w1.Scaling(p1.data);
	GenerateMock(w1, *pattern);
	double div1 = myData.KL_Div( mock->data );

	// div2 P matrix
	weightMatrix->TakeLog();

	WeightingMatrix p2(myControl.nObservations, rank);
	UpdateAmplitudeMatrix(p2);
	p2.TakeExponential();
	w2.Scaling(p2.data);
	GenerateMock(w2, *pattern);
	double div2 = myData.KL_Div(mock->data);

	if(div1 > div2) 
	{
		updateFactorA->data = p2.data;
	}
	else
	{
		updateFactorA->data = p1.data;
	}
}

void Div2::UpdateAmplitudeMatrix(WeightingMatrix &p)
{
	for (int j=0; j<rank; j++) 
	{
		for (int i=0; i<myControl.nObservations; i++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				weightedComb += weightMatrix->data[i][k] * pattern->data[j][k];
			}
			
			if(pattern->modes[j] == 0.0)
			{
				p.data[i][j] = 1.0;
			}
			else
			{
				p.data[i][j] = weightedComb / pattern->modes[j];
			}
		}
	}
}
*/

int Div2::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0)
	{
		converged = false;
		return 0;
	}
	else
	{
		errors = myData.Div2( mock->data );
		consecutiveError = fabs(errors - myData.Div2( prevMock ));
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

void Div2::UpdatePatternMatrix()
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
				weightedComb += log(mock->data[k][j]) * amplitude->data[k][i];
			}
			if(amplitude->modes[i] != 0.0)
			{
				pattern->data[i][j] = exp( weightedComb / amplitude->modes[i] );
			}
		}
	}
}

void Div2::UpdateAmplitudeMatrix()
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
				weightedComb += log( mock->data[i][k] ) * pattern->data[j][k];
			}
			
			if(pattern->modes[j] != 0.0)
			{
				amplitude->data[i][j] = exp( weightedComb / pattern->modes[j] );
			}
		}
	}
}
