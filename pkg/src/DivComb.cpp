#include "DivComb.h"

using namespace std;


void DivComb::UpdatePatternMatrix()
{
	// CalculateWeights();
	// weightMatrix = new WeightingMatrix(myData.nRows, myData.nColumns);
	myData.Dividing( mock->data, weightMatrix->data );

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
	double div2 = myData.Div2(mock->data);

	if(div1 > div2) 
	{
		pattern->Scaling( w2.data );
	}
	else
	{
		pattern->Scaling( w1.data ); 
	}

	// delete weightMatrix;
}

void DivComb::UpdateAmplitudeMatrix()
{
	// CalculateWeights();
	// weightMatrix = new WeightingMatrix(myData.nRows, myData.nColumns);
	myData.Dividing( mock->data, weightMatrix->data );
	pattern->GetModes();

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
	double div2 = myData.Div2(mock->data);

	if(div1 > div2) 
	{
		amplitude->Scaling( p2.data );
	}
	else
	{
		amplitude->Scaling( p1.data );
	}

	// delete weightMatrix;
}

void DivComb::UpdateAmplitudeMatrix(WeightingMatrix &p)
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

void DivComb::UpdatePatternMatrix(WeightingMatrix &w)
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


/*
void DivComb::UpdatePatternMatrix()
{
	// CalculateWeights();
	amplitude->GetModes();

	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double weightedComb1 = 0.0;
			double weightedComb2 = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				weightedComb1 += (myData.data[i][j] / mock->data[k][j]) * amplitude->data[k][i];
				weightedComb2 += (mock->data[k][j] / myData.data[i][j]) * amplitude->data[k][i];

			}
			if(amplitude->modes[i] != 0.0)
			{
				pattern->data[i][j] *= (weightedComb1 / amplitude->modes[i] + amplitude->modes[i] / weightedComb2) / 2.0;
			}
		}
	}
}

void DivComb::UpdateAmplitudeMatrix()
{
	// CalculateWeights();
	pattern->GetModes();

	for (int j=0; j<rank; j++) 
	{
		for (int i=0; i<myControl.nObservations; i++)
		{
			double weightedComb1 = 0.0;
			double weightedComb2 = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				weightedComb1 += (myData.data[i][j] / mock->data[i][j]) * pattern->data[j][k];
				weightedComb2 += (mock->data[i][j] / myData.data[i][j]) * pattern->data[j][k];
			}
			
			if(pattern->modes[j] != 0.0)
			{
				amplitude->data[i][j] *= (weightedComb1 / pattern->modes[j] + pattern->modes[j] / weightedComb2) / 2.0;
			}
		}
	}

}
*/

/*
void DivComb::UpdateAmplitudeMatrix( array2d &updateFactor )
{
	updateFactor.resize( myControl.nObservations, array1d(rank) );

	pattern->GetModes();
	// need to save mode for pattern and amplitude
	for (int j=0; j<rank; j++) 
	{
		for (int i=0; i<myControl.nObservations; i++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nParameters; k++)
			{
				weightedComb += weightMatrix->data[i][k] * pattern->data[j][k];
			}
			updateFactor[i][j] = pattern->modes[j] / weightedComb;
		}
	}
}


void DivComb::UpdatePatternMatrix( array2d &updateFactor )
{
	updateFactor.resize( rank, array1d(myControl.nParameters) );

	amplitude->GetModes();
	for (int i=0; i<rank; i++) 
	{
		for (int j=0; j<myControl.nParameters; j++)
		{
			double weightedComb = 0.0;
			for (int k=0; k<myControl.nObservations; k++)
			{
				weightedComb += weightMatrix->data[k][j] * amplitude->data[k][i];
			}
			updateFactor[i][j] = amplitude->modes[i] / weightedComb;
		}
	}
}

*/
