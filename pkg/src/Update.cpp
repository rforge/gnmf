#include "Update.h"

using namespace std;

// based on NMF divergence measurement, the update rule for pattern matrix is:
//
//             -                                      -
//             |    /         D(k, j)    \            |
//             |sum | ------------------ | * A(k, i)  |
//             | k  \ sum (A(k,l)P(l, j) /            |
//  P(i, j) *= | ------------------------------------ |
//             |                sum A(k, i)           |
//             |                 k                    |
//             -                                      -
//
// so the two "sum(k)" can be calculated in advance, because they only dependent on one of
// the index of the to-be-updated element of pattern matrix.
//
// the same thing is for amplitude matrix.
// 
// in the more generalized situations, the first sum(k) in the upper part is subject to
// a power of alpha.
//
// for updating, we need to calculate weightMatrix (the element of the first sum), the column
// (A) and row (P) mode.

void Update::Run( const ClustDataSet &myRefClust )
{
//Rprintf("    Update::Run() #1: entered function\n"); R_FlushConsole(); R_ProcessEvents();
	// cout << "RunType " << myControl.RunType << endl;
	if(myControl.target.find_first_of("ATTERN") == 1)
	{
//Rprintf("    Update::Run() #1: detected H\n"); R_FlushConsole(); R_ProcessEvents();
		HSimulation( myRefClust );
	}
	else
	{
//Rprintf("    Update::Run() #1: detected W\n"); R_FlushConsole(); R_ProcessEvents();
		WSimulation( myRefClust );
	}
}

void Update::HSimulation( const ClustDataSet &myRefClust )
{
//Rprintf("    Update::HSimulation() #1: entered function\n"); R_FlushConsole(); R_ProcessEvents();

	HEvaluation he( myControl, myRefClust );

	for (int i=0; i<myControl.nChains; i++)
	{
//Rprintf("    Update::HSimulation() #1: about to invoke Run()...\n"); R_FlushConsole(); R_ProcessEvents();
		Run(rank, i, alpha);
		if(myControl.RunType != "simulation")
		{
//Rprintf("    Update::HSimulation() #1: detected ! simulation\n"); R_FlushConsole(); R_ProcessEvents();
			he.targetMatrixes.push_back( pattern );
			amplitude->Transposing();
			he.amplitudes.push_back( amplitude );
		}
	}

	if(myControl.RunType != "simulation") he.Run(rank, alpha);
}

void Update::WSimulation( const ClustDataSet &myRefClust )
{
//Rprintf("    Update::WSimulation() #1: entered function\n"); R_FlushConsole(); R_ProcessEvents();
	WEvaluation we( myControl, myRefClust );

	for (int i=0; i<myControl.nChains; i++)
	{
		Run(rank, i, alpha);
		if(myControl.RunType != "simulation")
		{
			amplitude->Transposing();
			we.targetMatrixes.push_back( amplitude );
			we.patterns.push_back( pattern );
		}
	}

	if(myControl.RunType != "simulation") we.Run(rank, alpha);
}

void Update::Run(const int _rank, const int _chain, const double _alpha = 1.0 )
{
//Rprintf("    ##################################\n");
//Rprintf("    ##################################\n");
//Rprintf("    Update::Run() #2: entered function\n");
//Rprintf("    ##################################\n");
//Rprintf("    ##################################\n");
//R_FlushConsole();
//R_ProcessEvents();
	rank = _rank;
	chain = _chain;
	alpha = _alpha;
	alphaR = 1.0 / alpha;
//Rprintf("    Update::Run() #2: rank   = %d\n",rank);
//Rprintf("    Update::Run() #2: chain  = %d\n",chain);
//Rprintf("    Update::Run() #2: alphaR = %g\n",alphaR);
//Rprintf("    Update::Run() #2: Got this far #50\n");
//R_FlushConsole();
//R_ProcessEvents();

	// cout << "rank " << rank << "  alpha " << alpha << " china " << chain << endl;

	// copy initialized matrixes from myInitial to this object
	GetInitial();
//Rprintf("    Update::Run() #2: Got this far #60\n");
//R_FlushConsole();
//R_ProcessEvents();
	int lastRound = 0;
//Rprintf("    Update::Run() #2: Got this far #100\n");
//R_FlushConsole();
//R_ProcessEvents();

	// do updating simulation
	for (int i=0; i<myControl.nUpdateSteps; i++)
	{
		// cout << "step " << i << " concergeTest " << myControl.ConvergeTest << " testStep " << myControl.convergeTestStep << endl;
		if(myControl.ConvergeTest &&  (i + 1) % myControl.convergeTestStep == 0 ) prevMock = mock->data;
		// cout << "Here" << endl;
		if(i == myControl.nUpdateSteps - 2) prevMock = mock->data;

		if(myControl.theta != 0.0)
		{
			array2d tempW;
			amplitude->Multipling(thetaMatrix->data, tempW);
//Rprintf("    Update::Run() #2: tempW is %d by %d\n",tempW.size(),tempW[0].size());
//R_FlushConsole();
//R_ProcessEvents();
			amplitude->data = tempW;
		}
		// cout << "Update Pattern" << endl;
		UpdatePatternMatrix();
		// pattern->Scaling( updateFactorP->data );

		GenerateMock();

		if(myControl.theta != 0.0)
		{
			array2d tempH;
			thetaMatrix->Multipling(pattern->data, tempH);
			pattern->data = tempH;

		}
		// cout << "Update Amplitude" << endl;
		UpdateAmplitudeMatrix();
		// amplitude->Scaling( updateFactorA->data );

		GenerateMock();
		if(myControl.ConvergeTest && CheckConvergency( i + 1 ) == 1 ) break;
		Normalizing();

		lastRound = i;
	}
//Rprintf("    Update::Run() #2: Got this far #200\n");
//R_FlushConsole();
//R_ProcessEvents();

	if(! converged) CheckConvergency( lastRound );
	UpdateFeatures();

	if(myControl.RunType == "simulation")
	{
		SaveMatrixes();
		delete amplitude;
		delete pattern;
	}
//Rprintf("    Update::Run() #2: Got this far #300\n");
//R_FlushConsole();
//R_ProcessEvents();

	SaveIntermediaFiles();
	Closing();
//Rprintf("    Update::Run() #2: Got this far #400\n");
//R_FlushConsole();
//R_ProcessEvents();
	
}

void Update::SaveMatrixes() const 
{
/*
Rprintf("    Update::SaveMatrixes(): entered function...\n");
Rprintf("    Update::SaveMatrixes() #2: rank   = %d\n",rank);
Rprintf("    Update::SaveMatrixes() #2: chain  = %d\n",chain);
Rprintf("    Update::SaveMatrixes() #2: alphaR = %g\n",alphaR);
Rprintf("    Update::SaveMatrixes(): amplitude is %d by %d\n",amplitude->nRows,amplitude->nColumns);
Rprintf("    Update::SaveMatrixes(): pattern is %d by %d\n",pattern->nRows,pattern->nColumns);
*/
    int sizeAlphaBlock = 0;
    for (int j=myControl.startRank; j<=myControl.endRank; j++ )
        sizeAlphaBlock += j;
    sizeAlphaBlock *= myControl.nChains;
	int numHRows    = sizeAlphaBlock * myControl.nAlphas;
    int alphaOffset = myControl.alphaIndex * sizeAlphaBlock;
    int rankOffset = 0;
    for (int j=myControl.startRank; j<=rank-1; j++ )
        rankOffset += ( j * myControl.nChains );

        // Transfer data from combinedH to R matrix W.
/*
Rprintf("About to transfer data to W matrix...\n");
Rprintf("myData.nRows              = %d\n",myData.nRows);
Rprintf("sizeAlphaBlock            = %d\n",sizeAlphaBlock);
Rprintf("numHRows                  = %d\n",numHRows);
Rprintf("myControl.alphaIndex      = %d\n",myControl.alphaIndex);
Rprintf("myControl.nAlphas         = %d\n",myControl.nAlphas);
Rprintf("myControl.alphaSet.size() = %d\n",myControl.alphaSet.size());
Rprintf("alphaOffset               = %d\n",alphaOffset);
Rprintf("rankOffset                = %d\n",rankOffset);
R_FlushConsole();
R_ProcessEvents();
*/
    int rowIndex, colIndex;
/*
if ( chain == 0 )
{
    for (int j=0; j<amplitude->nRows; j++)
    {
        for (int k=0; k<amplitude->nColumns; k++)
        {
            Rprintf("%g ",amplitude->data[j][k]);
            if ( k > 0 ) Rprintf(" ");
        }
        Rprintf("\n");
    }
        
}
*/

    for (int j=0; j<amplitude->nRows; j++)
        for (int k=0; k<amplitude->nColumns; k++)
        {
//            rowIndex = alphaOffset + rankOffset + ( chain * amplitude->nRows ) + j;
//            colIndex = k;
            rowIndex = j;
            colIndex = alphaOffset + rankOffset + ( chain * amplitude->nColumns ) + k;
//Rprintf("            [ %d, %d ]\n",rowIndex,colIndex);
//            if ( ( j == 0 ) && ( k = 0 ) )
//Rprintf("            [ %d, %d ] gets mapped to [ %d, %d ]\n",j,k,rowIndex,colIndex);
//R_FlushConsole();
//R_ProcessEvents();
//            if ( chain == 0 )
            myControl.W[colIndex*myData.nRows + rowIndex] = amplitude->data[j][k];
        }
//Rprintf("    Update::SaveMatrixes(): wMatrixName.c_str() = %s\n",wMatrixName.data());
//Rprintf("    Update::SaveMatrixes(): amplitude is %d by %d\n",amplitude->nRows,amplitude->nColumns);
//Rprintf("                                %g   %g\n",amplitude->data[0][0],amplitude->data[0][1]);
//Rprintf("                                %g   %g\n",amplitude->data[1][0],amplitude->data[1][1]);
//R_FlushConsole();
//R_ProcessEvents();
//	ofstream aOUT(wMatrixName.c_str(), ios::out);
//	aOUT << *amplitude;
//	aOUT.close();

        // Transfer data from combinedH to R matrix H.
//Rprintf("About to transfer data to H matrix...\n");
//R_FlushConsole();
//R_ProcessEvents();
    for (int j=0; j<pattern->nRows; j++)
        for (int k=0; k<pattern->nColumns; k++)
        {
            rowIndex = alphaOffset + rankOffset + ( chain * pattern->nRows ) + j;
            colIndex = k;
//Rprintf("            [ %d, %d ]\n",rowIndex,colIndex);
//            if ( ( j == 0 ) && ( k == 0 ) ) {
//Rprintf("            [ %d, %d ] gets mapped to [ %d, %d ]\n",j,k,rowIndex,colIndex);
//R_FlushConsole();
//R_ProcessEvents();
//                  }
            myControl.H[colIndex*numHRows + rowIndex] = pattern->data[j][k];            
        }
//Rprintf("    Update::SaveMatrixes(): hMatrixName.c_str() = %s\n",hMatrixName.data());
//Rprintf("    Update::SaveMatrixes(): pattern is %d by %d\n",pattern->nRows,pattern->nColumns);
//Rprintf("                                %g   %g\n",pattern->data[0][0],pattern->data[0][1]);
//Rprintf("                                %g   %g\n",pattern->data[1][0],pattern->data[1][1]);
//R_FlushConsole();
//R_ProcessEvents();
//	ofstream pOUT(hMatrixName.c_str(), ios::out);
//	pOUT << *pattern;
//	pOUT.close();
//Rprintf("    Update::SaveMatrixes(): exiting function...\n");
//R_FlushConsole();
//R_ProcessEvents();
}

void Update::UpdateFeatures()
{
	SSMap properties;
	
	properties["Convergence"] = toString<bool>(converged);
	properties["TotalUpdateSteps"] = toString<int>(totalSteps);
	properties["Error"] = toString<double>(errors);
	properties["ConsecutiveError"] = toString<double>(consecutiveError);
	properties["Wsparseness"] = toString<double>(amplitude->GetSparseness());
	properties["Hsparseness"] = toString<double>(pattern->GetSparseness());

	pattern->SetProperties( properties );
	amplitude->SetProperties( properties );
}

void Update::Closing()
{
	delete mock;
	delete weightMatrix;
	prevMock.clear();
}

void Update::Normalizing()
{
	if(! myControl.Normalize) return;
	array1d mode;
	pattern->GetModes();
	pattern->Normalizing( pattern->modes );
	amplitude->Normalizing( pattern->modes );
}

void Update::GetInitial()
{
//Rprintf("    Update::GetInitial(): entered function...\n");
//R_FlushConsole();
//R_ProcessEvents();
	totalSteps = myControl.nUpdateSteps;
	converged = false;
//Rprintf("    Update::GetInitial(): Got this far #50...\n");
//Rprintf("    nColumns = %d\n",myData.nColumns);
//Rprintf("    baseSeed = %d\n",myControl.baseSeed);
//Rprintf("    rank = %d\n",rank);
//Rprintf("    chain = %d\n",chain);
//R_FlushConsole();
//R_ProcessEvents();
	pattern = new HMatrix(rank, myData.nColumns, Utility::GetSeed(myControl.baseSeed, rank, chain, 'P'));
//Rprintf("    Update::GetInitial(): Got this far #55...\n");
//R_FlushConsole();
//R_ProcessEvents();
	amplitude = new WMatrix(myData.nRows, rank, Utility::GetSeed(myControl.baseSeed, rank, chain, 'A'));
//Rprintf("    Update::GetInitial(): Got this far #60...\n");
//R_FlushConsole();
//R_ProcessEvents();
	mock = new VMatrix(myData.nRows, myData.nColumns);
	weightMatrix = new WeightingMatrix(myData.nRows, myData.nColumns);
//Rprintf("    Update::GetInitial(): Got this far #100...\n");
//R_FlushConsole();
//R_ProcessEvents();

	ostringstream outS;
	outS << myControl.rootName << "." << myControl.scheme << ".k" << rank;
	// outS << "/scratch/rootName." << myControl.scheme << ".k" << rank;
	if(myControl.scheme == "Renyi") outS << "_alpha" << alpha;
	outS << "_run" << chain;
	
	wMatrixName = outS.str() + ".W";
	hMatrixName = outS.str() + ".H";

	outS.clear();
//Rprintf("    Update::GetInitial(): Got this far #200...\n");
//R_FlushConsole();
//R_ProcessEvents();

	if(myControl.theta != 0.0) BuildThetaMatrix();

	// assign row and column names
	vector<string> metaGenes(pattern->nRows);
	for (int i=0; i<pattern->nRows; i++)
	{
		ostringstream columnName;
		columnName << "cluster_" << i+1;
        metaGenes[i] = columnName.str();
	}
//Rprintf("    Update::GetInitial(): Got this far #300...\n");
//R_FlushConsole();
//R_ProcessEvents();

	amplitude->rowNames = myData.rowNames;
	amplitude->columnNames = metaGenes;
	pattern->rowNames = metaGenes;
	pattern->columnNames = myData.columnNames;

	// try to bring mock as closer as possible to the original data
	if(myControl.InitialScaling)
	{
		GenerateMock();
		// myData.DScaling( mock->data );
		mock->DScaling(myData.data);
		mock->GetStatistical();
		amplitude->Scaling( 1.0 / mock->average );
	}
//Rprintf("    Update::GetInitial(): Got this far #400...\n");
//R_FlushConsole();
//R_ProcessEvents();

	// prepare the initial matrixes
	Normalizing();
	GenerateMock();
//Rprintf("    Update::GetInitial(): Got this far #500...\n");
//R_FlushConsole();
//R_ProcessEvents();
}

void Update::GenerateMock()
{
	if(myControl.theta == 0.0)
	{
		amplitude->Multipling(pattern->data, mock->data);
	}
	else
	{
		amplitude->Multipling(pattern->data, thetaMatrix->data, mock->data);
	}
}

void Update::GenerateMock(const WMatrix &w, const HMatrix &h)
{
	if(myControl.theta == 0.0)
	{
		w.Multipling(h.data, mock->data);
	}
	else
	{
		w.Multipling(h.data, thetaMatrix->data, mock->data);
	}
}

void Update::CalculateWeights()
{
	// myData.Dividing( mock->data, weightMatrix->data );
	// save some memory
	if(myData.HD)
	{
		for (unsigned i=0; i<myData.hdIndex.size(); i++)
		{
			for (unsigned j=0; j<myData.hdIndex[i].size(); j++)
			{
				mock->data[i][(int) myData.hdIndex[i][j]] = myData.hdValue[i][j] / mock->data[i][(int) myData.hdIndex[i][j]];
			}
		}
	}
	else
	{
		for (unsigned i=0; i<myData.data.size(); i++)
		{
			for (unsigned j=0; j<myData.data[i].size(); j++)
			{
				mock->data[i][j] = myData.data[i][j] / mock->data[i][j];
			}
		}
	}
}

int Update::CheckConvergency(int round)
{
	if( round % myControl.convergeTestStep != 0)
	{
		converged = false;
		return 0;
	}
	else
	{
		errors = myData.KL_Div( mock->data );
		consecutiveError = fabs(errors - myData.KL_Div( prevMock ));
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

void Update::IndMatrixEvaluation()
{
}

void Update::BuildThetaMatrix()
{
	thetaMatrix = new WMatrix(rank, rank);
	double offDia = myControl.theta / rank;
	for (unsigned i=0; i<(unsigned)thetaMatrix->nRows; i++)
	{
		thetaMatrix->data[i][i] = 1.0 - myControl.theta + offDia; 
		for (unsigned j=i+1; j<(unsigned)thetaMatrix->nColumns; j++)
		{
			thetaMatrix->data[i][j] = offDia;
			thetaMatrix->data[j][i] = offDia;
		}
	}
}

