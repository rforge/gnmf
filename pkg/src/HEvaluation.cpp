#include "HEvaluation.h"

// JMM (7/10/2016): Define R_NO_REMAP
// Reference: First few paragraphs in in Section 6 The R API: entry points for C code
// in the online documentation for "Writing R Extensions".
// Also see email from Prof. B. Ripley sent 7/10/2016
// at 3:58 AM entitled "CRAN packages failing to install with modern C++"
#define R_NO_REMAP
#include "R.h"

using namespace std;

void HEvaluation::ProcessingFiles()
{
Rprintf("HEvaluation::ProcessingFiles(): entered function...\n"); R_FlushConsole(); R_ProcessEvents();
Rprintf("HEvaluation::ProcessingFiles(): myControl.RunType = %s\n",myControl.RunType.data()); R_FlushConsole(); R_ProcessEvents();
	if(myControl.RunType == "evaluation")
	{
Rprintf("HEvaluation::ProcessingFiles(): Got here #10\n"); R_FlushConsole(); R_ProcessEvents();
		// cout << "Begin processing File" << endl;
		for(int chain=1; chain<=myControl.nChains; chain++)
		{	
			// cout << "chain " << chain << endl;
			ostringstream targetFile;
			targetFile << factorMatrixRootName << "_run" << chain << ".H";
			// cout << targetFile.str() << endl;
			ifstream IN;
//Rprintf("HEvaluation::ProcessingFiles(): targetFile.str() = %s\n",targetFile.str().data()); R_FlushConsole(); R_ProcessEvents();
			IN.open(targetFile.str().c_str(), ios::in);
Rprintf("HEvaluation::ProcessingFiles(): Got here #15\n"); R_FlushConsole(); R_ProcessEvents();
//			if(IN.fail())
//			{
Rprintf("HEvaluation::ProcessingFiles(): Got here #16\n"); R_FlushConsole(); R_ProcessEvents();
				// cout << "Read assembled file" << endl;
//				ReadinAssembledFiles();
Rprintf("HEvaluation::ProcessingFiles(): Got here #17\n"); R_FlushConsole(); R_ProcessEvents();
//				return;
//			}
Rprintf("HEvaluation::ProcessingFiles(): Got here #18\n"); R_FlushConsole(); R_ProcessEvents();
	    		HMatrix *pm = new HMatrix(nRows, nColumns);
			// cout << "allocate mem for pm " << nRows << " " << nColumns << endl;
			IN >> *pm;
			IN.clear();
			IN.close();


			// delete original file
//			remove( targetFile.str().c_str() );
//			targetFile.clear();
			
			targetMatrixes.push_back( pm );
Rprintf("HEvaluation::ProcessingFiles(): Got here #30\n"); R_FlushConsole(); R_ProcessEvents();

//			CopySaveProperties( pm->properties, chain );
            int numRanks    = myControl.endRank - myControl.startRank + 1;
            int numChains   = myControl.nChains;
            int alphaOffset = numRanks * numChains;
            int rankIndex   = rank - myControl.startRank;
			int j = ( myControl.alphaIndex * alphaOffset ) + ( rankIndex * numChains ) + ( chain - 1 );
Rprintf("HEvaluation::ProcessingFiles(): j = %d\n",j); R_FlushConsole(); R_ProcessEvents();
			SSMap::const_iterator ip      = pm->properties.begin();
			myControl.consecutiveError[j] = 2.6; ip++;
			myControl.converged[j]        = atof(ip->second.data()); ip++;
			myControl.error[j]            = atof(ip->second.data()); ip++;
			myControl.Hsparseness[j]      = atof(ip->second.data()); ip++;
			myControl.totalSteps[j]       = atof(ip->second.data()); ip++;
			myControl.Wsparseness[j]      = atof(ip->second.data());
Rprintf("HEvaluation::ProcessingFiles(): Got here #58\n"); R_FlushConsole(); R_ProcessEvents();
		}
Rprintf("HEvaluation::ProcessingFiles(): Got here #59\n"); R_FlushConsole(); R_ProcessEvents();
	}
	else
	{
Rprintf("HEvaluation::ProcessingFiles(): Got here #60\n"); R_FlushConsole(); R_ProcessEvents();
        int numRanks    = myControl.endRank - myControl.startRank + 1;
        int numChains   = myControl.nChains;
        int alphaOffset = numRanks * numChains;
        int rankIndex   = rank - myControl.startRank;
		for (unsigned i=0; i<targetMatrixes.size(); i++)
		{
//			CopySaveProperties( targetMatrixes[i]->properties, i );
			int j = ( myControl.alphaIndex * alphaOffset ) + ( rankIndex * numChains ) + i;
Rprintf("HEvaluation::ProcessingFiles(): j = %d\n",j); R_FlushConsole(); R_ProcessEvents();
			SSMap::const_iterator ip      = targetMatrixes[i]->properties.begin();
Rprintf("HEvaluation::ProcessingFiles(): Got here #70\n"); R_FlushConsole(); R_ProcessEvents();
			myControl.consecutiveError[j] = 2.6; ip++;
Rprintf("HEvaluation::ProcessingFiles(): Got here #71\n"); R_FlushConsole(); R_ProcessEvents();
			myControl.converged[j]        = atof(ip->second.data()); ip++;
Rprintf("HEvaluation::ProcessingFiles(): Got here #72\n"); R_FlushConsole(); R_ProcessEvents();
			myControl.error[j]            = atof(ip->second.data()); ip++;
Rprintf("HEvaluation::ProcessingFiles(): Got here #73\n"); R_FlushConsole(); R_ProcessEvents();
			myControl.Hsparseness[j]      = atof(ip->second.data()); ip++;
Rprintf("HEvaluation::ProcessingFiles(): Got here #74\n"); R_FlushConsole(); R_ProcessEvents();
			myControl.totalSteps[j]       = atof(ip->second.data()); ip++;
Rprintf("HEvaluation::ProcessingFiles(): Got here #75\n"); R_FlushConsole(); R_ProcessEvents();
			myControl.Wsparseness[j]      = atof(ip->second.data());
Rprintf("HEvaluation::ProcessingFiles(): Got here #76\n"); R_FlushConsole(); R_ProcessEvents();
		}
	}
Rprintf("HEvaluation::ProcessingFiles(): Got here #80\n"); R_FlushConsole(); R_ProcessEvents();
	// order matrixes make all runs comparable
	OrderMatrixes();
	CombineMatrixes();
}
void HEvaluation::GetSigmas()
{
	average.resize(nRows, array1d(nColumns));
	sigma.resize(nRows, array1d(nColumns));

	for(int chain=0; chain<myControl.nChains; chain++)
	{	
		for (int i=0; i<nRows; i++)
		{
			for (int j=0; j<nColumns; j++)
			{
				average[i][j] += targetMatrixes[chain]->data[i][j];
				sigma[i][j] += targetMatrixes[chain]->data[i][j] * targetMatrixes[chain]->data[i][j];
			}
		}
	}

	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			sigma[i][j] = (sigma[i][j] - average[i][j] * average[i][j] / myControl.nChains) / (myControl.nChains - 1);
			average[i][j] /= myControl.nChains;
		}
	}
}

void HEvaluation::CombineMatrixes()
{
	// vector< WMatrix * > amplitudes;
	if( myControl.RunType == "evaluation" )
	{
		for (int i=0; i<myControl.nChains; i++)
		{
			ostringstream targetFile;
			targetFile << factorMatrixRootName << "_run" << i + 1 << ".W";
			ifstream IN(targetFile.str().c_str(), ios::in);

	    		WMatrix *am = new WMatrix(myControl.nObservations, rank);
			IN >> *am;
			IN.close();
			
			// delete original file
			remove( targetFile.str().c_str() );
			targetFile.clear();

			am->Transposing();
			amplitudes.push_back(am);
		}
	}

	// re-ordering
	for (int i=1; i<myControl.nChains; i++)
	{
		array2d newAData;
		array2d newPData;
		for (int j=0; j<rank; j++)
		{
			newPData.push_back( targetMatrixes[i]->data[ maps4chains[i-1][j] ] );
			newAData.push_back( amplitudes[i]->data[ maps4chains[i-1][j] ] );
		}
		targetMatrixes[i]->data = newPData;
		amplitudes[i]->data = newAData;
	}

	// integrate W and H from all runs into single files
	for (int i=0; i<rank; i++)
	{
//		ostringstream wf;
//		ostringstream hf;

//		wf << factorMatrixRootName << ".W" << i;
//		hf << factorMatrixRootName << ".H" << i;

//		ofstream wO(wf.str().c_str(), ios::out);
//		ofstream hO(hf.str().c_str(), ios::out);

		DataMatrix combinedW(myControl.nChains, myControl.nObservations);
		DataMatrix combinedH(myControl.nChains, myControl.nParameters);
		combinedH.columnNames = targetMatrixes[0]->columnNames;
		combinedW.columnNames = amplitudes[0]->columnNames;
		for (int j=0; j<myControl.nChains; j++)
		{
			ostringstream temp;
			temp << "Run" << j;
			combinedW.rowNames[j] = temp.str();
			combinedH.rowNames[j] = combinedW.rowNames[j];
		}

		// first run's data, they are the based order
		combinedH.data[0] = targetMatrixes[0]->data[i];
		combinedW.data[0] = amplitudes[0]->data[i];

		for (int j=1; j<myControl.nChains; j++)
		{
			combinedH.data[j] = targetMatrixes[j]->data[i];
			combinedW.data[j] = amplitudes[j]->data[i];
		}
Rprintf("        HEvaluation: i = %d\n",i);
Rprintf("        combinedW is %d by %d\n",combinedW.nRows,combinedW.nColumns);
Rprintf("        combinedH is %d by %d\n",combinedH.nRows,combinedH.nColumns);
Rprintf("        combinedW.data[2][3] = %g\n",combinedW.data[2][3]);
Rprintf("        combinedH.data[2][3] = %g\n",combinedH.data[2][3]);
Rprintf("        alphaIndex is %d\n",myControl.alphaIndex);
R_FlushConsole();
R_ProcessEvents();

//		wO << combinedW;
//		hO << combinedH;

        int sizeAlphaBlock = 0;
        for (int j=myControl.startRank; j<=myControl.endRank; j++ )
            sizeAlphaBlock += j;
        sizeAlphaBlock *= myControl.nChains;
        int alphaOffset = myControl.alphaIndex * sizeAlphaBlock;
        int rankOffset = 0;
        for (int j=myControl.startRank; j<=rank-1; j++ )
            rankOffset += ( j * myControl.nChains );
Rprintf("        sizeAlphaBlock is %d\n",sizeAlphaBlock);
Rprintf("        startRank      is %d\n",myControl.startRank);
Rprintf("        endRank        is %d\n",myControl.endRank);
Rprintf("        rank           is %d\n",rank);
Rprintf("        rankOffset     is %d\n",rankOffset);

        // Transfer data from combinedH to R matrix H.
Rprintf("About to transfer data to H matrix...\n");
R_FlushConsole();
R_ProcessEvents();
        int rowIndex, colIndex;
        for (int j=0; j<combinedH.nRows; j++)
            for (int k=0; k<combinedH.nColumns; k++)
            {
                rowIndex = alphaOffset + rankOffset + ( j * myControl.endRank ) + i;
                colIndex = k;
//Rprintf("            [ %d, %d ]\n",rowIndex,colIndex);
                myControl.H[rowIndex*combinedH.nColumns + colIndex] = combinedH.data[j][k];
            }
            
        // Transfer data from combinedH to R matrix W.
Rprintf("About to transfer data to W matrix...\n");
R_FlushConsole();
R_ProcessEvents();
        for (int j=0; j<combinedW.nRows; j++)
            for (int k=0; k<combinedW.nColumns; k++)
            {
                rowIndex = alphaOffset + rankOffset + ( j * myControl.endRank ) + i;
                colIndex = k;
Rprintf("            [ %d, %d ]\n",rowIndex,colIndex);
                myControl.W[rowIndex*combinedW.nColumns + colIndex] = combinedW.data[j][k];
            }

//		wO.close();
//		hO.close();
	}

	// for (unsigned i=0; i<amplitudes.size(); i++) delete amplitudes[i];
}

void HEvaluation::ReadinAssembledFiles()
{
	// targetMatrixes.resize( myControl.nChains );
	for (int i=0; i<myControl.nChains; i++)
	{
		HMatrix *pm = new HMatrix(rank, myControl.nParameters);
		// pm->data.resize(rank, array1d(myControl.nParameters));
		// pm->columnNames.resize(myControl.nParameters);

		targetMatrixes.push_back( pm );
	}

	for (int i=0; i<rank; i++)
	{
		ostringstream hf;
		hf << factorMatrixRootName << ".H" << i;
		ifstream IN(hf.str().c_str(), ios::in);

		DataMatrix assembled(myControl.nChains, myControl.nParameters);
		IN >> assembled;

		for (int j=0; j<myControl.nChains; j++)
		{
			targetMatrixes[j]->data[i] = assembled.data[j];
			targetMatrixes[j]->columnNames = assembled.columnNames;
		}

		IN.close();
	}

}

void HEvaluation::SetParameters()
{
	fileType = ".H";
	nRows = rank;
	nColumns = myControl.nParameters;
	// cout << "Allocate memory for cMatrix " << nColumns << endl;
	ccMatrix = new CMatrix(nColumns);

	// open log file for recording message
	ostringstream logFile;
	logFile << factorMatrixLogRootName << fileType << "log";
	outLOG.open(logFile.str().c_str(), ios::out);
	// cout << "#Simulation results for rank(K)=" << rank << " alpha=" << alpha << " " << logFile.str() << endl;
	/*
	logFile << ".xml";
	xmlOutLog = xmlNewTextWriterFilename( logFile.str().c_str(), 0 );
	if(xmlOutLog == NULL)
	{
		std::cerr << "textXmlWriterFilename: Error creating the xml writer!" << endl;
	}

	xmlTextWriterStartDocument(xmlOutLog, NULL, NULL, NULL);

	xmlTextWriterStartElement(xmlOutLog, BAD_CAST "Simulation");
	xmlTextWriterWriteFormatAttribute(xmlOutLog, BAD_CAST "rank", "%d", rank);
	xmlTextWriterWriteFormatAttribute(xmlOutLog, BAD_CAST "alpha", "%5.3f", alpha);
	*/
	logFile.clear();

}

