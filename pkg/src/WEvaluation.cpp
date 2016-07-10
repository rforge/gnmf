#include "WEvaluation.h"

// JMM (7/10/2016): Do not include R include files in GNMF header files, because sometimes
// GNMF header files include other GNMF header files. This makes it difficult to ensure
// that R header files are included after GNMF header files. It is probably best to just
// include R header files directly in the CPP files where they are used.
#define R_NO_REMAP
#include "R.h"           // R functions

using namespace std;

void WEvaluation::ProcessingFiles()
{
Rprintf("WEvaluation::ProcessingFiles(): entered function...\n"); R_FlushConsole(); R_ProcessEvents();
Rprintf("WEvaluation::ProcessingFiles(): myControl.RunType = %s\n",myControl.RunType.data()); R_FlushConsole(); R_ProcessEvents();
	if(myControl.RunType == "evaluation") // or if(targetMatrixes.size() == 0
	{
Rprintf("WEvaluation::ProcessingFiles(): Got here #10\n"); R_FlushConsole(); R_ProcessEvents();
		for(int chain=1; chain<=myControl.nChains; chain++)
		{	
			ostringstream targetFile;
			targetFile << factorMatrixRootName << "_run" << chain << ".W";
			ifstream IN;
			IN.open(targetFile.str().c_str(), ios::in);
			if(IN.fail())
			{
				ReadinAssembledFiles();
				return;
			}

			HMatrix *am = new HMatrix(nRows, nColumns);
			IN >> *am;
			IN.close();

			// delete original file
			remove( targetFile.str().c_str() );
			targetFile.clear();
			
			am->Transposing();
			targetMatrixes.push_back( am );
Rprintf("WEvaluation::ProcessingFiles(): Got here #30\n"); R_FlushConsole(); R_ProcessEvents();

			CopySaveProperties( am->properties, chain );
            int numRanks    = myControl.endRank - myControl.startRank + 1;
            int numChains   = myControl.nChains;
            int alphaOffset = numRanks * numChains;
            int rankIndex   = rank - myControl.startRank;
			int j = ( myControl.alphaIndex * alphaOffset ) + ( rankIndex * numChains ) + ( chain - 1 );
Rprintf("WEvaluation::ProcessingFiles(): j = %d\n",j); R_FlushConsole(); R_ProcessEvents();
			SSMap::const_iterator ip      = am->properties.begin();
			myControl.consecutiveError[j] = 2.6; ip++;
			myControl.converged[j]        = atof(ip->second.data()); ip++;
			myControl.error[j]            = atof(ip->second.data()); ip++;
			myControl.Hsparseness[j]      = atof(ip->second.data()); ip++;
			myControl.totalSteps[j]       = atof(ip->second.data()); ip++;
			myControl.Wsparseness[j]      = atof(ip->second.data());
		}
	}
	else
	{
Rprintf("WEvaluation::ProcessingFiles(): Got here #60\n"); R_FlushConsole(); R_ProcessEvents();
        int numRanks    = myControl.endRank - myControl.startRank + 1;
        int numChains   = myControl.nChains;
        int alphaOffset = numRanks * numChains;
        int rankIndex   = rank - myControl.startRank;
		for (unsigned i=0; i<targetMatrixes.size(); i++)
		{
			CopySaveProperties( targetMatrixes[i]->properties, i );
			int j = ( myControl.alphaIndex * alphaOffset ) + ( rankIndex * numChains ) + i;
Rprintf("WEvaluation::ProcessingFiles(): j = %d\n",j); R_FlushConsole(); R_ProcessEvents();
			SSMap::const_iterator ip      = targetMatrixes[i]->properties.begin();
			myControl.consecutiveError[j] = 2.6; ip++;
			myControl.converged[j]        = atof(ip->second.data()); ip++;
			myControl.error[j]            = atof(ip->second.data()); ip++;
			myControl.Hsparseness[j]      = atof(ip->second.data()); ip++;
			myControl.totalSteps[j]       = atof(ip->second.data()); ip++;
			myControl.Wsparseness[j]      = atof(ip->second.data());
		}
	}

	// order matrixes make all runs comparable
	OrderMatrixes();
	CombineMatrixes();
}

void WEvaluation::GetSigmas()
{
	// am in targetMatrixes is transposed
	average.resize(nColumns, array1d(nRows));
	sigma.resize(nColumns, array1d(nRows));

	for(int chain=0; chain<myControl.nChains; chain++)
	{	
		for (int i=0; i<nRows; i++)
		{
			for (int j=0; j<nColumns; j++)
			{
				average[i][j] += targetMatrixes[chain]->data[j][i]; // transposed before
				sigma[i][j] += targetMatrixes[chain]->data[j][i] * targetMatrixes[chain]->data[j][i];
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


void WEvaluation::CombineMatrixes()
{
	// vector< HMatrix *> patterns;
	if( myControl.RunType == "evaluation" )
	{
		for (int i=0; i<myControl.nChains; i++)
		{
			ostringstream targetFile;
			targetFile << factorMatrixRootName << "_run" << i + 1 << ".H";
			ifstream IN(targetFile.str().c_str(), ios::in);

	    	HMatrix *pm = new HMatrix(rank, myControl.nParameters);
			IN >> *pm;
			IN.close();
			
			// delete original file
			remove( targetFile.str().c_str() );
			targetFile.clear();

			patterns.push_back(pm);
		}
	}

	// re-ordering
	for (int i=1; i<myControl.nChains; i++)
	{
		array2d newAData;
		array2d newPData;
		for (int j=0; j<rank; j++)
		{
			newAData.push_back( targetMatrixes[i]->data[ maps4chains[i-1][j] ] );
			newPData.push_back( patterns[i]->data[ maps4chains[i-1][j] ] );
		}
		targetMatrixes[i]->data = newAData;
		patterns[i]->data = newPData;
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
		combinedW.columnNames = targetMatrixes[0]->columnNames; // already transposed
		combinedH.columnNames = patterns[0]->columnNames;
		for (int j=0; j<myControl.nChains; j++)
		{
			ostringstream temp;
			temp << "Run" << j;
			combinedW.rowNames[j] = temp.str();
			combinedH.rowNames[j] = combinedW.rowNames[j];
		}

		// first run's data, they are the based order
		combinedW.data[0] = targetMatrixes[0]->data[i]; // after transposed
		combinedH.data[0] = patterns[0]->data[i];

		for (int j=1; j<myControl.nChains; j++)
		{
			combinedW.data[j] = targetMatrixes[j]->data[i]; // after transposed
			combinedH.data[j] = patterns[j]->data[i];
		}
Rprintf("        WEvaluation: i = %d\n",i);
Rprintf("        combinedW is %d by %d\n",combinedW.nRows,combinedW.nColumns);
Rprintf("        combinedH is %d by %d\n",combinedH.nRows,combinedH.nColumns);
Rprintf("        combinedW.data[2][3] = %g\n",combinedW.data[2][3]);
Rprintf("        combinedH.data[2][3] = %g\n",combinedH.data[2][3]);
Rprintf("        alphaIndex is %d\n",myControl.alphaIndex);

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

	/* do NOT need to, because transposed one is easier to process 
	// transpose back A matrixes
	for (int i=1; i<myControl.nChains; i++)
	{
		amplitudes[i].Transposing();
	}
	*/
}

void WEvaluation::ReadinAssembledFiles()
{
	// targetMatrixes.resize( myControl.nChains );
	for (int i=0; i<myControl.nChains; i++)
	{
		HMatrix *am = new HMatrix(rank, myControl.nObservations);
		// am->data.resize(rank, array1d(myControl.nObservations));
		targetMatrixes.push_back( am );
	}

	for (int i=0; i<rank; i++)
	{
		ostringstream af;
		af << factorMatrixRootName << ".W" << i;
		ifstream IN(af.str().c_str(), ios::in);

		DataMatrix assembled(myControl.nChains, myControl.nObservations);
		IN >> assembled;

		for (int j=0; j<myControl.nChains; j++)
		{
			targetMatrixes[j]->data[i] = assembled.data[j];
		}

		IN.close();
	}
}

void WEvaluation::SetParameters()
{
	fileType = ".W";
	nRows = myControl.nObservations;
	nColumns = rank;
	ccMatrix = new CMatrix(nRows);

	// open log file for recording message
	ostringstream logFile;
	logFile << factorMatrixLogRootName << fileType << "log";
	outLOG.open(logFile.str().c_str(), ios::out);

	// outLOG << "#Simulation results for rank(K)=" << rank << " alpha=" << alpha << endl;
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


