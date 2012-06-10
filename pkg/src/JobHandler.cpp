#include "JobHandler.h"

using namespace std;



JobHandler::JobHandler(double *matrixX,
                       int    *numRows,
                       int    *numCols,
                       char   **scheme,
                       int    *nsteps,
                       int    *repeats,
                       int    *rankrange,
                       int    *numRankRange,
                       char   **cltarget,
                       char   **clscheme,
                       char   **reffile,
                       char   **scaling,
                       char   **normalizing,
                       double *alphas,
                       int    *nalphas,
                       char   **runtype,
                       int    *cstepsize,
                       double *idealization,
                       double *matrixH,
                       double *matrixW,
                       double *converged,
                       double *totalSteps,
                       double *error,
                       double *consecutiveError,
                       double *Wsparseness,
                       double *Hsparseness)
{
//Rprintf("JobHandler::JobHandler(): entered function...\n");
//R_FlushConsole();
//R_ProcessEvents();

	myControl = new ParaControler();

myControl->rootName = "gnmf";

    // Parse rankrange parameter.
    myControl->startRank = rankrange[0];
	if( *numRankRange == 1 )
	{
	    myControl->endRank = rankrange[0];
	} else {
	    myControl->endRank = rankrange[1];
    }

    // Load alpha values, possibly multiple.
    myControl->alphaSet.clear();
    for (int i=0; i< *nalphas; i++)
        myControl->alphaSet.push_back(alphas[i]);
//    myControl->nAlphas = *nalphas;
    myControl->nAlphas = myControl->alphaSet.size();
//Rprintf("JobHandler::JobHandler(): Got this far #100...\n");
//R_FlushConsole();
//R_ProcessEvents();

    // Load other values.
	myControl->scheme         = scheme[0];
	myControl->nUpdateSteps   = *nsteps;
	myControl->nChains        = *repeats;
	myControl->idealization   = *idealization;    
    myControl->target         = cltarget[0];
    myControl->InitialScaling = ( scaling[0][0] == 'T' ) ? true:false;
    myControl->refClusterFile = reffile[0];
	myControl->mcFlag         = false; // Turn ref file off for now
	myControl->Normalize      = ( normalizing[0][0] == 'T' ) ? true:false;
	myControl->RunType        = runtype[0];
//Rprintf("JobHandler::JobHandler(): runtype = %s...\n",myControl->RunType.data());
//R_FlushConsole();
//R_ProcessEvents();
	myControl->HDFile         = false; // Turn HD file off for now
	myControl->clustScheme    = clscheme[0];
//  myControl->uncertFile     = value;

	thisrank = -1;

	myData = new VMatrix();
    myData->data.resize(*numRows);
    for (int i=0; i<*numRows; i++)
        myData->data[i].resize(*numCols);
    myData->nRows    = *numRows;
    myData->nColumns = *numCols;
    // Now myData is numRows rows by numCols columns.

	myData->HD = myControl->HDFile;
//Rprintf("JobHandler::JobHandler(): Got this far #200...\n");
//Rprintf("myData is %d by %d\n",myData->nRows,myData->nColumns);
//R_FlushConsole();
//R_ProcessEvents();

	// TRANSFER INPUT DATA FROM R TO MYDATA.
    for (int i=0; i< *numRows; i++)
        for (int j=0; j< *numCols; j++)
        {
	        myData->data[i][j] = matrixX[(j*(*numRows))+i];
        }

//Rprintf("JobHandler::JobHandler(): Got this far #240...\n");
//R_FlushConsole();
//R_ProcessEvents();
	myControl->nObservations = myData->nRows;
	myControl->nParameters   = myData->nColumns;
//Rprintf("JobHandler::JobHandler(): Got this far #250...\n");
//R_FlushConsole();
//R_ProcessEvents();

	// myData->Idealize();
	// For sparse matrix
	myData->Idealize( myControl->idealization );
//Rprintf("JobHandler::JobHandler(): Got this far #275...\n");
//R_FlushConsole();
//R_ProcessEvents();

	// uncertainty is only used for LSNMF
	if(myControl->uncertFile != "")
	{
		mySTD = new VMatrix();
		mySTD->HD = myControl->HDFile;
		ifstream iUncert(myControl->uncertFile.c_str(), ios::in);
		iUncert >> *mySTD;
		iUncert.close();

		mySTD->Idealize( 0.1 );

		// scale data according to the uncertainty
		myData->DScaling(mySTD->data);
		// myControl->scheme = "KL";
	}
	else
	{
		mySTD = 0;
	}
//Rprintf("JobHandler::JobHandler(): Got this far #300...\n");
//R_FlushConsole();
//R_ProcessEvents();

	// reference clust data
	if(myControl->mcFlag)
	{
		refClustData = new ClustDataSet();
		ifstream iRef(myControl->refClusterFile.c_str(), ios::in);
		iRef >> *refClustData;
		iRef.close();
	}
	else
	{
		refClustData = 0;
	}
	
	// Pointers to matrices in R memory.
//Rprintf("JobHandler::JobHandler: setting up pointers...\n"); R_FlushConsole(); R_ProcessEvents();
	myControl->H                = matrixH;
    myControl->W                = matrixW;
    myControl->converged        = converged;
    myControl->totalSteps       = totalSteps;
    myControl->error            = error;
    myControl->consecutiveError = consecutiveError;
    myControl->Wsparseness      = Wsparseness;
    myControl->Hsparseness      = Hsparseness;
    
//Rprintf("JobHandler::JobHandler(): exiting function...\n");
//R_FlushConsole();
//R_ProcessEvents();

}


void JobHandler::VerifyArgumentList()
{
	jobStatus = 1;
	if(myCpuIndex == 0)
	{
		if(nArguments < 3)
		{
			PrintUsage();
			jobStatus = -1;
		}
	}

#ifdef MPI_PARALLEL
	MPI_Bcast(&jobStatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(jobStatus == -1) Finalizing();
#else
	if(jobStatus == -1)
//    exit(0);
    ;
#endif
}

void JobHandler::PrintUsage()
{
Rprintf("JobHandler::PrintUsage()\n");
R_FlushConsole();
R_ProcessEvents();
/*
	cout << "Usage: GNMF -i input_data_matrix" << endl;
	cout << endl;
	cout << "General suggestion: For an input matrix containing no zero entries, it's recommended that idealization be set to 0 (default)," << endl;
	cout << "  scaling to 'T' and normalizing to 'T'. For an input matrix containing zero entries, idealization should be set to 0.1," << endl;
	cout << "  scaling to 'F' and normalizing to 'F'. If idealization is set to 0.1, setting both scaling and normalizing to 'T' " << endl;
	cout << "  is recommended for faster convergence" << endl;
	cout << endl;
	cout << "Other options:" << endl;
	cout << "  -s: scheme, 'KL', 'Renyi', 'ED', 'GammaJD', 'GammaKL', 'DivComb', 'Div2', 'Pareto', 'InverseLink_GammaKL', 'BD', 'NBD', 'ODP'" << endl;
	// cout << "  -s: scheme, 'KL', 'Renyi', 'ED'" << endl;
	cout << "  -n: update steps, default 2000" << endl;
	cout << "  -c: repeats, default 20" << endl;
	cout << "  -k: rank range, default 2-2" << endl;
	cout << "  -t: clustering target, default 'PATTERN (H matrix)'" << endl;
	cout << "  -cs: clustering scheme, default 'Binary', could be 'PearsonHC'" << endl;
	cout << "  -r: reference file, string" << endl;
	cout << "  -scaling: initial scaling for faster convergence, default 'F'" << endl;
	cout << "  -normalizing: H matrix normalization, default 'F'" << endl;
	cout << "  -alphas: csv string, valid only for scheme Renyi and Pareto, default 1.0" << endl;
	// cout << "  -alphas: csv string, valid only for scheme Renyi, default 1.0" << endl;
	cout << "  -RunType: default 'whole', could be 'simulation' or 'evaluation'" << endl;
	// cout << "  -std_file: expression variance file, can be used to weight primary input matrix elements" << endl; 
	// cout << "  -hd: sparse matrix iput? default 'F'" << endl;
	cout << "  -cts: set convergence test step size, default 20" << endl;
	cout << "  -idealization, default 0, can be set to 0.1" << endl;
	*/
}

void JobHandler::GetInput()
{
	myControl = new ParaControler();

	for (int i=1; i<nArguments; i+=2)
	{
		string option = string(argumentList[i]);
		string value = string(argumentList[i + 1]);
		if(option == "-i")
		{
			myControl->sourceFile = value;
			myControl->rootName = value.substr(0, value.rfind('.'));
		}
		else if(option == "-s")
		{
			myControl->scheme = value;
		}
		else if(option == "-n")
		{
			myControl->nUpdateSteps = atoi( value.c_str() );
		}
		else if(option == "-cts")
		{
			myControl->convergeTestStep = atoi( value.c_str() );
		}
		else if(option == "-c")
		{
			myControl->nChains = atoi( value.c_str() );
		}
		else if(option == "-idealization")
		{
			myControl->idealization = atof( value.c_str() );
		}
		else if(option == "-k")
		{
			string::size_type loc = value.find('-');
			if(loc == string::npos)
			{
//				cout << "Wrong format for -k!" << endl;
				if(myCpuIndex == 0) PrintUsage();
//				exit(0);
			}
			myControl->startRank = atoi( value.substr(0, loc).c_str() );
			myControl->endRank = atoi( value.substr(loc+1).c_str() );
		}
		else if(option == "-alphas")
		{
			char alpha[20];
			stringstream alphas;
			alphas << value;
			myControl->alphaSet.clear();
			while( alphas.getline(alpha, 20, ',') )
			{
				myControl->alphaSet.push_back( atof(alpha) );
			}
			myControl->nAlphas = myControl->alphaSet.size();
		}
		else if(option == "-t")
		{
			myControl->target = value;
		}
		else if(option == "-r")
		{
			myControl->refClusterFile = value;
			myControl->mcFlag = true;
		}
		else if(option == "-scaling")
		{
			myControl->InitialScaling = ( value == "T" ) ? true:false;
		}
		else if(option == "-normalizing")
		{
			myControl->Normalize = ( value == "T" ) ? true:false;
		}
		else if(option == "-RunType")
		{
			myControl->RunType = value;
		}
		else if(option == "-hd")
		{
			myControl->HDFile = ( value == "T" ) ? true:false;
		}
		else if(option == "-cs")
		{
			myControl->clustScheme = value;
		}
		else if(option == "-std_file")
		{
			myControl->uncertFile = value;
		}
		else
		{
			if(myCpuIndex == 0) PrintUsage();
//			exit(0);
		}
			
	}

	thisrank = -1;

	myData = new VMatrix();
	myData->HD = myControl->HDFile;
	ifstream iData(myControl->sourceFile.c_str(), ios::in);
	if(! iData )
	{
Rprintf("JobHandler::GetInput(): myControl->sourceFile is not opened\n");
R_FlushConsole();
R_ProcessEvents();
//		cout << myControl->sourceFile << " is not opened" << endl;
	}
	
	iData >> *myData;
	iData.close();

	myControl->nObservations = myData->nRows;
	myControl->nParameters = myData->nColumns;

	// myData->Idealize();
	// For sparse matrix
	myData->Idealize( myControl->idealization );

	// uncertainty is only used for LSNMF
	if(myControl->uncertFile != "")
	{
		mySTD = new VMatrix();
		mySTD->HD = myControl->HDFile;
		ifstream iUncert(myControl->uncertFile.c_str(), ios::in);
		iUncert >> *mySTD;
		iUncert.close();

		mySTD->Idealize( 0.1 );

		// scale data according to the uncertainty
		myData->DScaling(mySTD->data);
		// myControl->scheme = "KL";
	}
	else
	{
		mySTD = 0;
	}

	// reference clust data
	if(myControl->mcFlag)
	{
		refClustData = new ClustDataSet();
		ifstream iRef(myControl->refClusterFile.c_str(), ios::in);
		iRef >> *refClustData;
		iRef.close();
	}
	else
	{
		refClustData = 0;
	}
}

void JobHandler::PrintInput() const
{
	// only master node print out checkpoint file
	if(myCpuIndex > 0) return;

	string checkFile = myControl->rootName + ".chk";
	ofstream checkOUT(checkFile.c_str(), ios::out);

	checkOUT << "Control parameters:" << endl;
	checkOUT << *myControl;
	checkOUT << "\nOriginal input data:" << endl;
	checkOUT << *myData;

	if(mySTD != 0)
	{
		checkOUT << "\nUncertainty input data:" << endl;
		checkOUT << *mySTD;
	}

	checkOUT.close();
}

#ifdef MPI_PARALLEL
void JobHandler::RunSimulation()
{
	if(myCpuIndex == 0)
	{
		DistributeSimulations();
	}

	ProcessingSimulation();
}

void JobHandler::RunEvaluation()
{
	if(myCpuIndex == 0)
	{
		DistributeEvaluations();
	}

	ProcessingEvaluation();
}

void JobHandler::DistributeSimulations()
{
	// feed node, send one job to each node each time
	int rank = myControl->startRank;
	int chain = 1;
	int alphaIndex = 0;
	double alpha = myControl->alphaSet[0];
	
	// cout << "rank: " << rank << " alpha: " << myControl->alphaSet[0] << endl;
	
	for(int cpuIndex=CPU0; cpuIndex<nProcessors; cpuIndex++) 
	{
		// Send job to each rank
Rprintf("JobHandler::DistributeSimulations(): Initially send cpu\n");
R_FlushConsole();
R_ProcessEvents();
//		cout << "\tInitially send cpu" << cpuIndex << " the job unit of " << "rank:" << rank << ", chain:" << chain << ", alpha:" << alpha << endl;

		MPI_Send(&rank, 1, MPI_INT, cpuIndex, 100, MPI_COMM_WORLD);
		MPI_Send(&chain, 1, MPI_INT, cpuIndex, 200, MPI_COMM_WORLD);
		MPI_Send(&alpha, 1, MPI_DOUBLE, cpuIndex, 300, MPI_COMM_WORLD);

		// update tested rank and chain

		// first alpha set
		alphaIndex ++;
		if(alphaIndex >= myControl->nAlphas)
		{
			// second chain
			alphaIndex = 0;
			chain ++;
			if(chain > myControl->nChains)
			{
				// third rank
				chain = 1;
				rank ++;

				if(rank > myControl->endRank)
				{
					break;
				}
			}
		}

		alpha = myControl->alphaSet[alphaIndex];
	}

	// loop over getting new job requests until there is no more job to be done
//	cout << "nChains " <<  myControl->nChains << " rank " << rank << " alpha " << alpha <<  " chain " << chain << endl;
	while(rank <= myControl->endRank) 
	{
		// Receive results from a note
		int resultIndicator;
		MPI_Status status;
		MPI_Recv(&resultIndicator, 1, MPI_INT, MPI_ANY_SOURCE, 400, MPI_COMM_WORLD, &status);
		int cpuIndex = status.MPI_SOURCE;
		// if(cpuIndex == 0) continue; // this way the first CPU is only used as a distributor
//		cout << "\t\tReceive request from cpu" << cpuIndex << endl;

		// Send the node new job
//		cout << "\tSend cpu" << cpuIndex << " the job unit of " << "rank:" << rank << ", chain:" << chain << ", alpha:" << alpha << endl;
		MPI_Send(&rank, 1, MPI_INT, cpuIndex, 100, MPI_COMM_WORLD);
		MPI_Send(&chain, 1, MPI_INT, cpuIndex, 200, MPI_COMM_WORLD);
		MPI_Send(&alpha, 1, MPI_DOUBLE, cpuIndex, 300, MPI_COMM_WORLD);

		alphaIndex ++;
		if(alphaIndex >= myControl->nAlphas)
		{
			// second chain
			alphaIndex = 0;
			chain ++;
			if(chain > myControl->nChains)
			{
				// third rank
				chain = 1;
				rank ++;
			}
		}
		alpha = myControl->alphaSet[alphaIndex];
	}

	// Tell all the nodes to exit by sending an empty message with DIETAG
	SendEndMessage(100);
}

// processing job on each node
void JobHandler::ProcessingSimulation()
{
	int    rank, chain;
	double alpha;
	bool   FLAG = true;
	MPI_Status status;

	// this cpu wasn't designed to work on jobs
	if(myCpuIndex < CPU0)
	{
		return;
	}

	while( FLAG ) 
	{
		// Receive a message from the master
		MPI_Recv(&rank, 1, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
		if(rank < 0) 
		{
			FLAG = false;
			return;
		} // signal to terminate the process

		MPI_Recv(&chain, 1, MPI_INT, 0, 200, MPI_COMM_WORLD, &status);
		MPI_Recv(&alpha, 1, MPI_DOUBLE, 0, 300, MPI_COMM_WORLD, &status);

		// cout << "\tcpu" << myCpuIndex << " working on job unit of " << "rank:" << rank << " chain:" << chain << " alpha:" << alpha << endl;

		// do simulation
		// Simulation mySimulation(*myData, *myControl);
		// mySimulation.Run(rank, chain, alpha);
		Simulation *mySimulation= new Simulation(*myData, *myControl);
		mySimulation->SetRank( rank );
		mySimulation->SetAlpha( alpha );
		mySimulation->SetChain( chain );
		mySimulation->Run();
		delete mySimulation;

		int jobStatus = 1;
		MPI_Send(&jobStatus, 1, MPI_INT, 0, 400, MPI_COMM_WORLD);
	}
}

// processing job on each node
void JobHandler::ProcessingEvaluation()
{
	int    rank;
	double alpha;
	bool   FLAG = true;
	MPI_Status status;

	// this cpu wasn't designed to work on jobs
	if(myCpuIndex < CPU0)
	{
		return;
	}

	while( FLAG ) 
	{
		// Receive a message from the master
		MPI_Recv(&rank, 1, MPI_INT, 0, 101, MPI_COMM_WORLD, &status);
		if(rank < 0) 
		{
			FLAG = false;
			return;
		} // signal to terminate the process

		MPI_Recv(&alpha, 1, MPI_DOUBLE, 0, 300, MPI_COMM_WORLD, &status);

//		cout << "\tcpu" << myCpuIndex << " working on job unit of " << "rank:" << rank << " alpha:" << alpha << endl;

		// do simulation
		if(myControl->target.find_first_of("ATTERN") == 1)
		{
			HEvaluation he(*myControl, *refClustData);
			he.Run(rank, alpha);
		}
		else
		{
			WEvaluation we(*myControl, *refClustData);
			we.Run(rank, alpha);
		}

		int jobStatus = 1;
		MPI_Send(&jobStatus, 1, MPI_INT, 0, 400, MPI_COMM_WORLD);
	}
}

void JobHandler::Finalizing() const
{
	MPI_Finalize();
//	exit(0);
}

void JobHandler::Barrier() const
{
	MPI_Barrier(MPI_COMM_WORLD);
}

void JobHandler::SendEndMessage(const int signalIndex)
{
	int rank = -1;
	for(int cpuIndex=CPU0; cpuIndex<nProcessors; cpuIndex++)
	{
		MPI_Send(&rank, 1, MPI_INT, cpuIndex, signalIndex, MPI_COMM_WORLD);
	}
}


void JobHandler::ProcessingSimulation(const int rank, const int chain, const double alpha)
{
	if(rank < 0) 
	{
		return;
	}

	// cout << "\tWorking on job unit of " << "rank:" << rank << " chain:" << chain << " alpha:" << alpha << endl;

	// do simulation
	Simulation *mySimulation= new Simulation(*myData, *myControl);
	mySimulation->SetRank( rank );
	mySimulation->SetAlpha( alpha );
	mySimulation->SetChain( chain );
	mySimulation->Run();
	delete mySimulation;
}

void JobHandler::DistributeEvaluations()
{
	// feed node, send one job to each node each time
	int rank = myControl->startRank;
	int alphaIndex = 0;
	double alpha = myControl->alphaSet[0];
	
	// cout << rank << " alpha: " << alpha << endl;
	for(int cpuIndex=CPU0; cpuIndex<nProcessors; cpuIndex++) 
	{
		// Send job to each rank
//		cout << "Send " << cpuIndex << " ==> " << "rank:" << rank << ", alpha:" << alpha << endl;

		MPI_Send(&rank, 1, MPI_INT, cpuIndex, 101, MPI_COMM_WORLD);
		MPI_Send(&alpha, 1, MPI_DOUBLE, cpuIndex, 300, MPI_COMM_WORLD);

		// update tested rank and alpha

		// first alpha set
		alphaIndex ++;
		if(alphaIndex >= myControl->nAlphas)
		{
			// second rank
			alphaIndex = 0;
			rank ++;

			if(rank > myControl->endRank)
			{
				break;
			}
		}

		alpha = myControl->alphaSet[alphaIndex];
	}

	// loop over getting new job requests until there is no more job to be done
	// cout << "alpha: " << alpha << endl;
	while(rank <= myControl->endRank) 
	{
		// Receive results from a note
		int resultIndicator;
		MPI_Status status;
		MPI_Recv(&resultIndicator, 1, MPI_INT, MPI_ANY_SOURCE, 400, MPI_COMM_WORLD, &status);
		int cpuIndex = status.MPI_SOURCE;
		// if(cpuIndex == 0) continue; // this way the first CPU is only used as a distributor
//		cout << "  Receive  " << cpuIndex << " <== " << resultIndicator << endl;

		// Send the node new job
//		cout << "Send " << cpuIndex << " ==> " << "rank:" << rank << ", alpha:" << alpha << endl;
		MPI_Send(&rank, 1, MPI_INT, cpuIndex, 101, MPI_COMM_WORLD);
		MPI_Send(&alpha, 1, MPI_DOUBLE, cpuIndex, 300, MPI_COMM_WORLD);

		alphaIndex ++;
		if(alphaIndex >= myControl->nAlphas)
		{
			// second rank
			alphaIndex = 0;
			rank ++;
		}
		alpha = myControl->alphaSet[alphaIndex];
	}

	// Tell all the nodes to exit by sending an empty message with DIETAG
	SendEndMessage(101);
}


#endif

void JobHandler::Run()
{
#ifdef MPI_PARALLEL
	RunSimulation();
	Barrier();
	RunEvaluation();
	Finalizing();
#else
//printf("JobHandler::Run() #1: entered function...\n"); R_FlushConsole(); R_ProcessEvents();
	for (int rank=myControl->startRank; rank<=myControl->endRank; rank++)
	{
		for (unsigned i=0; i<myControl->alphaSet.size(); i++)
		{
//Rprintf("JobHandler::Run() #1: rank = %d, i = %d\n",rank,i); R_FlushConsole(); R_ProcessEvents();
            myControl->alphaIndex = i;
			Run(rank, myControl->alphaSet[i]);
		}
	}
#endif
}

void JobHandler::Run(int rank, double alpha)
{
	if(myControl->RunType == "evaluation")
	{
		ProcessingEvaluation(rank, alpha);
	}
	else
	{
		// do simulation
//Rprintf("JobHandler::Run() #2: about to instantiate Simulation object...\n"); R_FlushConsole(); R_ProcessEvents();
		Simulation mySimulation(*myData, *myControl, refClustData);
//		Simulation mySimulation(*myData, *myControl, refClustData, H, W);
//Rprintf("JobHandler::Run() #2: about to invoke SetRank()...\n"); R_FlushConsole(); R_ProcessEvents();

		mySimulation.SetRank( rank );
//Rprintf("JobHandler::Run() #2: about to invoke SetAlpha()...\n"); R_FlushConsole(); R_ProcessEvents();
		mySimulation.SetAlpha( alpha );
//Rprintf("JobHandler::Run() #2: about to invoke Run()...\n"); R_FlushConsole(); R_ProcessEvents();
		mySimulation.Run();
//Rprintf("JobHandler::Run() #2: returned from Run()...\n"); R_FlushConsole(); R_ProcessEvents();
	}
}

void JobHandler::ProcessingEvaluation(const int rank, const double alpha)
{
	if(rank < 0) 
	{
		return;
	}

	// do simulation
	if(myControl->target.find_first_of("ATTERN") == 1)
	{
		HEvaluation he(*myControl, *refClustData);
		he.Run(rank, alpha);
	}
	else
	{
		WEvaluation we(*myControl, *refClustData);
		we.Run(rank, alpha);
	}
}
