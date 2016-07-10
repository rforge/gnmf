#ifndef _JOBHANDLER_H_
#define _JOBHANDLER_H_

#ifndef WIN32
#include <unistd.h>
#else
#include <direct.h> // Windows lib for getcwd
#endif

#include <fstream>
#include "Parameters.h"
#include "DataMatrix.h"
#include "VMatrix.h"
#include "ParaControler.h"
#include "Simulation.h"
#include "HEvaluation.h"
#include "WEvaluation.h"
#include "ClustDataSet.h"

// JMM (7/10/2016): Do not include R include files in GNMF header files, because sometimes
// GNMF header files include other GNMF header files. This makes it difficult to ensure
// that R header files are included after GNMF header files. It is probably best to just
// include R header files directly in the CPP files where they are used.
//#define R_NO_REMAP
//#include "R.h"           // R functions

#ifdef MPI_PARALLEL
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
// #include </usr/local/mpich2/include/mpi.h>
#include </apps/mpich2/include/mpi.h>
#endif

class JobHandler
{
public:
    
	JobHandler(int argc, char **argv) : 
		nArguments(argc), 
		argumentList(argv)
	{ 
#ifdef MPI_PARALLEL
		MPI_Init(&nArguments, &argumentList);
		MPI_Comm_rank(MPI_COMM_WORLD, &myCpuIndex);
		MPI_Comm_size(MPI_COMM_WORLD, &nProcessors); 
#else
		nProcessors = 1;
		myCpuIndex = 0;
#endif
	}

    // New constructor that gets input data from R rather than from
    // files on disk and from the command line.
	JobHandler(double *matrixX,
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
               double *Hsparseness);

	~JobHandler() 
	{
		delete myData;
		delete myControl;
		if( mySTD ) delete mySTD;
		if( refClustData ) delete refClustData;
	}

	void Run();
	void RunSimulation();
	void RunEvaluation();
	void VerifyArgumentList();
	void GetInput();
	void GetInputText();
	void PrintInput() const;
	void Finalizing() const;
	void Barrier() const;

private:
	int nArguments;
	char **argumentList;
	std::string mode;
	std::string batchType;

	int jobStatus; // default is fine, if something goes wrong, it will be changed to -1
	int nProcessors;
	int myCpuIndex;
	int thisrank;
	double thisalpha;
	// std::string rootName;
	ParaControler *myControl;
	VMatrix *myData, *mySTD;
	ClustDataSet *refClustData;

	void Run(int, double);
	void PrintUsage();
	void DistributeSimulations();
	void ProcessingSimulation();
	void ProcessingSimulation(const int rank, const int chain, const double alpha);
	void DistributeEvaluations();
	void ProcessingEvaluation();
	void ProcessingEvaluation(const int rank, const double alpha);
	void SendEndMessage(const int signalIndex);
	
//	double *H;
//	double *W;
};

#endif
