#ifndef EVALUATION_H_
#define EVALUATION_H_

#include "R.h"           // R functions

#include "Parameters.h"
#include "muscle.h"
#include "clust.h"
#include "ClustSetNMF.h"
#include "Utility.h"
#include "DataMatrix.h"
#include "CMatrix.h"
#include "WMatrix.h"
#include "HMatrix.h"
#include "ParaControler.h"
#include "ClustDataSet.h"
#include <limits>
#include <fstream>
#include <sstream>
// #include <libxml/encoding.h>
// #include <libxml/xmlwriter.h>

#ifndef WIN32
#include <unistd.h> // check whether a file is existed or not, getcwd
#include <sys/stat.h> // getting information about a file
#include <sys/types.h>
#else
#include <stdlib.h> 
#include <io.h> // required by _access function in Windows
#include <direct.h>
// #include <objbase.h>
// #include <msxml6.h>
#endif

class Evaluation
{
public:
	std::vector< FactorMatrix *> targetMatrixes; // since W matrix will be transposed
	std::vector< FactorMatrix * > patterns;
	std::vector< FactorMatrix * > amplitudes;
	// the output files from simulation are named as rootName_rank_chain_alpha.A (P)
	Evaluation( 
		const ParaControler &_control,
		const ClustDataSet &_refClust
		) :
        myControl(_control),
		refClustData(_refClust)
	{
		ccMatrix = 0;
	}

	~Evaluation() 
	{
		if(ccMatrix != 0) delete ccMatrix;
	}

	void   Run();
	void   Run(const int rank, const double alpha);
	void   Run(const int rank, const double alpha, const double *matrixH, const double *matrixW);
	
protected:
	int  nRows, nColumns;
	int  rank;
	double alpha;
	const ParaControler &myControl;
	std::vector< IIMap > maps4chains;
	array2d average;
	array2d sigma;
	std::string factorMatrixRootName;
	std::string factorMatrixLogRootName;

	std::string fileType;
	CMatrix *ccMatrix;
	// xmlTextWriterPtr xmlOutLog;

	std::ofstream outLOG;

	virtual void   ProcessingFiles() = 0;
	virtual void   CombineMatrixes() = 0;
	virtual void   ReadinAssembledFiles() = 0;
	virtual void   SetParameters() = 0;
	virtual void   GetSigmas() = 0;
	void   CopySaveProperties( SSMap &, int);
	void   OrderMatrixes();
	IIMap GetMapping( 
		array1d scores4Ranking, 
		const int dim
		);

private:
	double kupaOverall;
	double cophenetic;
	double spearman;
	double nu;
	double MC;
	double adjustedRandIndex;
	double nmi;
	double tPointScatter;
	double wPointScatter; // overall within clusters
	double indNMI;
	double indRandIndex;
	double indnu;
	array1d kupa;
	Clust myClust;
	ClustDataSet ccClustData, bdClustData;
	const ClustDataSet &refClustData;
	std::vector<ClustNode*> selectedNodes;
	std::vector<int> clustSizeInPairs;
	std::vector<int> orderedSequence;

	void   HierarchCustering();
	void   SaveClusterInfo() const;
	void   DoClusterBasedTest();
	void   Closing();
	double MisMatchedRate(ClustDataSet &target); 
	double MCRate();
	void   IndHC( CMatrix *c);
	void   DoConsensusMatrix();
	void   SelectClusters(const int nClustNodes);
	void   ClustConsensusTest();
	void   OverallConsensusTest();
	void   SetClustDataSetSigma(ClustDataSet &toBeSet);
	double ClustPointScatter();
	double OverallPointScatter();
	double ClustPointScatter(const ClustDataSet &targetClust);
	// void   GetPatternShapes();
	void   MCRateTest();
	void   ReorderConsensusMatrix();
	void   BDClustering();
	void   GetCopheneticFactor();
	void   RandIndex(ClustDataSet &testClustData, array1d &enrichment);
	void   MutualInformationIndex(ClustDataSet &testClustData, array1d &enrichment);
	void   GetOrderedIndex(ClustNode &Parent);

};

#endif

