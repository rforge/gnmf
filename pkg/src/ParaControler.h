#ifndef PARACONTROLER_H_
#define PARACONTROLER_H_

#include "Parameters.h"
#include <iomanip>
#include <fstream>
#include <map>
// #include <libxml/parser.h>
// #include <libxml/tree.h>
/*
#ifndef WIN32
#else
#include <objbase.h>
#include <msxml6.h>
#endif
*/

class ParaControler
{
	friend std::istream &operator >> (std::istream &in, ParaControler &myControler);
	friend std::ostream &operator << (std::ostream &out, const ParaControler &myControler);

public:
	std::string rootName;
	std::string workingDir;
	std::string sourceFile;
	std::string refClusterFile;
	std::string uncertFile;
	std::string scheme;
	std::string target;
	std::string clustScheme;
	std::string RunType;
	double idealization;
	SSMap parameters;
	int    nObservations, nParameters, nChains, startRank, endRank;
	int    nUpdateSteps, nAlphas;
	long   baseSeed;
	bool   preEvaluation, mcFlag;
	bool   oneStepFlag;
	bool   Normalize;
	bool   InitialScaling;
	bool   ConvergeTest;
	bool   HDFile;
	double zscoreThreshold, proportionCut, stdShift;
	double theta;
	int convergeTestStep;
	array1d alphaSet;

	// Pointers to memory for two output R matrices.
	double *H;
	double *W;

	// Pointers to memory for output diagnostic parameters (to return to R).
	double *converged;
	double *totalSteps;
	double *error;
	double *consecutiveError;
	double *Wsparseness;
	double *Hsparseness;
	
    int    alphaIndex; // Index indicating which alpha is currently being processed.

	ParaControler()
	{
		// /*
		// assign defaults 
		this->startRank = 2;
		this->endRank = 2;
		this->scheme = scheme_;
		this->target = target_;
		this->Normalize = Normalization_;
		this->InitialScaling = InitialScaling_;
		this->nAlphas = nAlphas_;
		this->nChains = nChains_;
		this->nUpdateSteps = nSteps_;
		this->ConvergeTest = ConvergeTest_;
		this->clustScheme = clustScheme_;
		this->RunType = RunType_;
		this->mcFlag = false;
		this->preEvaluation = false;
		this->nAlphas = 1;
		this->alphaSet.push_back(1.0);
		this->baseSeed = SEED;
		this->theta = THETA;
		this->convergeTestStep = 20;
		this->HDFile = false;
		this->uncertFile = "";
		this->refClusterFile = "";
		this->proportionCut = -1.0;
		this->stdShift = -1.0;
		this->idealization = 1.0;
		
		this->H                = 0;
		this->W                = 0;
		this->converged        = 0;
		this->totalSteps       = 0;
		this->error            = 0;
		this->consecutiveError = 0;
		this->Wsparseness      = 0;
		this->Hsparseness      = 0;
	}
	// */

	~ParaControler() {}

	// necessary parameter setting
	void setObservatonNum(const int &_nObser) { nObservations = _nObser; }
	void setParameterNum(const int &_nPara) { nParameters = _nPara; }
	void setChainNum(const int &_nChains) { nChains = _nChains; }
	void setUpdateSteps(const int &_steps) { nUpdateSteps = _steps; }
	void setConvergeTestStep(const int &_n) { convergeTestStep = _n; }
	void setRankRange(const int &_start, const int &_end) 
	{ 
		startRank = _start;
		endRank = _end;
	}
	void GetParameters( std::string xmlFile );

	// optional setting, these have default values
	/*
	void setEPS(const double &_eps) { EPS = _eps; }
	void setWeight(const double &_weight) { weight = _weight; }
	void setUpperBound(const double &_ub) { upperBound = _ub; }
	void setAcceptedDiff(const double &_af) { ACCEPTEDDIFF = _af; }
	*/
	// void setMethod(const int &_methodFlag) { simulationScheme = methodFlag; }

	// parameters for evaluation
	void setPreEvaluationFlag(bool &_preEvaluation) { preEvaluation = _preEvaluation; }
	void setMethod(std::string &_method) { scheme = _method; }
	void setZCutoff(double &_cutoff) { zscoreThreshold = _cutoff; } // zscore cutoff 2.0
	void setPCutoff(double &_ratio) { proportionCut = _ratio; }    // proportion cutoff 0.2
	void setSigmaCutoff(double &_std) { stdShift = _std; }         // sigma cutoff 3.0
};

#endif

