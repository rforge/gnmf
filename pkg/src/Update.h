#ifndef UPDATE_H_
#define UPDATE_H_

#include "HMatrix.h"
#include "WMatrix.h"
#include "VMatrix.h"
#include "WeightingMatrix.h"
#include "Parameters.h"
#include "HEvaluation.h"
#include "WEvaluation.h"
// #include "Utility.h"
#include "ParaControler.h"
#include "ClustDataSet.h"
#include <fstream>
#include <sstream>
#include <vector>

class Update
{
public:
	int    rank, chain;
	double alpha, alphaR;
	double average, sd;
	const VMatrix &myData;
	const ParaControler &myControl;
	HMatrix *pattern;
	WMatrix *amplitude;
	VMatrix *mock;
	WeightingMatrix *weightMatrix, *updateFactorA, *updateFactorP;

	Update(	
		const VMatrix &_data, 
		const ParaControler &_control
		) :	
		myData(_data), 
		myControl(_control) 
	{
		thetaMatrix = 0;
	}

	~Update()
	{
		if(thetaMatrix != 0) delete thetaMatrix;
	}

	void Run(const int rank, const int chain, const double _alpha);
	void Run( const ClustDataSet &_refClustData );

protected:
	double errors;
	array2d prevMock;
	double consecutiveError;
	bool   converged;
	int    totalSteps;
	WMatrix *thetaMatrix; // could be HMatrix too, doesn't matter

	void   GetInitial();
	void   Normalizing();

	void   Closing();
	void   UpdateFeatures();
	void   SaveMatrixes() const;
	void GenerateMock();
	void GenerateMock(const WMatrix &w, const HMatrix &h);
	void CalculateWeights();
	virtual void   UpdateAmplitudeMatrix() = 0;
	virtual void   UpdatePatternMatrix() = 0;
	virtual int    CheckConvergency( int round = 0 );
	virtual void   SaveIntermediaFiles() = 0;
	void   IndMatrixEvaluation();
	void   HSimulation( const ClustDataSet &_refClustData );
	void   WSimulation( const ClustDataSet &_refClustData );
	std::string wMatrixName, hMatrixName;
private:
	template <class T> std::string toString( T inputValue )
	{
		std::ostringstream ss;
		ss << inputValue;

		return ss.str();
	}
	void BuildThetaMatrix();
};

#endif

