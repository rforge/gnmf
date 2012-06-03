#pragma once
#include "Evaluation.h"
#include "WMatrix.h"

class WEvaluation :
	public Evaluation
{
public:
/*
	WEvaluation(
		const ParaControler &_control,
		const ClustDataSet &_refClust,
		const double *matrixH,
		const double *matrixW
		)
		: Evaluation( _control, _refClust ),
		H(matrixH),
		W(matrixW)
		{}
*/
	WEvaluation(
		const ParaControler &_control,
		const ClustDataSet &_refClust
		)
		: Evaluation( _control, _refClust ) 
		{}

	~WEvaluation(void) {}

private:
//	const double *H;
//	const double *W;

	void   ProcessingFiles();
	void   CombineMatrixes();
	void   ReadinAssembledFiles();
	void   SetParameters();
	void   GetSigmas();
};
