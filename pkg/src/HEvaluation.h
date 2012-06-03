#pragma once
#include "Evaluation.h"
#include "HMatrix.h"

class HEvaluation :
	public Evaluation
{
public:
/*
	HEvaluation(
		const ParaControler &_control,
		const ClustDataSet &_refClust,
		const double *matrixH,
		const double *matrixW
		)
		:
		Evaluation( _control, _refClust ),
		H(matrixH),
		W(matrixW)
		{}
*/
	HEvaluation(
		const ParaControler &_control,
		const ClustDataSet &_refClust
		)
		:
		Evaluation( _control, _refClust ) 
		{}

	~HEvaluation(void) {};

private:
//	const double *H;
//	const double *W;

	void   ProcessingFiles();
	void   CombineMatrixes();
	void   ReadinAssembledFiles();
	void   SetParameters();
	void   GetSigmas();
};
