#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "Parameters.h"
#include "VMatrix.h"
#include "ED.h"
#include "KL.h"
#include "GammaJD.h"
#include "GammaKL.h"
#include "Renyi.h"
#include "DivComb.h"
#include "Pareto.h"
#include "Div2.h"
#include "InverseLink_GammaKL.h"
#include "ClustDataSet.h"
#include "BD.h"
#include "Gamma.h"
#include "NBD.h"
#include "ODP.h"

class Simulation
{
public:
/*
	Simulation(	const VMatrix &data, 
				const ParaControler &control, 
				const ClustDataSet *refClustData,
				const double *matrixH,
				const double *matrixW
				) 
				: 
				myData(data), 
				myControl(control),
				myRefClust( refClustData ),
				H(matrixH),
				W(matrixW)
				{}
*/				
	Simulation(	const VMatrix &data, 
				const ParaControler &control, 
				const ClustDataSet *refClustData
				) 
				: 
				myData(data), 
				myControl(control),
				myRefClust( refClustData )
				{}

	Simulation(	const VMatrix &data, 
				const ParaControler &control 
				) 
				: 
				myData(data), 
				myControl(control),
				myRefClust( NULL )
				{}

	~Simulation() 
	{
	}

	void Run();
	void SetRank(const int rank_) { rank = rank_; }
	void SetAlpha(const double alpha_) { alpha = alpha_; }
	void SetChain(const int chain_) { chain = chain_; }

private:
	const VMatrix &myData;
	const ParaControler &myControl;
	const ClustDataSet *myRefClust;

	int rank, chain;
	double alpha;

//	const double *H;
//	const double *W;

	void Run_ED();
	void Run_KL();
	void Run_Renyi();
	void Run_GammaJD();
	void Run_GammaKL();
	void Run_InverseLink_GammaKL();
	void Run_DivComb();
	void Run_Div2();
	void Run_Pareto();
	void Run_BD();
	void Run_NBD();
	void Run_Gamma();
	void Run_ODP();
};

#endif

