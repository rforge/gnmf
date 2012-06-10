#include "Simulation.h"
#include "R.h"           // R functions

using namespace std;

void Simulation::Run()
{
//Rprintf("Simulation::Run(): entered function...\n"); R_FlushConsole(); R_ProcessEvents();
	if(myControl.scheme == "ED")
	{
		Run_ED();
	}
	else if(myControl.scheme == "KL")
	{
//Rprintf("Simulation::Run(): Detected KL...\n"); R_FlushConsole(); R_ProcessEvents();
		Run_KL();
	}
	else if(myControl.scheme == "Renyi")
	{
//Rprintf("Simulation::Run(): detected Renyi...\n"); R_FlushConsole(); R_ProcessEvents();
		Run_Renyi();
	}
	else if(myControl.scheme == "GammaJD")
	{
		Run_GammaJD();
	}
	else if(myControl.scheme == "GammaKL")
	{
		Run_GammaKL();
	}
	else if(myControl.scheme == "Gamma")
	{
		Run_Gamma();
	}
	else if(myControl.scheme == "DivComb")
	{
		Run_DivComb();
	}
	else if(myControl.scheme == "Div2")
	{
		Run_Div2();
	}
	else if(myControl.scheme == "InverseLink_GammaKL")
	{
		Run_InverseLink_GammaKL();
	}
	else if(myControl.scheme == "Pareto")
	{
		Run_Pareto();
	}
	else if(myControl.scheme == "BD")
	{
		Run_BD();
	}
	else if(myControl.scheme == "NBD")
	{
		Run_NBD();
	}
	else if(myControl.scheme == "ODP")
	{
		Run_ODP();
	}
	else
	{
Rprintf("Simulation::Run(): error detected...\n");
R_FlushConsole();
R_ProcessEvents();
//		std::cerr << "Not a valid scheme: " << myControl.scheme << "!" << std::endl;
	}
}

void Simulation::Run_ED()
{
	ED myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_KL()
{
//Rprintf("Simulation::Run_KL(): entered function...\n"); R_FlushConsole(); R_ProcessEvents();
	KL myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
//Rprintf("Simulation::Run_KL(): about to invoke myUpdate.Run()...\n"); R_FlushConsole(); R_ProcessEvents();
	myUpdate.Run( *myRefClust );
}


void Simulation::Run_BD()
{
	BD myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}


void Simulation::Run_NBD()
{
	NBD myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_Pareto()
{
	Pareto myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_InverseLink_GammaKL()
{
	InverseLink_GammaKL myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_Gamma()
{
	Gamma myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_Renyi()
{
//Rprintf("Simulation::Run_Renyi(): entered function...\n"); R_FlushConsole(); R_ProcessEvents();
	// cout << "Run_Renyi" << endl;
	Renyi myUpdate(myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
//Rprintf("Simulation::Run_Renyi(): about to invoke myUpdate.Run()...\n"); R_FlushConsole(); R_ProcessEvents();
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_GammaJD()
{
	GammaJD myUpdate(myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_GammaKL()
{
	GammaKL myUpdate(myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_DivComb()
{
	DivComb myUpdate(myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_Div2()
{
	Div2 myUpdate(myData, myControl );
	myUpdate.rank = rank;
	// myUpdate.alpha = alpha;
	myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}

void Simulation::Run_ODP()
{
	ODP myUpdate( myData, myControl );
	myUpdate.rank = rank;
	myUpdate.alpha = alpha;
	// myUpdate.chain = chain;
	if(! myControl.oneStepFlag) myUpdate.chain = chain;
	myUpdate.Run( *myRefClust );
}
