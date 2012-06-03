#include "ParaControler.h"

using namespace std;

ostream& operator << (std::ostream &out, const ParaControler &myControl)
{
	out << "working dir: " << myControl.workingDir << endl; 
	out << "source data file: " << myControl.sourceFile << endl;
	out << "uncertainty file: " << myControl.uncertFile << endl;
	out << "reference cluster file: " << myControl.refClusterFile << endl;
	out << "Total genes (observations): " << std::setw(10) << myControl.nObservations << endl;
	out << "Total samples (parameters): " << std::setw(10) << myControl.nParameters << endl;
	out << "Simulation will try rank from " << myControl.startRank << " to " << myControl.endRank << endl;
	out << "Each rank will run " << myControl.nChains << " times with different initialized matrixes" << endl;
	out << "Each simulation will run up to " << myControl.nUpdateSteps << " rounds" << endl << endl;
	out << "Simulation will test " << myControl.nAlphas << " alphas: ";
	for (unsigned i=0; i<myControl.alphaSet.size(); i++)
	{
		out << myControl.alphaSet[i] << " ";
	}
	out << endl;
	out << "Step normalizing " << myControl.Normalize << endl;
	out << "Initial Scaling " << myControl.InitialScaling << endl;;
	out << "Measure scheme: " << myControl.scheme << "  Target matrix: " << myControl.target << std::endl;
	out << "rootName: " << myControl.rootName << std::endl;
	out << "base seed for initializing matrixes: " << myControl.baseSeed << std::endl;
	// out << "Some control parameter setting:\n";
	// out << "    UPPERBOUND   " << myControl.UPPERBOUND << "\n";
	// out << "    ACCEPTEDDIFF " << std::scientific << myControl.ACCEPTEDDIFF << "\n";
	// out << "    EPS          " << std::scientific << myControl.eps << "\n";
	// out << "    weight (N/A) " << myControl.weight << std::endl;

	out << "Parameters for evaluation setting:\n";
	out << "PreEvaluation flag: " << myControl.preEvaluation << std::endl;
	out << "Zscore threshold: " << myControl.zscoreThreshold << "  ";
	out	<< "Proportion cutoff: " << myControl.proportionCut << "  ";
	out	<< "Sigma threhold: " << myControl.stdShift << std::endl;
	out << "nu test? " << myControl.mcFlag << std::endl;

	return (out);
}

// read in parameters for controling the simulations
istream& operator >> (std::istream &in, ParaControler &myControl)
{
	string termName;
	in >> termName >> myControl.workingDir;
	in >> termName >> myControl.sourceFile;
	in >> termName >> myControl.uncertFile;
	in >> termName >> myControl.refClusterFile;
	in >> termName >> myControl.rootName;
	in >> termName >> myControl.nObservations;
	in >> termName >> myControl.nParameters;
	in >> termName >> myControl.nChains;
	in >> termName >> myControl.startRank;
	in >> termName >> myControl.endRank;
	in >> termName >> myControl.nUpdateSteps;
	in >> termName >> myControl.baseSeed;
	in >> termName >> myControl.scheme;
	in >> termName >> myControl.target;
	in >> termName >> myControl.nAlphas;
	myControl.alphaSet.resize(myControl.nAlphas);
	for (int i=0; i<myControl.nAlphas; i++)
	{
		in >> myControl.alphaSet[i];
	}
	in >> termName >> myControl.preEvaluation;
	in >> termName >> myControl.zscoreThreshold;
	in >> termName >> myControl.proportionCut;
	in >> termName >> myControl.stdShift;
	in >> termName >> myControl.mcFlag;
	in >> termName >> myControl.Normalize;
	in >> termName >> myControl.InitialScaling;
	in >> termName >> myControl.clustScheme;
	in >> termName >> myControl.theta;
	in >> termName >> myControl.ConvergeTest;

	return (in);
}


void ParaControler::GetParameters( string xmlFile )
{
	/*
	// cout << xmlFile.c_str() << endl;
	xmlDocPtr doc = xmlReadFile(xmlFile.c_str(), NULL, 0);
	if(doc == NULL)
	{
		std::cerr << "Failed to parse " << xmlFile << endl;
		return;
	}

	xmlNodePtr root_node = xmlDocGetRootElement(doc);
	// cout << "Here starting!" << endl;
	for (xmlNodePtr cur_node = root_node->children; cur_node != NULL; cur_node=cur_node->next)
	{
		if(cur_node->type != XML_ELEMENT_NODE) continue;
		// cout << cur_node->name << endl;

		if( xmlStrEqual(cur_node->name, BAD_CAST "alphas") == 1 )
		{
			// cout << "Here for alphas!" << endl;
			for (xmlNodePtr alpha_node = cur_node->children; alpha_node != NULL; alpha_node = alpha_node->next)
			{
				if(alpha_node->type != XML_ELEMENT_NODE) continue;
				// cout << (char *)alpha_node->name << " " << (char*) alpha_node->children->content << endl;
				alphaSet.push_back( atof( (char *) alpha_node->children->content ) );
			}
		}
		else
		{
			// may need to cast: reinterpret_case<const char *>
			if(cur_node->children != NULL) parameters[ string( (char *) cur_node->name) ] = string( (char *) cur_node->children->content );
		}
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	rootName = parameters["rootName"];
	workingDir = parameters["workingDir"];
	sourceFile = parameters["sourceFile"];
	convergeTestStep = 20;
	refClusterFile = parameters["refClustFile"];
	uncertFile = parameters["uncertFile"];
	scheme = parameters["scheme"];
	target = parameters["target"];
	nObservations = atoi( parameters["nRows"].c_str() );
	nParameters = atoi( parameters["nColumns"].c_str() );
	nChains = atoi( parameters["nChains"].c_str() );
	startRank = atoi( parameters["startRank"].c_str() );
	endRank = atoi( parameters["endRank"].c_str() );
	nUpdateSteps = atoi( parameters["nSteps"].c_str() );
	nAlphas = atoi( parameters["nAlphas"].c_str() );
	baseSeed = atoi( parameters["baseSeed"].c_str() );
	preEvaluation = atoi( parameters["preEvaluation"].c_str() );
	mcFlag = atoi( parameters["mcrate_flag"].c_str() );
	zscoreThreshold = atof( parameters["zscore"].c_str() );
	proportionCut = atof( parameters["proportion"].c_str() );
	stdShift = atof( parameters["bd"].c_str() );
	Normalize = atoi( parameters["Normalize"].c_str() );
	InitialScaling = atoi( parameters["InitialScaling"].c_str() );
	idealization = atof( parameters["idealization"].c_str() );
	theta = atof( parameters["theta"].c_str() );
	clustScheme = parameters["clustScheme"];
	ConvergeTest = atoi(parameters["ConvergeTest"].c_str());

  	// cout << "Finished reading and casting" << endl;
	*/
}
