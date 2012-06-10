#include "Evaluation.h"
#include "R.h"           // R functions

using namespace std;

// this class is designed to be used after the simulations, so all factorized matrixes with the
// same rank and alpha would be read in, and the evaluation will be done in the sense of how
// consistent the multiple runs performance for the specific set of parameters.
// 

void Evaluation::Run()
{
	for (int rank=myControl.startRank; rank<=myControl.endRank; rank++)
	{
		for (unsigned alphaIndex=0; alphaIndex<myControl.alphaSet.size(); alphaIndex++)
		{
			Run(rank, myControl.alphaSet[alphaIndex]);

		}
	}
}

void Evaluation::Run(const int _rank, const double _alpha)
{
	rank = _rank;
	alpha = _alpha;
/*
	ostringstream outMatrix;
	outMatrix << myControl.rootName << "." << myControl.scheme << ".k" << rank;
	if(myControl.scheme == "Renyi") outMatrix << "_alpha" << alpha;
	factorMatrixRootName = outMatrix.str();
*/
//	ostringstream outLog;
//	outLog << myControl.rootName << "." << myControl.scheme << "_" << myControl.clustScheme << ".k" << rank;
//	if(myControl.scheme == "Renyi") outLog << "_alpha" << alpha;
//	factorMatrixLogRootName = outLog.str();

	//cout << "Evaluation::Run(): SetParameters... rank " << rank << " alpha " << alpha << endl;
 	SetParameters();

	// read in individul matrix
	//cout << "Evaluation::Run(): ProcessingFiles" << endl;
Rprintf("Evaluation::Run: Processing files...\n"); R_FlushConsole(); R_ProcessEvents();
	ProcessingFiles();

	// do hierarchical clustering
	//cout << "Evaluation::Run(): HierarchCustering" << endl;
	HierarchCustering();

	// do evaluations based on selected cluster
	//cout << "Evaluation::Run(): DoClusterBasedTest" << endl;
	DoClusterBasedTest();

	// MC rate test
	//cout << "Evaluation::Run(): MCRateTest" << endl;
	MCRateTest();

	//cout << "Evaluation::Run(): Closing" << endl;
	Closing();

}

void Evaluation::Closing()
{
	ReorderConsensusMatrix();
/*
	outLOG << "rank: " << rank;
	outLOG << " alpha: " << alpha;
	outLOG << " cophenetic: " << cophenetic;
	outLOG << " spearman: " << spearman;
	outLOG << " adjustedRandIndex: " << adjustedRandIndex;
	outLOG << " MutualInformationIndex: " << nmi;
	outLOG << " kupaOverall: " << kupaOverall;
	outLOG << " tPointScatter: " << tPointScatter;
	outLOG << " wPointScatter: " << wPointScatter;
	outLOG << " nu: " << nu;
	outLOG << " MC: " << MC;
	outLOG << " indnu: " << this->indnu;
	outLOG << " indNMI: " << this->indNMI;
	outLOG << " indRandIndex: " << this->indRandIndex;
	outLOG << " Kupas => ";
	for (unsigned i=0; i<kupa.size(); i++)	
	{
		outLOG << " clt" << i << " " << kupa[i];
	}
	outLOG << endl;

	outLOG.close();
*/
	/*
	// start children element for evaluation parameters
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "cophenetic", "%8.6f", cophenetic);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "spearman", "%8.6f", spearman);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "kupaOverall", "%8.6f", kupaOverall);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "tPointScatter", "%8.6f", tPointScatter);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "wPointScatter", "%8.6f", wPointScatter);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "adjustedRandIndex", "%8.6f", adjustedRandIndex);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "MutualInformationIndex", "%8.6f", nmi);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "nu", "%8.6f", nu);
	xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "MC", "%8.6f", MC);

	xmlTextWriterStartElement(xmlOutLog, BAD_CAST "kupas");
	for (unsigned i=0; i<kupa.size(); i++)
	{
		xmlTextWriterWriteFormatElement(xmlOutLog, BAD_CAST "kupa", "%8.6f", kupa[i]);
	}
	xmlTextWriterEndElement(xmlOutLog);

	xmlTextWriterEndElement(xmlOutLog);

	xmlTextWriterEndDocument(xmlOutLog);
	xmlFreeTextWriter(xmlOutLog);
	*/
}

void Evaluation::CopySaveProperties( SSMap &features, int run)
{
/*
	outLOG << "run" << run;
	for (SSMap::const_iterator ip=features.begin(); ip!=features.end(); ip++)
	{
		outLOG << " " << ip->first << " " << ip->second; 
	}
	outLOG << endl;
*/
	/*
	xmlTextWriterStartElement(xmlOutLog, BAD_CAST "IndividualInfo");
	xmlTextWriterWriteFormatAttribute(xmlOutLog, BAD_CAST "run", "%d", run);
	for (SSMap::const_iterator ip=features.begin(); ip!=features.end(); ip++)
	{
		xmlTextWriterWriteElement(xmlOutLog, BAD_CAST ip->first.c_str(), BAD_CAST ip->second.c_str());
	}
	xmlTextWriterEndElement(xmlOutLog);
	*/
}

void Evaluation::SaveClusterInfo() const
{
/*     
	ostringstream clustFile;
	clustFile << factorMatrixLogRootName << fileType << "clt";
	ofstream outCLT(clustFile.str().c_str(), ios::out);
	outCLT << ccClustData;
	outCLT.close();
	*/
}

void Evaluation::DoConsensusMatrix()
{
	indNMI = 0.0;
	indRandIndex = 0.0;
	indnu = 0.0;
	// cout << "Starting DoConsensusMatrix..." << endl;
	ccMatrix->rowNames = targetMatrixes[0]->columnNames;
	ccMatrix->columnNames = ccMatrix->rowNames;

	for (int i=0; i<myControl.nChains; i++)
	{
		// cout << "run " << i << endl;
		if(myControl.clustScheme == "Binary")
		{
			targetMatrixes[i]->GetConnectivity();
			ccMatrix->Adding( targetMatrixes[i]->connectivity->data );
			delete targetMatrixes[i]->connectivity;
		}
		else
		{
			targetMatrixes[i]->GetConnectivityReal( myControl.clustScheme );

			if(myControl.clustScheme.find("HC", 0) != myControl.clustScheme.npos)
			{
				IndHC( targetMatrixes[i]->realConnet );
				if(&refClustData != 0)
				{
					indnu += MisMatchedRate( ccClustData );
					indNMI += nmi;
					indRandIndex += adjustedRandIndex;
				}
				// build connectivity matrix and calculate consensus matrix
				for (int i=0; i<ccClustData.nClusters; i++)
				{
					for (unsigned j=0; j<ccClustData.members[i].size(); j++)
					{
						for (unsigned k=0; k<ccClustData.members[i].size(); k++)
						{
							ccMatrix->data[ ccClustData.members[i][j]->index ] [ ccClustData.members[i][k]->index ] ++;
						}
					}
				}

			}
			else
			{
				ccMatrix->Adding( targetMatrixes[i]->realConnet->data );
			}

			delete targetMatrixes[i]->realConnet;
		}

	}

	indnu /= myControl.nChains;
	indNMI /= myControl.nChains;
	indRandIndex /= myControl.nChains;
	ccMatrix->Scaling( 1.0 / myControl.nChains );

}

void Evaluation::IndHC( CMatrix *indC )
{
        ClustSetNMF nmfSet(*indC);

        // build hirarchical cluster
        myClust.Create(nmfSet, CLUSTER_UPGMA);

        // cluster consensus test for k clusters, k should be rank
        SelectClusters(rank);
}

void Evaluation::DoClusterBasedTest()
{
	// Calculate cluster consensus for each given cluster, it's defined as the mean consensus index between
	// all pairs of samples belonging to the same cluster.
	ClustConsensusTest();
	/*
	cout << "\n>>Cluster consensus (kupas):\n";
	for (unsigned i=0; i<kupa.size(); i++)
	{
		cout << "\tclustter" << i+1 << ": " << kupa[i];
	}
	cout << endl;
	*/

	// the overall cluster consensus is defined as the mean consensus index between all pairs if samples
	// belonging to the same cluster and averaged across the rank clusters
	OverallConsensusTest();
	// cout << "\n>>Overall cluster consensus: " << kupaOverall << endl;
	// the total point scatter measures the degree to diversity of all samples

	tPointScatter = OverallPointScatter();
	// cout << "\n>>Total point scatter: " << tPointScatter << endl;
	// overall within-cluster point scatter measures the degree to which samples assigned to the same cluster
	// are close to one another, averaged across the rank clusters.

	wPointScatter = ClustPointScatter();
	// cout << "\n>>Overall within-cluster point scatter: " << wPointScatter << endl;
}


void Evaluation::BDClustering()
{
	GetSigmas();
	SetClustDataSetSigma( bdClustData );
	// bdClustData.CompleteClustSet( ccMatrix->rowNames );
	// outLOG << "\n>>Clustering result for BD-style, based on sigma:\n";
	// outLOG << bdClustData;
}


void Evaluation::MCRateTest()
{
	// do MCRate evaluation if applied
	// cout << "mcFlag " <<  myControl.mcFlag << endl;
	if(! myControl.mcFlag)
	{
		adjustedRandIndex = -1;
		nmi = -1;
		nu = -1.0;
		MC = -1.0;
		return;
	}

	// Now only test for the right rank
	// cout << refClustData;
	// cout << "refClustData.nClusters " << refClustData.nClusters << endl;
	MC = MisMatchedRate( ccClustData );
	if(refClustData.nClusters == ccClustData.nClusters)
	{
		nu = MCRate();
	}
	else
	{
		nu = -1.0;
	}

	if(myControl.stdShift > 0.0)
	{
		//do BD-style clustering
		BDClustering();
		MC = MisMatchedRate( bdClustData );
		// cout << "MC " << MC << endl;
	}

	return;
}

void Evaluation::ReorderConsensusMatrix()
{
	/*
	ostringstream ccFile0;
	ccFile0 << factorMatrixLogRootName << fileType << "cc0";
	ofstream outCC0(ccFile0.str().c_str(), ios::out);
	outCC0 << *ccMatrix;
	outCC0.close();
	*/

	// re-order the consensus matrix for analysis and plot purpose
	// make sure the index (name) is also reordered correspondingly
	ClustNode &Parent = myClust.m_Nodes[myClust.GetFirstCluster()];
	GetOrderedIndex(Parent);
	// output data points in order in the hierarchical cluster structure
/*
	ostringstream ordFile;
	ordFile << factorMatrixLogRootName << fileType << "ord";
	ofstream outORD(ordFile.str().c_str(), ios::out);

	outORD << "#Ordered sequence of data points for rank(K)=" << rank << " alpha=" << alpha << ":" << endl;
	for (unsigned i=0; i<orderedSequence.size(); i++)
	{
		outORD << ccMatrix->columnNames[orderedSequence[i]] << endl;
	}

	outORD.close();
	ordFile.clear();
*/
	// re-order rowName
	Utility::ReOrderVector<string>(ccMatrix->rowNames, orderedSequence);

	// re-order columnName
	Utility::ReOrderVector<string>(ccMatrix->columnNames, orderedSequence);

	// re-order data; first reorder columns
	for (unsigned i=0; i<ccMatrix->data.size(); i++)
	{
		Utility::ReOrderVector<double>(ccMatrix->data[i], orderedSequence);
	}

	// second reorder rows
	ccMatrix->Transposing();
	for (unsigned i=0; i<ccMatrix->data.size(); i++)
	{
		Utility::ReOrderVector<double>(ccMatrix->data[i], orderedSequence);
	}
/*
	// output re-ordered consensusMatrix
	ostringstream ccFile;
	ccFile << factorMatrixLogRootName << fileType << "cc";
	ofstream outCC(ccFile.str().c_str(), ios::out);
	outCC << *ccMatrix;
	outCC.close();

	// out put a data file for gnuplot
	ostringstream plotFile;
	plotFile << factorMatrixLogRootName << fileType << "plot";

	ofstream outPlot(plotFile.str().c_str(), ios::out);
	for (unsigned i=0; i<ccMatrix->data.size(); i++)
	{
		for (unsigned j=0; j<ccMatrix->data[i].size(); j++)
		{
			outPlot << i+1 << "\t" << j+1 << "\t" << ccMatrix->data[i][j] << endl;
		}
		outPlot << "\n";
	}
	outPlot.close();
	plotFile.clear();
*/
}

void Evaluation::GetOrderedIndex(ClustNode &Parent)
{
	ClustNode &Left = myClust.m_Nodes[Parent.m_ptrLeft->m_uIndex];
	ClustNode &Right = myClust.m_Nodes[Parent.m_ptrRight->m_uIndex];

	if( myClust.IsLeaf(Left.m_uIndex) )
	{
		orderedSequence.push_back(Left.m_uLeafIndexes[0]);
	}
	else
	{
		GetOrderedIndex(Left);
	}

	if( myClust.IsLeaf(Right.m_uIndex) )
	{
		orderedSequence.push_back(Right.m_uLeafIndexes[0]);
	}
	else
	{
		GetOrderedIndex(Right);
	}
}


void Evaluation::HierarchCustering()
{
	// most evaluation parameters should be come from consensus matrix
	DoConsensusMatrix();

	// define dataset to be clustered
	// ClustSetNMF nmfSet(ccMatrix->data);
	ClustSetNMF nmfSet(*ccMatrix);

	// build hirarchical cluster
	myClust.Create(nmfSet, CLUSTER_UPGMA);

	// calculate cophenetic coefficient factor
	GetCopheneticFactor();

	// cluster consensus test for k clusters, k should be rank
	SelectClusters(rank);
	SaveClusterInfo();

	// ReorderConsensusMatrix();
}


// cophenetic correlation for cluster tree is defined as the linear correlation coefficient
// between the cophenetic distances obtained from the tree, and the original distances 
// (or dissimilarities) used to construct the tree. Thus, it is a measure of how faithfully 
// the tree represents the dissimilarities among observations
//
// The cophenetic distance between two observations is represented in a dendrogram by the
// height of the link at which those two observations are first joined. That height is the 
// distance between the two subclusters that are merged by that link
void Evaluation::GetCopheneticFactor()
{
	// get original distances used to construct the tree
	const unsigned nPairCount = myClust.m_uLeafCount * (myClust.m_uLeafCount - 1) / 2;

	// convert clust double array dist to double vector dist
	array1d origDist(nPairCount);
	for (unsigned i=0; i<nPairCount; i++)
	{
		origDist[i] = myClust.m_dDist[i];
	}

	// get cophenetic distances from the tree, from top clust to end, get leaf index and calculate vector index
	array1d cophDist(nPairCount);
	for (unsigned uNodeIndex = myClust.GetRootNodeIndex(); uNodeIndex > myClust.GetLeafCount()-1; 
		uNodeIndex--)
	{
		ClustNode &Parent = myClust.m_Nodes[uNodeIndex];
		ClustNode &Left = myClust.m_Nodes[Parent.m_ptrLeft->m_uIndex];
		ClustNode &Right = myClust.m_Nodes[Parent.m_ptrRight->m_uIndex];

		for (unsigned i=0; i<Left.m_uSize; i++)
		{
			for (unsigned j=0; j<Right.m_uSize; j++)
			{
				cophDist[ myClust.VectorIndex(Left.m_uLeafIndexes[i], Right.m_uLeafIndexes[j]) ] =
					myClust.GetDist(Left.m_uIndex, Right.m_uIndex);
			}
		}
	}

	// using pearson correlation scheme to return cophenetic factor
	cophenetic = Utility::PearsonCorrelation(origDist, cophDist);
	// outLOG << ">>Cophenetic correlation coefficient: " << cophenetic << endl;

	vector<int> temp( origDist.size() );
    for (unsigned i=0; i < origDist.size(); i++) temp[i] = i;
	sort(temp.begin(), temp.end(), SortVectorIndex<double>(origDist));

	vector<int> ranks1( origDist.size() );
	for (unsigned i=0; i < temp.size(); i++) ranks1[ temp[i] ] = i;

    for (unsigned i=0; i < cophDist.size(); i++) temp[i] = i;
	sort(temp.begin(), temp.end(), SortVectorIndex<double>(cophDist));

	vector<int> ranks2( cophDist.size() );
	for (unsigned i=0; i<temp.size(); i++) ranks2[ temp[i] ] = i;
	spearman = 0.0;
	for (unsigned i=0; i<temp.size(); i++)
	{
		double diff = ranks1[i] - ranks2[i];
		spearman += diff * diff;
	}

	spearman *= 6.0;
	double ooo = ranks1.size() * (ranks1.size() * ranks1.size() - 1);
	spearman = 1 - spearman / ooo;
	// outLOG << ">>Spearman's rank coefficient: " << spearman << endl;
}


void Evaluation::ClustConsensusTest()
{
	// cout << "before test" << selectedNodes.size() << endl;
	kupa.resize( selectedNodes.size() );
	clustSizeInPairs.resize( selectedNodes.size() );
	
	for (unsigned i=0; i<selectedNodes.size(); i++)
	{
		ClustNode *node = selectedNodes[i];
		// cout << "i = " << i <<  " " <<  node->m_uSize << endl;

		if(node->m_uSize == 1)
		{
			kupa[i] = 0.0;
			continue;
		}

		kupa[i] = 0.0;
		for (unsigned j=0; j<node->m_uSize; j++)
		{
			for (unsigned k=j+1; k<node->m_uSize; k++)
			{
				kupa[i] += ccMatrix->data[node->m_uLeafIndexes[j]][node->m_uLeafIndexes[k]];
			}
		}

		clustSizeInPairs[i] = static_cast<int>(node->m_uSize * (node->m_uSize - 1) / 2);

		kupa[i] /= clustSizeInPairs[i];
	}
}

void Evaluation::OverallConsensusTest()
{
	// first need ClusterConsensus
	kupaOverall = 0.0;
	int denominator = 0; 
	for (unsigned i=0; i<kupa.size(); i++)
	{
		kupaOverall += clustSizeInPairs[i] * kupa[i];
		denominator += clustSizeInPairs[i];
	}

	kupaOverall /= denominator;
}

double Evaluation::OverallPointScatter()
{
	// this value is independent of cluster assignment
	double pointScatter = 0.0;
	for (unsigned i=0; i<ccMatrix->data.size(); i++)
	{
		for (unsigned j=i+1; j<ccMatrix->data.size(); j++)
		{
			pointScatter += 1.0 - ccMatrix->data[i][j];
		}
	}

	return(pointScatter);
}

void Evaluation::SelectClusters(const int nClustNodes)
{	
	// take the fact that the upper cluster nodes must hold the more diversed groups, all
	// nodes with smaller index than endNodeIndex are considered cluster entity to be evaluated
	unsigned endNodeIndex = myClust.GetRootNodeIndex() - nClustNodes + 2;

	// selectedNodes = new ClustNode[nClustNodes];
	selectedNodes.resize(nClustNodes);

	int selected = 0;
	for ( unsigned uNodeIndex=myClust.GetRootNodeIndex(); uNodeIndex>=endNodeIndex;
		uNodeIndex-- )
	{
		ClustNode &Parent = myClust.m_Nodes[uNodeIndex];
		unsigned leftIndex = Parent.m_ptrLeft->m_uIndex;
		ClustNode &Left = myClust.m_Nodes[leftIndex];
		unsigned rightIndex = Parent.m_ptrRight->m_uIndex;
		ClustNode &Right = myClust.m_Nodes[rightIndex];

		if( Left.m_uIndex < endNodeIndex )
		{
			selectedNodes[selected] = &Left;
			selected ++;
		}

		if( Right.m_uIndex < endNodeIndex )
		{
			selectedNodes[selected] = &Right;
			selected ++;
		}
	}

	if(selected != nClustNodes)
	{
Rprintf("Evaluation::SelectClusters: selected nodes != designed number\n");
R_FlushConsole();
R_ProcessEvents();
//		cout << "Evaluation::SelectClusters: selected nodes " << selected << 
//			" != designed number " << nClustNodes << "!" << endl;
//		JOEexit(-1);
        return;
	}

	// extract clust data set 
	ccClustData.nClusters = nClustNodes;
	ccClustData.members.resize( nClustNodes );
	struct member sp;
	for (unsigned i=0; i<selectedNodes.size(); i++)
	{
		ClustNode *node = selectedNodes[i];
		ccClustData.members[i].resize( node->m_uSize );

		// first sort indexes
		for (unsigned j=0; j<node->m_uSize; j++)
		{
			ccClustData.members[i][j] = new member(node->m_uLeafIndexes[j], ccMatrix->rowNames[node->m_uLeafIndexes[j]]);
		}

		sort(ccClustData.members[i].begin(), ccClustData.members[i].end(), sp);
	}
}

double Evaluation::ClustPointScatter()
{
	double pointScatter = 0.0;
	int denominator = 0; 
	for (unsigned i=0; i<selectedNodes.size(); i++)
	{
		// ClustNode node = myClust.m_Nodes[selectedNodes[i]->m_uIndex];
		ClustNode *node = selectedNodes[i];
		for (unsigned j=0; j<node->m_uSize; j++)
		{
			for (unsigned k=j+1; k<node->m_uSize; k++)
			{
				pointScatter += 1.0 - ccMatrix->data[node->m_uLeafIndexes[j]][node->m_uLeafIndexes[k]];
			}
		}
		denominator += clustSizeInPairs[i];
	}

	pointScatter /= (denominator * tPointScatter);

	return(pointScatter);
}

double Evaluation::ClustPointScatter(const ClustDataSet &targetClust)
{
	// ofstream temp("temp", ios::out);
	double pointScatter = 0.0;
	double fullValue = 0;
	for (int i=0; i<targetClust.nClusters; i++)
	{
		if(targetClust.members[i].size() <= 1) continue;

		fullValue += targetClust.members[i].size() * (targetClust.members[i].size() - 1) / 2;
		for (unsigned j=0; j<targetClust.members[i].size(); j++)
		{
			for (unsigned k=j+1; k<targetClust.members[i].size(); k++)
			{
				pointScatter += 1.0 - ccMatrix->data[targetClust.members[i][j]->index][targetClust.members[i][k]->index];
			}
		}
	}

	if(fullValue == 0.0)
	{
		return -1.0;
	}
	else
	{
		return(pointScatter/fullValue);
	}
}

// !!! must first do mapping between repeats
void Evaluation::SetClustDataSetSigma(ClustDataSet &toBeSet)
{
	toBeSet.nClusters = rank;
	toBeSet.members.resize( toBeSet.nClusters );
	for (unsigned i=0; i<average.size(); i++)
	{
		for (unsigned j=0; j<average[i].size(); j++)
		{
			double stdev = average[i][j] / sigma[i][j];
			if( stdev > myControl.stdShift)
			{
				toBeSet.members[i].push_back( new member(j, "") );
			}
		}
	}
}

// calculate the NMI
void Evaluation::MutualInformationIndex( ClustDataSet &testClustData, array1d &enrichment )
{
	double rowSumComb = 0;
	for (unsigned i=0; i<(unsigned)refClustData.nClusters; i++)
	{
		double nSize = static_cast<double> (refClustData.members[i].size());
		rowSumComb += nSize * log ( nSize / ccMatrix->nRows );
	}

	double columnSumComb = 0;
	for (unsigned i=0; i<(unsigned)testClustData.nClusters; i++)
	{
		double nSize = static_cast<double> (testClustData.members[i].size());
		columnSumComb += nSize * log ( nSize / ccMatrix->nRows );
	}

	double confusionComb = 0;
	for (int i=0; i<refClustData.nClusters; i++)
	{
		for (int j=0; j<testClustData.nClusters; j++)
		{
			double common = enrichment[i * testClustData.nClusters + j];
			if(common > 0)
			{
				confusionComb += common * log ( common * ccMatrix->nRows / (refClustData.members[i].size() * testClustData.members[j].size()) );
			}
		}
	}

	nmi = confusionComb / sqrt( rowSumComb * columnSumComb );
}


// calculate the adjusted rand index
void Evaluation::RandIndex( ClustDataSet &testClustData, array1d &enrichment )
{
	double rowSumComb = 0;
	for (unsigned i=0; i<(unsigned)refClustData.nClusters; i++)
	{
		rowSumComb += (double)Utility::Combination42( refClustData.members[i].size() );
	}

	double columnSumComb = 0;
	for (unsigned i=0; i<(unsigned)testClustData.nClusters; i++)
	{
		columnSumComb += (double)Utility::Combination42( testClustData.members[i].size() );
	}

	double confusionComb = 0;
	for (int i=0; i<refClustData.nClusters; i++)
	{
		for (int j=0; j<testClustData.nClusters; j++)
		{
			confusionComb += (double) Utility::Combination42( static_cast<int>( enrichment[i * testClustData.nClusters + j] ) );
		}
	}

	double totalComb = (double) Utility::Combination42( ccMatrix->nRows );
	double rowColumnSum = columnSumComb + rowSumComb;
	double rowColumnMul = rowSumComb * columnSumComb / totalComb;
	// cout << "totalComb " << totalComb << " rowColumnSum " << rowColumnSum << " rowColumnMul " << rowColumnMul << endl;
	// cout << "confusionComb " <<  confusionComb << endl;

	adjustedRandIndex = (confusionComb - rowColumnMul) / (0.5*rowColumnSum - rowColumnMul);

	// outLOG << "The adjusted Rand Index: " << adjustedRandIndex << endl;
}

double Evaluation::MisMatchedRate(ClustDataSet &testClustData)
{
	// get enrichment
	array1d enrichment;
	// cout << "Build enrichment" << endl;
	for (int i=0; i<refClustData.nClusters; i++)
	{
		for (int j=0; j<testClustData.nClusters; j++)
		{
			double enriched = 0.0;
			int chkPoint = 0;
			for (int m=0; m<(int)refClustData.members[i].size(); m++)
			{
				for (int n=chkPoint; n<(int)testClustData.members[j].size(); n++)
				{
					if(refClustData.members[i][m]->index == testClustData.members[j][n]->index)
					{
						enriched ++;
						chkPoint = n + 1;
						break;
					}
					else if(refClustData.members[i][m]->index < testClustData.members[j][n]->index)
					{
						chkPoint = n;
						break;
					}
				}
			}
			enrichment.push_back(enriched);
		}
	}

	// cout << "Getting map" << endl;
	IIMap rtMap = GetMapping( enrichment, testClustData.nClusters );
	// cout << "Map is done" << endl;

	// calculate misclassification (MC) rate and print out the cluster mapping
	// int totalPoints = 0;
	double matches = 0;
	for (int i=0; i<refClustData.nClusters; i++)
	{
		IIMap::const_iterator it = rtMap.find(i);
		if(it == rtMap.end()) continue;
		double commons = enrichment[i * testClustData.nClusters + rtMap[i] ];
 		// int countSize = refClustData.members[i].size() + testClustData.members[ rtMap[i] ].size() - static_cast<int>(commons);
		// mismatches += countSize - commons;
		matches += commons;
		// totalPoints += countSize;
/*
		outLOG << "   refClustData " << i << "(" << refClustData.members[ i ].size() << ")"
			<<  " ==> testCluster " << rtMap[i] << "(" << testClustData.members[rtMap[i]].size() << ")"
			<<  " matched " << commons << endl;
*/
	}

	// cout << "invoke RandIndex" << endl;
	RandIndex( testClustData, enrichment );
	MutualInformationIndex( testClustData, enrichment );
	// cout << "mismatches " << mismatches << " totalPoints " << totalPoints << endl;

	return (1.0 - matches / ccMatrix->nRows);
}


void Evaluation::OrderMatrixes()
{
	// correlation coefficients
	for (unsigned i=1; i<targetMatrixes.size(); i++)
	{
		array1d scores4Ranking;
		for (int j=0; j<rank; j++)
		{
			for (int k=0; k<rank; k++)
			{
				double coor = Utility::PearsonCorrelation( targetMatrixes[i]->data[j], targetMatrixes[0]->data[k] ); 
				scores4Ranking.push_back( coor );
			}
		}
		
		maps4chains.push_back( GetMapping( scores4Ranking, rank ) );
	}
}

IIMap Evaluation::GetMapping( array1d scores4Ranking, const int dim )
{
	vector<int> origIndex( scores4Ranking.size() );
	for (unsigned i=0; i<origIndex.size(); i++) origIndex[i] = i;
	sort(origIndex.begin(), origIndex.end(), SortVectorIndex<double>(scores4Ranking));

	IIMap map0;
	IIMap map1;
	for (unsigned i=0; i<origIndex.size(); i++)
	{
		int index1 = (int) (origIndex[i] / dim);
		int index2 = origIndex[i] - index1 * dim;

		if(map0.find( index1 ) == map0.end() && map1.find( index2 ) == map1.end() )
		{
			map0[index1] = index2;
			map1[index2] = index1;
		}
	}

	return map0;
}

double Evaluation::MCRate()
{
	return ClustPointScatter( refClustData );
}


