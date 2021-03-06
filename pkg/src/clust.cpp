#include <iostream>
#include <fstream>
#include "muscle.h"
// #include "clustset.h"
#include "clust.h"

// JMM (7/11/2016): Include R header file AFTER system headers. 
// Reference: First few paragraphs in in Section 6 The R API: entry points for C code
// in the online documentation for "Writing R Extensions".
// Also see email from Prof. B. Ripley sent 7/10/2016
// at 3:58 AM entitled "CRAN packages failing to install with modern C++"
// Note: Inclusion of R headers used to be done in the GNMF header files,
// e.g. in
#define R_NO_REMAP
#include "R.h"           // R functions

using namespace std;

Clust::Clust()
	{
	m_Nodes = 0;
	m_uNodeCount = 0;
	m_uLeafCount = 0;
	m_uClusterCount = 0;
	// m_JoinStyle = JOIN_Undefined;
	m_JoinStyle = JOIN_NearestNeighbor; // change undefined to 0
	m_dDist = 0;
	m_uLeafCount = 0;
	m_ptrSet = 0;
	}

Clust::~Clust()
	{
	delete[] m_Nodes;
	delete[] m_dDist;
	// delete[] m_ClusterIndexToNodeIndex;
	}

// root function to generate a hierarchical cluster
void Clust::Create(ClustSet &Set, CLUSTER Method)
{
	m_ptrSet = &Set;

	SetLeafCount(Set.GetLeafCount());

	switch (Method)
	{
	case CLUSTER_UPGMA:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Avg;
		break;

	case CLUSTER_UPGMAMax:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Max;
		break;

	case CLUSTER_UPGMAMin:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Min;
		break;

	case CLUSTER_UPGMB:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Biased;
		break;

	case CLUSTER_NeighborJoining:
		m_JoinStyle = JOIN_NeighborJoining;
		m_CentroidStyle = LINKAGE_NeighborJoining;
		break;

	default:
		Rprintf("Clust::Create: Method noit handled, invalid method\n");
		R_FlushConsole();
		R_ProcessEvents();
        return;	
	}

	if (m_uLeafCount <= 1)
	{
		Rprintf("Clust::Create: m_uLeafCount <= 1, no leaves\n");
		R_FlushConsole();
		R_ProcessEvents();
        return;
	}
	// Quit("Clust::Create: no leaves");

	// how many possible nodes including the leaves
	// m leaves can have m-1 new nodes, plus m leaves, 2m-1
	m_uNodeCount = 2*m_uLeafCount - 1; 

	m_Nodes = new ClustNode[m_uNodeCount];
	// m_ClusterIndexToNodeIndex = new unsigned[m_uLeafCount];

	// initializing node list pointer
	m_ptrClusterList = 0;
	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
	{
		// shouldn't be: ClustNode &Node; define Node fitst and assign mNodes[uNodeIndex] = &Node ?
		ClustNode &Node = m_Nodes[uNodeIndex];
		Node.m_uIndex = uNodeIndex;

		// assign each data points (leaf) to a clust node
		if (uNodeIndex < m_uLeafCount)
		{
			Node.m_uSize = 1;
			Node.m_uLeafIndexes = new unsigned[1];
			Node.m_uLeafIndexes[0] = uNodeIndex;
			AddToClusterList(uNodeIndex);
		}
		// initializing other clust node container for the possible new-merged nodes
		else
			Node.m_uSize = 0;
	}

	// Compute initial distance matrix between leaves
	// SetProgressDesc("Build dist matrix");
	unsigned uPairIndex = 0;

	// total pairwise distances
//	const unsigned uPairCount = (m_uLeafCount*(m_uLeafCount - 1))/2;
	for (unsigned i = 0; i < m_uLeafCount; ++i)
	{
		for (unsigned j = 0; j < i; ++j)
		{
			// this ComputeDist is different than the one in this class (Clust)
			// ClustSet.ComputeDist is more specific for leaves distance, and
			// Clust.ComputeDist is only for nodes distance. The schemes for distance
			// calculation in these two functions are different. Clust.ComputeDist is
			// more general, it only uses the value of leaves distance but doesn't care 
			// how it was calculated. ClustSet.ComputeDist is the one to implement the
			// specific algorithm for data point distance calculation
			const double dDist = m_ptrSet->ComputeDist(*this, i, j);
			SetDist(i, j, dDist);
			/*
			if (0 == uPairIndex%10000)
			Progress(uPairIndex, uPairCount);
			*/
			++uPairIndex;
		}
	}
	// ProgressStepsDone();

	// Call CreateCluster once for each internal node in the tree
	// SetProgressDesc("Build guide tree");
	m_uClusterCount = m_uLeafCount;
//	const unsigned uInternalNodeCount = m_uNodeCount - m_uLeafCount;
	for (unsigned uNodeIndex = m_uLeafCount; uNodeIndex < m_uNodeCount; ++uNodeIndex)
	{
		// unsigned i = uNodeIndex + 1 - m_uLeafCount;
		// Progress(i, uInternalNodeCount);
		// cluster re-orginization is done in CreateCluster/JoinNodes
		CreateCluster();
	}
	// ProgressStepsDone();

}
// forming a new cluster
void Clust::CreateCluster()
	{
	// select the closest pair of nodes (once a pair selected, they form a new node, and
	// the next time selection would only deal with the new node, not the previous pair anymore
	unsigned uLeftNodeIndex;
	unsigned uRightNodeIndex;
	double dLeftLength;
	double dRightLength;
	ChooseJoin(&uLeftNodeIndex, &uRightNodeIndex, &dLeftLength, &dRightLength);

	const unsigned uNewNodeIndex = m_uNodeCount - m_uClusterCount + 1;

	// assign a new node by merging the selected pair and re-orginize cluster
	JoinNodes(uLeftNodeIndex, uRightNodeIndex, dLeftLength, dRightLength,
	  uNewNodeIndex);

	// Compute distances to other clusters
	//
	// the total cluster number is the total leaves number-1, once a new cluster has formed,
	// the leftover cluster number should be decreased by 1, and also this number is used
	// to calculate the node index, so the next node index would increased by 1.
	--m_uClusterCount;
	for (unsigned uNodeIndex = GetFirstCluster(); uNodeIndex != uInsane;
	  uNodeIndex = GetNextCluster(uNodeIndex))
		{
		// these two have just merged to a new node
		if (uNodeIndex == uLeftNodeIndex || uNodeIndex == uRightNodeIndex)
			continue;

		// only distances of the new node to other nodes are need to be re-calculated
		// so new node itself is skipped too
		if (uNewNodeIndex == uNodeIndex)
			continue;

		// again calculate distance of new node to all other nodes
		const double dDist = ComputeDist(uNewNodeIndex, uNodeIndex);
		SetDist(uNewNodeIndex, uNodeIndex, dDist);
		}

	// this is just like the above section, but is for metric instead of distance
	for (unsigned uNodeIndex = GetFirstCluster(); uNodeIndex != uInsane;
	  uNodeIndex = GetNextCluster(uNodeIndex))
		{
		if (uNodeIndex == uLeftNodeIndex || uNodeIndex == uRightNodeIndex)
			continue;

		if (uNewNodeIndex == uNodeIndex)
			continue;

#if	REDLACK
		const double dMetric = ComputeMetric(uNewNodeIndex, uNodeIndex);
		InsertMetric(uNewNodeIndex, uNodeIndex, dMetric);
#endif
		}
	}

// Choose which pair of nodes need to be merged next
void Clust::ChooseJoin(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
  double *ptrdLeftLength, double *ptrdRightLength)
	{
	switch (m_JoinStyle)
		{
	case JOIN_NearestNeighbor:
		ChooseJoinNearestNeighbor(ptruLeftIndex, ptruRightIndex, ptrdLeftLength,
		  ptrdRightLength);
		return;
	case JOIN_NeighborJoining:
		ChooseJoinNeighborJoining(ptruLeftIndex, ptruRightIndex, ptrdLeftLength,
		  ptrdRightLength);
		return;
		}
	Rprintf("Clust::ChooseJoin: m_JoinStyle not handled, Invalid join style\n");
	R_FlushConsole();
	R_ProcessEvents();
    return;
	}

// Choose node pair to be merged in next step by Nearest Neighbor scheme
// the difference between NearestNeeoghbora and Neighbor is that NearestNeighbor
// doesn't care the difference of either the parters with other nodes, it only uses
// the distance between these two specific nodes and the heights of the two
void Clust::ChooseJoinNearestNeighbor(unsigned *ptruLeftIndex,
  unsigned *ptruRightIndex, double *ptrdLeftLength, double *ptrdRightLength)
	{
//	const unsigned uClusterCount = GetClusterCount();

	unsigned uMinLeftNodeIndex;
	unsigned uMinRightNodeIndex;
	GetMinMetric(&uMinLeftNodeIndex, &uMinRightNodeIndex);

	double dMinDist = GetDist(uMinLeftNodeIndex, uMinRightNodeIndex);

	const double dLeftHeight = GetHeight(uMinLeftNodeIndex);
	const double dRightHeight = GetHeight(uMinRightNodeIndex);

	*ptruLeftIndex = uMinLeftNodeIndex;
	*ptruRightIndex = uMinRightNodeIndex;
	*ptrdLeftLength = dMinDist/2 - dLeftHeight;
	*ptrdRightLength = dMinDist/2 - dRightHeight;
	}

// Choose node pair to be merged in next step by Neighbor scheme
void Clust::ChooseJoinNeighborJoining(unsigned *ptruLeftIndex,
  unsigned *ptruRightIndex, double *ptrdLeftLength, double *ptrdRightLength)
	{
//	const unsigned uClusterCount = GetClusterCount();
	unsigned uMinLeftNodeIndex;
	unsigned uMinRightNodeIndex;
	GetMinMetric(&uMinLeftNodeIndex, &uMinRightNodeIndex);

	const double dDistLR = GetDist(uMinLeftNodeIndex, uMinRightNodeIndex);
	const double rL = Calc_r(uMinLeftNodeIndex);
	const double rR = Calc_r(uMinRightNodeIndex);

	const double dLeftLength = (dDistLR + rL - rR)/2;
	const double dRightLength = (dDistLR - rL + rR)/2;

	*ptruLeftIndex = uMinLeftNodeIndex;
	*ptruRightIndex = uMinRightNodeIndex;
	*ptrdLeftLength = dLeftLength;
	*ptrdRightLength = dRightLength;
	}

// merging two nodes into a new cluster and re-orginize cluster
void Clust::JoinNodes(unsigned uLeftIndex, unsigned uRightIndex, double dLeftLength,
  double dRightLength, unsigned uNodeIndex)
	{
	ClustNode &Parent = m_Nodes[uNodeIndex];
	ClustNode &Left = m_Nodes[uLeftIndex];
	ClustNode &Right = m_Nodes[uRightIndex];

	Left.m_dLength = dLeftLength;
	Right.m_dLength = dRightLength;

	Parent.m_ptrLeft = &Left;
	Parent.m_ptrRight = &Right;

	Left.m_ptrParent = &Parent;
	Right.m_ptrParent = &Parent;

	const unsigned uLeftSize = Left.m_uSize;
	const unsigned uRightSize = Right.m_uSize;
	const unsigned uParentSize = uLeftSize + uRightSize;
	Parent.m_uSize = uParentSize;

	// assert(0 == Parent.m_uLeafIndexes);
	Parent.m_uLeafIndexes = new unsigned[uParentSize];

	// memcpy is more efficient than hard copy
	const unsigned uLeftBytes = uLeftSize*sizeof(unsigned);
	const unsigned uRightBytes = uRightSize*sizeof(unsigned);
	memcpy(Parent.m_uLeafIndexes, Left.m_uLeafIndexes, uLeftBytes);
	memcpy(Parent.m_uLeafIndexes + uLeftSize, Right.m_uLeafIndexes, uRightBytes);

	// delete children node, re-orginzing cluster
	DeleteFromClusterList(uLeftIndex);
	DeleteFromClusterList(uRightIndex);
	AddToClusterList(uNodeIndex);
	}

// Calc_r calculate the average distance of one node to all other nodes
double Clust::Calc_r(unsigned uNodeIndex) const
	{
	const unsigned uClusterCount = GetClusterCount();
	
	// when nClusterCount = 2, only one cluster is there, i.e. uNodeIndex refer to the root node, 
	if (2 == uClusterCount)
		return 0;

	double dSum = 0;
	for (unsigned i = GetFirstCluster(); i != uInsane; i = GetNextCluster(i))
		{
		if (i == uNodeIndex)
			continue;
		dSum += GetDist(uNodeIndex, i);
		}
	return dSum/(uClusterCount - 2);
	}

// function to calculate the distance between a pair of nodes
double Clust::ComputeDist(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	switch (m_CentroidStyle)
		{
	case LINKAGE_Avg:
		return ComputeDistAverageLinkage(uNewNodeIndex, uNodeIndex);

	case LINKAGE_Min:
		return ComputeDistMinLinkage(uNewNodeIndex, uNodeIndex);

	case LINKAGE_Max:
		return ComputeDistMaxLinkage(uNewNodeIndex, uNodeIndex);

	case LINKAGE_Biased:
		return ComputeDistMAFFT(uNewNodeIndex, uNodeIndex);

	case LINKAGE_NeighborJoining:
		return ComputeDistNeighborJoining(uNewNodeIndex, uNodeIndex);
		}
	Rprintf("Clust::ComputeDist: m_CentroidStyle not handled, invalid centroid style\n");
	R_FlushConsole();
	R_ProcessEvents();
    return 0.0;
	}

double Clust::ComputeDistMinLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const double dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const double dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	return (dDistL < dDistR ? dDistL : dDistR);
	}

double Clust::ComputeDistMaxLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const double dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const double dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	return (dDistL > dDistR ? dDistL : dDistR);
	}

double Clust::ComputeDistAverageLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const double dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const double dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	// changed Sep. 28, 2007 to reflect the complete average
	// return (dDistL + dDistR)/2;
	int lSize = m_Nodes[uLeftNodeIndex].m_uSize;
	int rSize = m_Nodes[uRightNodeIndex].m_uSize;
	return (dDistL*lSize + dDistR*rSize)/(lSize+rSize);
	}

double Clust::ComputeDistNeighborJoining(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const double dDistLR = GetDist(uLeftNodeIndex, uRightNodeIndex);
	const double dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const double dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	const double dDist = (dDistL + dDistR - dDistLR)/2;
	return dDist;
	}

// This is a mysterious variant of UPGMA reverse-engineered from MAFFT source.
double Clust::ComputeDistMAFFT(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);

//	const double dDistLR = GetDist(uLeftNodeIndex, uRightNodeIndex);
	const double dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const double dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	const double dMinDistLR = (dDistL < dDistR ? dDistL : dDistR);
	const double dSumDistLR = dDistL + dDistR;
	const double dDist = dMinDistLR*(1 - g_dSUEFF) + dSumDistLR*g_dSUEFF/2;
	return dDist;
	}

unsigned Clust::GetClusterCount() const
	{
	return m_uClusterCount;
	}

const ClustNode &Clust::GetNode(unsigned uNodeIndex) const
	{
	if (uNodeIndex >= m_uNodeCount)
	{
		Rprintf("ClustNode::GetNode()?\n");
		R_FlushConsole();
		R_ProcessEvents();
        return m_Nodes[0];
	}
	return m_Nodes[uNodeIndex];
	}

bool Clust::IsLeaf(unsigned uNodeIndex) const
	{
	return uNodeIndex < m_uLeafCount;
	}

unsigned Clust::GetClusterSize(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	return Node.m_uSize;
	}

unsigned Clust::GetLeftIndex(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	if (0 == Node.m_ptrLeft)
	{
		Rprintf("Clust::GetLeftIndex: 0 == Node.m_ptrLeft, no leaf\n");
		R_FlushConsole();
		R_ProcessEvents();
        return 0;
	}
	return Node.m_ptrLeft->m_uIndex;
	}

unsigned Clust::GetRightIndex(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	if (0 == Node.m_ptrRight)
	{
		Rprintf("Clust::GetRightIndex: 0 == Node.m_ptrRight, no leaf\n");
		R_FlushConsole();
		R_ProcessEvents();
        return 0;
	}
	return Node.m_ptrRight->m_uIndex;
	}

double Clust::GetLength(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	return Node.m_dLength;
	}

void Clust::SetLeafCount(unsigned uLeafCount)
	{
	if (uLeafCount <= 1)
	{
		Rprintf("Clust::SetLeafCount(): uLeafCount <= 1, non-positive leaf count\n");
		R_FlushConsole();
		R_ProcessEvents();
        return;
	}

	m_uLeafCount = uLeafCount;
	const unsigned uNodeCount = GetNodeCount();

// Triangular matrix size excluding diagonal (all zeros in our case).
	m_uTriangularMatrixSize = (uNodeCount*(uNodeCount - 1))/2;
	m_dDist = new double[m_uTriangularMatrixSize];
	}

unsigned Clust::GetLeafCount() const
	{
	return m_uLeafCount;
	}

// Get vector index 
unsigned Clust::VectorIndex(unsigned uIndex1, unsigned uIndex2) const
	{
	const unsigned uNodeCount = GetNodeCount();
	if (uIndex1 >= uNodeCount || uIndex2 >= uNodeCount)
	{
		Rprintf("Clust::VectorIndex: uIndex1 >= uNodeCount or uIndex2 >= uNodeCount\n");
		R_FlushConsole();
		R_ProcessEvents();
        return 0;
	}

	unsigned v;
	if (uIndex1 >= uIndex2)
		v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
	else
		v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
	// assert(v < m_uTriangularMatrixSize);
	return v;
	}

// Get distance from vector 
double Clust::GetDist(unsigned uIndex1, unsigned uIndex2) const
	{
	unsigned v = VectorIndex(uIndex1, uIndex2);
	return m_dDist[v];
	}

// Set distance into vector index
void Clust::SetDist(unsigned uIndex1, unsigned uIndex2, double dDist)
	{
	unsigned v = VectorIndex(uIndex1, uIndex2);
	m_dDist[v] = dDist;
	}

double Clust::GetHeight(unsigned uNodeIndex) const
	{
	if (IsLeaf(uNodeIndex))
		return 0;

	const unsigned uLeftIndex = GetLeftIndex(uNodeIndex);
	const unsigned uRightIndex = GetRightIndex(uNodeIndex);
	const double dLeftLength = GetLength(uLeftIndex);
	const double dRightLength = GetLength(uRightIndex);
	const double dLeftHeight = dLeftLength + GetHeight(uLeftIndex);
	const double dRightHeight = dRightLength + GetHeight(uRightIndex);
	return (dLeftHeight + dRightHeight)/2;
	}

const char *Clust::GetNodeName(unsigned uNodeIndex) const
	{
	if (!IsLeaf(uNodeIndex))
	{
		Rprintf("Clust::GetNodeName: uNodeIndex is not leaf\n");
		R_FlushConsole();
		R_ProcessEvents();
        return m_ptrSet->GetLeafName(uNodeIndex);
	}
	return m_ptrSet->GetLeafName(uNodeIndex);
	}

unsigned Clust::GetNodeId(unsigned uNodeIndex) const
	{
	if (uNodeIndex >= GetLeafCount())
		return 0;
	return m_ptrSet->GetLeafId(uNodeIndex);
	}

unsigned Clust::GetLeaf(unsigned uNodeIndex, unsigned uLeafIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	const unsigned uLeafCount = Node.m_uSize;
	if (uLeafIndex >= uLeafCount)
	{
		Rprintf("Clust::GetLeaf: uLeafIndex >= uLeafCount, invalid index\n");
		R_FlushConsole();
		R_ProcessEvents();
        return 0;
	}

	const unsigned uIndex = Node.m_uLeafIndexes[uLeafIndex];
	if (uIndex >= m_uNodeCount)
	{
		Rprintf("Clust::GetLeaf: uIndex >= m_uNodeCount, index out of range\n");
		R_FlushConsole();
		R_ProcessEvents();
        return 0;
	}
	return uIndex;
	}

unsigned Clust::GetFirstCluster() const
	{
	if (0 == m_ptrClusterList)
		return uInsane;
	return m_ptrClusterList->m_uIndex;
	}

unsigned Clust::GetNextCluster(unsigned uIndex) const
	{
	ClustNode *ptrNode = &m_Nodes[uIndex];
	if (0 == ptrNode->m_ptrNextCluster)
		return uInsane;
	return ptrNode->m_ptrNextCluster->m_uIndex;
	}

void Clust::DeleteFromClusterList(unsigned uNodeIndex)
	{
	// assert(uNodeIndex < m_uNodeCount);
	ClustNode *ptrNode = &m_Nodes[uNodeIndex];
	ClustNode *ptrPrev = ptrNode->m_ptrPrevCluster;
	ClustNode *ptrNext = ptrNode->m_ptrNextCluster;

	if (0 != ptrNext)
		ptrNext->m_ptrPrevCluster = ptrPrev;
	if (0 == ptrPrev)
		{
		// assert(m_ptrClusterList == ptrNode);
		m_ptrClusterList = ptrNext;
		}
	else
		ptrPrev->m_ptrNextCluster = ptrNext;

	ptrNode->m_ptrNextCluster = 0;
	ptrNode->m_ptrPrevCluster = 0;
	}

void Clust::AddToClusterList(unsigned uNodeIndex)
	{
	// assert(uNodeIndex < m_uNodeCount);
	ClustNode *ptrNode = &m_Nodes[uNodeIndex];

	// insert node before the last node on the list
	if (0 != m_ptrClusterList)
		m_ptrClusterList->m_ptrPrevCluster = ptrNode;

	ptrNode->m_ptrNextCluster = m_ptrClusterList;
	ptrNode->m_ptrPrevCluster = 0;

	m_ptrClusterList = ptrNode;
	}

double Clust::ComputeMetric(unsigned uIndex1, unsigned uIndex2) const
	{
	switch (m_JoinStyle)
		{
	case JOIN_NearestNeighbor:
		return ComputeMetricNearestNeighbor(uIndex1, uIndex2);

	case JOIN_NeighborJoining:
		return ComputeMetricNeighborJoining(uIndex1, uIndex2);
		}
	Rprintf("Clust::ComputeMetric: m_JoinStyle not handled\n");
	R_FlushConsole();
	R_ProcessEvents();
	return 0;
	}

double Clust::ComputeMetricNeighborJoining(unsigned i, unsigned j) const
	{
	double ri = Calc_r(i);
	double rj = Calc_r(j);
	double dij = GetDist(i, j);
	double dMetric = dij - (ri + rj);
	return dMetric;
	}

double Clust::ComputeMetricNearestNeighbor(unsigned i, unsigned j) const
	{
	return GetDist(i, j);
	}

double Clust::GetMinMetricBruteForce(unsigned *ptruIndex1, unsigned *ptruIndex2) const
	{
	unsigned uMinLeftNodeIndex = uInsane;
	unsigned uMinRightNodeIndex = uInsane;
	double dMinMetric = PLUS_INFINITY;
	for (unsigned uLeftNodeIndex = GetFirstCluster(); uLeftNodeIndex != uInsane;
	  uLeftNodeIndex = GetNextCluster(uLeftNodeIndex))
		{
		for (unsigned uRightNodeIndex = GetNextCluster(uLeftNodeIndex);
		  uRightNodeIndex != uInsane;
		  uRightNodeIndex = GetNextCluster(uRightNodeIndex))
			{
			double dMetric = ComputeMetric(uLeftNodeIndex, uRightNodeIndex);
			if (dMetric < dMinMetric)
				{
				dMinMetric = dMetric;
				uMinLeftNodeIndex = uLeftNodeIndex;
				uMinRightNodeIndex = uRightNodeIndex;
				}
			}
		}
	*ptruIndex1 = uMinLeftNodeIndex;
	*ptruIndex2 = uMinRightNodeIndex;
	return dMinMetric;
	}

double Clust::GetMinMetric(unsigned *ptruIndex1, unsigned *ptruIndex2) const
	{
	return GetMinMetricBruteForce(ptruIndex1, ptruIndex2);
	}

