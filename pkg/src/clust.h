#ifndef Clust_h
#define Clust_h

#include "muscle.h"
#include <map>
#include <vector>
#include "ClustDataSet.h"

class Clust;
class ClustNode;
class ClustSet;

enum CLUSTER
{
	CLUSTER_UPGMA,
	CLUSTER_UPGMAMax,
	CLUSTER_UPGMAMin,
	CLUSTER_UPGMB,
	CLUSTER_NeighborJoining
};

enum JOIN
{
	JOIN_NearestNeighbor,
	JOIN_NeighborJoining
};

enum LINKAGE
{
	LINKAGE_Min,
	LINKAGE_Avg,
	LINKAGE_Max,
	LINKAGE_NeighborJoining,
	LINKAGE_Biased
};

class ClustNode
	{
public:
	ClustNode()
		{
		m_uIndex = uInsane;
		m_uSize = uInsane;
		m_dLength = (double) dInsane;
		m_ptrLeft = 0;
		m_ptrRight = 0;
		m_ptrParent = 0;
		m_ptrNextCluster = 0;
		m_ptrPrevCluster = 0;
		m_uLeafIndexes = 0;
		}
	~ClustNode()
		{
		delete[] m_uLeafIndexes;
		}
	unsigned m_uIndex;
	unsigned m_uSize;
	double m_dLength;
	ClustNode *m_ptrLeft;
	ClustNode *m_ptrRight;
	ClustNode *m_ptrParent;
	ClustNode *m_ptrNextCluster;
	ClustNode *m_ptrPrevCluster;
	unsigned *m_uLeafIndexes;
	};

class Clust
	{
public:
	Clust();
	virtual ~Clust();

	void Create(ClustSet &Set, CLUSTER Method);

	unsigned GetLeafCount() const;

	unsigned GetClusterCount() const;
	unsigned GetClusterSize(unsigned uNodeIndex) const;
	unsigned GetLeaf(unsigned uClusterIndex, unsigned uLeafIndex) const;

	unsigned GetNodeCount() const { return 2*m_uLeafCount - 1; }
	const ClustNode &GetRoot() const { return m_Nodes[GetRootNodeIndex()]; }
	unsigned GetRootNodeIndex() const { return m_uNodeCount - 1; }

	const ClustNode &GetNode(unsigned uNodeIndex) const;
	bool IsLeaf(unsigned uNodeIndex) const;
	unsigned GetLeftIndex(unsigned uNodeIndex) const;
	unsigned GetRightIndex(unsigned uNodeIndex) const;
	double GetLength(unsigned uNodeIndex) const;
	double GetHeight(unsigned uNodeIndex) const;
	const char *GetNodeName(unsigned uNodeIndex) const;
	unsigned GetNodeId(unsigned uNodeIndex) const;

	JOIN GetJoinStyle() const { return m_JoinStyle; }
	LINKAGE GetCentroidStyle() const { return m_CentroidStyle; }

	void SetDist(unsigned uIndex1, unsigned uIndex2, double dDist);
	double GetDist(unsigned uIndex1, unsigned uIndex2) const;

	// void ToPhylip(Phylip &tree);

//private:
	void SetLeafCount(unsigned uLeafCount);

	void CreateCluster();
	void JoinNodes(unsigned uLeftNodeIndex, unsigned uRightNodeIndex, 
	  double dLeftLength, double dRightLength, unsigned uNewNodeIndex);

	void ChooseJoin(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
	  double *ptrdLeftLength, double *ptrdRightLength);
	void ChooseJoinNeighborJoining(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
	  double *ptrdLeftLength, double *ptrdRightLength);
	void ChooseJoinNearestNeighbor(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
	  double *ptrdLeftLength, double *ptrdRightLength);

	double ComputeDist(unsigned uNewNodeIndex, unsigned uNodeIndex);
	double ComputeDistAverageLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex);
	double ComputeDistMinLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex);
	double ComputeDistMaxLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex);
	double ComputeDistNeighborJoining(unsigned uNewNewIndex, unsigned uNodeIndex);
	double ComputeDistMAFFT(unsigned uNewNewIndex, unsigned uNodeIndex);

	double Calc_r(unsigned uNodeIndex) const;

	unsigned VectorIndex(unsigned uIndex1, unsigned uIndex2) const;

	unsigned GetFirstCluster() const;
	unsigned GetNextCluster(unsigned uNodeIndex) const;

	double ComputeMetric(unsigned uIndex1, unsigned uIndex2) const;
	double ComputeMetricNearestNeighbor(unsigned i, unsigned j) const;
	double ComputeMetricNeighborJoining(unsigned i, unsigned j) const;

	void InitMetric(unsigned uMaxNodeIndex);
	void InsertMetric(unsigned uIndex1, unsigned uIndex2, double dMetric);
	double GetMinMetric(unsigned *ptruIndex1, unsigned *ptruIndex2) const;
	double GetMinMetricBruteForce(unsigned *ptruIndex1, unsigned *ptruIndex2) const;
	void DeleteMetric(unsigned uIndex);
	void DeleteMetric(unsigned uIndex1, unsigned uIndex2);
	void ListMetric() const;

	void DeleteFromClusterList(unsigned uNodeIndex);
	void AddToClusterList(unsigned uNodeIndex);

	/*
	void RBDelete(unsigned RBNode);
	unsigned RBInsert(unsigned i, unsigned j, double fMetric);

	unsigned RBNext(unsigned RBNode) const;
	unsigned RBPrev(unsigned RBNode) const;
	unsigned RBMin(unsigned RBNode) const;
	unsigned RBMax(unsigned RBNode) const;

	void ValidateRB(const char szMsg[] = 0) const;
	void ValidateRBNode(unsigned Node, const char szMsg[]) const;
	*/

//private:
	JOIN m_JoinStyle;
	LINKAGE m_CentroidStyle;
	ClustNode *m_Nodes;
	unsigned *m_ClusterIndexToNodeIndex;
	unsigned *m_NodeIndexToClusterIndex;
	unsigned m_uLeafCount;
	unsigned m_uNodeCount;
	unsigned m_uClusterCount;
	unsigned m_uTriangularMatrixSize;
	double *m_dDist;
	ClustSet *m_ptrSet;
	ClustNode *m_ptrClusterList;
	};

class ClustSet
	{
public:
	virtual unsigned GetLeafCount() = 0;
	virtual double ComputeDist(const Clust &C, unsigned uNodeIndex1,
	  unsigned uNodeIndex2) = 0;
	// An object of class ClustSetNMF is instantiated in Evaluation::IndHC() and in Evaluation::HierarchCustering()
	// CLANG is complaining that an Abstract Class cannot be instantiated.
	// To avoid this problem, comment out the declaration of the virtual function JoinNodes().
	// This will prevent ClustSetNMF from being an Abstract Class.
//	virtual void JoinNodes(const Clust &C, unsigned uLeftNodeIndex,
//	  unsigned uRightNodeIndex, unsigned uJoinedNodeIndex,
//	  double *ptrdLeftLength, double *ptrdRightLength) = 0;
	virtual const char *GetLeafName(unsigned uNodeIndex) = 0;
	virtual unsigned GetLeafId(unsigned uNodeIndex) = 0;
	};

#endif // Clust_h

