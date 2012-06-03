#ifndef _CLUSTDATASET_H
#define _CLUSTDATASET_H

#include "Parameters.h"
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>

struct member
{
	member() : index(0), name("")
	{
	}

	member(int _index, std::string _name)
	{
		index = _index;
		name = _name;
	}

	int index;
	std::string name;

	bool operator()( member *m1, member *m2 )
	{
		return (m1->index < m2->index);
	}
};

class ClustDataSet
{
	friend std::istream &operator >> (std::istream &in, ClustDataSet &myClust);
	friend std::ostream &operator << (std::ostream &out, const ClustDataSet &myClust);

public:
	int nClusters;
	std::vector< std::vector<member *> > members;

	ClustDataSet() {}
	ClustDataSet( const std::vector< std::vector<int> > &_members, const std::vector< std::string > &names );

	~ClustDataSet() 
	{
		members.clear();
		/*
		for (int i=0; i<members.size(); i++)
		{
			delete members[i];
		}
		*/
	}
};

#endif
