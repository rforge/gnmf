#include "ClustDataSet.h"

using namespace std;

ClustDataSet::ClustDataSet(
		const vector< vector<int> > &indexes, 
		const vector< string > &names
		)
{
	nClusters = indexes.size();
	members.resize(nClusters);
	for (int i=0; i<nClusters; i++)
	{
		for (int j=0; j<(int)indexes[i].size(); j++)
		{
			members[i].push_back( new member(indexes[i][j], names[indexes[i][j]]) );
		}
	}
}

istream &operator >> (std::istream &in, ClustDataSet &myClust)
{
	// new input
	int clustID, dataIndex;
	std::string dataName;
	myClust.members.resize(10);
	myClust.nClusters = 0;
	string line;
	// while(! in.eof())
	while(getline(in, line))
	{
		if(line.length() < 2) break;
		istringstream iss(line);
		iss >> dataName >> dataIndex >> clustID;
		// in >> dataName >> dataIndex >> clustID;
		// std::cout << dataName << " " << dataIndex << " " << clustID << endl;

		if(myClust.nClusters < clustID) myClust.nClusters = clustID;
		clustID --; // usually the input is using index starting from 1
		dataIndex --;

		if(clustID >= (int)myClust.members.size()) myClust.members.resize( clustID * 2 );
		myClust.members[clustID].push_back( new member(dataIndex, dataName) );
		// cout << myClust.nClusters  << "  clt" << clustID << " size " << myClust.members[clustID].size() << endl;
	}
	// cout << myClust.nClusters << endl;
	myClust.members.resize( myClust.nClusters );

	// sort group memeber
	struct member sp;
	for (unsigned i=0; i<myClust.members.size(); i++)
	{
		// sort(myClust.members[i].begin(), myClust.members[i].end(), sp);
		sort(myClust.members[i].begin(), myClust.members[i].end(), sp);
	}

	return (in);
}

ostream &operator << (std::ostream &out, const ClustDataSet &myClust)
{
	// new output
	for (unsigned i=0; i<(unsigned)myClust.nClusters; i++)
	{
		for (unsigned j=0; j<myClust.members[i].size(); j++)
		{
			out << myClust.members[i][j]->name << "\t" << myClust.members[i][j]->index + 1 << "\t" << i+1 << std::endl;
		}
	}

	return (out);
}
