#ifndef UTILITY_H_
#define UTILITY_H_

#include "Parameters.h"
#include <math.h>
#include <functional>
#include <limits> // get the largest and smallest possible values for a numeric type in the platform
#include <exception>

template <class T>
class SortVectorIndex : public std::binary_function<int, int, bool>
{
	const std::vector<T> &myVector;

public:
	// constructor which take a reference to a vector
	SortVectorIndex(const std::vector<T> &_vector) : myVector(_vector) {}

	// comparison operator. this can be called by std:sort with two integers
	bool operator() (int i, int j) const
	{
		return myVector[j] < myVector[i];
	}
}; // class SortVectorIndex;

class Utility
{
public:
	Utility() {}
	~Utility() {}

	static double PowerUp(double base, double power);
	static double ran2(long *);
	static double PearsonCorrelation(const array1d &, 
		const array1d &
		);
	static int    Combination42(const int size);
	static void	  Tokenize(const std::string &str, 
		std::vector<std::string> &tokens, 
		const std::string &delimiters = " "
		);
	static long GetSeed(long baseSeed, int rank, int chain, char fileType);
	static void   GetDataInfo(const array1d prf, 
		double &sum, double &var
		);
	static int GetZscores(const array1d &orig, array1d &zscore);

	template <class T>
	static void GetVectorMaxValue(const std::vector<T> &vec, T &ele, int &index)
	{
		// ele = static_cast<T>(-1000000.0);
		ele = vec[0];
		// ele = std::numeric_limits<T>::max();
		for (unsigned i=0; i<vec.size(); i++)
		{
			if(vec[i] >= ele)
			{
				ele = vec[i];
				index = i;
			}
		}
	}

	template <class T>
	static void GetVectorMaxValue(const std::vector<T> &vec, int &index)
	{
		// T ele = std::numeric_limits<T>::max();
		// T ele = static_cast<T>(-1000000.0);
		T ele = vec[0];
		for (unsigned i=0; i<vec.size(); i++)
		{
			if(vec[i] >= ele)
			{
				ele = vec[i];
				index = i;
			}
		}
	}

	/*
	template <class T>
	static void ReOrderVector(std::vector<T> &vec, const std::vector<int> &sequence)
	{
		// generate temp vector based on which vec would be sorted
		std::vector<int> temp( sequence.size() );
        for (int i=0; i<sequence.size(); i++) temp[sequence[i]] = i;

		// sort vec based on temp, e.g. temp[5] = 1, temp[1] = 5, then vec[1] = vec0[5], vec[5] = vec0[1]
		std::sort(vec.begin(), vec.end(), SortVectorIndex(temp));
	}

	template <class T>
	static void ReOrderVector(std::vector< vector<T> > &vec, const std::vector<int> &sequence)
	{
		// generate temp vector based on which vec would be sorted
		std::vector<int> temp( sequence.size() );
        for (int i=0; i<sequence.size(); i++) temp[sequence[i]] = i;

		// sort vec's row based on temp, e.g. temp[5] = 1, temp[1] = 5, then vec[1] = vec0[5], vec[5] = vec0[1]
		std::sort(vec.begin(), vec.end(), SortVectorIndex(temp));

		// change the value in upper-right triangle of vec correspondingly
		for (unsigned i=0; i<vec.size(); i++)
		{
			for (unsigned j=i; j<vec.size(); j++)
			{
				vec[i][j] = vec[j][i];
			}
		}
	}
	*/
	
	template <class T>
	static void ReOrderVector(std::vector<T> &vec, const std::vector<int> &sequence)
	{
		// create a temp vector holding original value of vec
		std::vector<T> temp( vec.size() );
		temp.assign( vec.begin(), vec.end() );

		// fill in value to vec according to sequence
		for (unsigned i=0; i<sequence.size(); i++)
		{
			vec[i] = temp[ sequence[i] ];
		}
	}

	/*
	template <class T>
	static void ReOrderVector(std::vector< vector<T> > &vec, const std::vector<int> &sequence)
	{
		// generate temp vector holding original value of vec
		std::vector< vector<T> > temp;
		temp.resize( vec.size() );
		for (unsigned i=0; i<temp.size(); i++)
		{
			temp[i].resize( vec[i].size() );
			temp[i].assign( vec[i].begin(), vec[i].end() );
		}

        for (unsigned i=0; i<sequence.size(); i++)
		{
			unsigned newIndex = sequence[i];
			vec[i].assign( temp[newIndex].begin(), temp[newIndex].end() );
		}

		// change the value in upper-right triangle of vec correspondingly
		for (unsigned i=0; i<vec.size(); i++)
		{
			for (unsigned j=i; j<vec.size(); j++)
			{
				vec[i][j] = vec[j][i];
			}
		}
	}
	*/
};



#endif

