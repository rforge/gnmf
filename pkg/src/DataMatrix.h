#ifndef DATAMATRIX_H_
#define DATAMATRIX_H_

#include "Parameters.h"
#include <fstream>
#include <sstream>
#include <math.h>
#include <limits>
#include <vector>
#include "Utility.h"
#include <map>

class DataMatrix
{
	friend std::istream &operator >> (std::istream &in, DataMatrix &myMatrix);
	friend std::ostream &operator << (std::ostream &out, const DataMatrix &myMatrix);
public:
	int    nRows, nColumns;
	bool   HD;
	array2d data;
	array2d hdIndex;
	array2d hdValue;
	std::vector< std::string > rowNames, columnNames;
	double average, stdd, average2;
	SSMap properties;

	DataMatrix() 
	{
		HD = false;
	}
	DataMatrix(const int nR, const int nC);
	DataMatrix(const DataMatrix &_data);

	int  GetZscores(); 
	int  GetStatistical();
	void Transposing(array2d &);
	void Transposing();
	void Scaling( const double );
	int  Scaling( const array2d & );
	int  DScaling( const array2d & );
	int TakePower(const double alpha);
	void Adding( const array2d &);
	void SetProperties( SSMap & );
	SSMap & GetProperties();


	~DataMatrix();

private:
	// const int BUFFSIZE = 20000;
};


#endif

