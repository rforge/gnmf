#pragma once
#include "DataMatrix.h"

class CMatrix :
	public DataMatrix
{
public:
	CMatrix(void) {}
	CMatrix(int nComponents) 
		: DataMatrix(nComponents, nComponents) 
	{}

	~CMatrix(void) {}

	void ReadTriangle( std::string fileName);
};
