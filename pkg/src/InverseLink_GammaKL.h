#pragma once
#include "Update.h"

class InverseLink_GammaKL :
	public Update
{
public:

	InverseLink_GammaKL(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~InverseLink_GammaKL(void) {}

private:
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	int    CheckConvergency( int round = 0 );
	void   SaveIntermediaFiles() {}
};
