#pragma once
#include "Update.h"

class GammaKL :
	public Update
{
public:

	GammaKL(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~GammaKL(void) {}

private:
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	int    CheckConvergency( int round = 0 );
	void   SaveIntermediaFiles() {}
};
