#pragma once
#include "Update.h"

class Renyi :
	public Update
{
public:
	Renyi(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{
	}

	~Renyi(void) {}

private:
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	int    CheckConvergency( int round = 0 );
	void   SaveIntermediaFiles() {}
};
