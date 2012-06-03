#pragma once
#include "Update.h"

class Gamma :
	public Update
{
public:

	Gamma(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~Gamma(void) {}

private:
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	int    CheckConvergency( int round = 0 );
	void   SaveIntermediaFiles() {}

};
