#pragma once
#include "Update.h"

class ED :
	public Update
{
public:

	ED(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~ED(void) {}

private:
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	int    CheckConvergency( int round = 0 );
	void   SaveIntermediaFiles() {}

};
