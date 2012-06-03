#pragma once
#include "Update.h"

class KL :
	public Update
{
public:

	KL(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~KL(void) {}

private:
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	void   SaveIntermediaFiles() {}

};
