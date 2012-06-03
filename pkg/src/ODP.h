#pragma once
#include "Update.h"

class ODP :
	public Update
{
public:

	ODP(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~ODP(void) {}

private:
	std::vector< double > phi;
	int CheckConvergency(int round = 0);
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	void   SaveIntermediaFiles();
};
