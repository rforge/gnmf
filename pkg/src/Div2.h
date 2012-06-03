#pragma once
#include "Update.h"

class Div2 :
	public Update
{
public:
	Div2(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~Div2(void) {}

private:
	void UpdatePatternMatrix(WeightingMatrix &);
	void UpdateAmplitudeMatrix(WeightingMatrix &);
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	void   SaveIntermediaFiles() {}
	int    CheckConvergency( int round = 0 );
};
