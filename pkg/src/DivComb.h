#pragma once
#include "Update.h"

class DivComb :
	public Update
{
public:
	DivComb(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{}

	~DivComb(void) {}



private:
	void UpdateAmplitudeMatrix();
	void UpdatePatternMatrix();
	void SaveIntermediaFiles() {}
	void UpdatePatternMatrix( WeightingMatrix &  );
	void UpdateAmplitudeMatrix( WeightingMatrix & );
};
