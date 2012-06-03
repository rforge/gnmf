#pragma once
#include "Update.h"

class NBD :
	public Update
{
public:

	NBD(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{
		DataMatrix tt(_data);
		tt.GetStatistical();
		double tmp = tt.average * tt.average;
		bdn = tmp / ( tt.average2 - tt.average - tmp);
	}

	~NBD(void) {}


private:
	double bdn;
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	int    CheckConvergency( int round = 0 );
	void   SaveIntermediaFiles() {}

};
