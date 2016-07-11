#pragma once
#include "Update.h"

// JMM (7/11/2016): Include R header file AFTER system headers. 
// Reference: First few paragraphs in in Section 6 The R API: entry points for C code
// in the online documentation for "Writing R Extensions".
// Also see email from Prof. B. Ripley sent 7/10/2016
// at 3:58 AM entitled "CRAN packages failing to install with modern C++"
// Note: Inclusion of R headers used to be done in the GNMF header files,
// e.g. in 
#define R_NO_REMAP
#include "R.h"           // R functions

class BD :
	public Update
{
public:

	BD(	
		const VMatrix &_data, 
		const ParaControler &_control
		) 
		:	
		Update(_data, _control)
	{
		DataMatrix tt(_data);
		tt.GetStatistical();
		double tmp = tt.average * tt.average;
		bdn = tmp / ( tt.average + tmp - tt.average2);
	}

	~BD(void) {}

private:
	double bdn;
	void   UpdateAmplitudeMatrix();
	void   UpdatePatternMatrix();
	int    CheckConvergency( int round = 0 );
	void   SaveIntermediaFiles() {}

};
