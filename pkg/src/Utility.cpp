#include "Utility.h"

using namespace std;

long Utility::GetSeed(long baseSeed, int rank, int chain, char fileType)
{
	return baseSeed + rank + chain + (int) fileType;
}

double Utility::PowerUp(double base, double power)
{
	if(power == 1.0)
	{
		return base;
	}
	else if(power == 2.0)
	{
		return base * base;
	}
	else if(power == 0.5)
	{
		return sqrt(base);
	}
	else
	{
		return exp( power * log(base) );
	}
}

void Utility::GetDataInfo(const array1d prf, double &sum, double &var)
{
	sum = 0.0;
	var = 0.0;
	for (unsigned i=0; i<prf.size(); i++)
	{
		sum += prf[i];
		var += prf[i] * prf[i];
	}
	var = sqrt( prf.size() * var - sum * sum );
}

double Utility::PearsonCorrelation(const array1d &vector1, const array1d &vector2)
{
	if(vector1.size() != vector2.size())
	{
//		cout << "Two vectors have different sizes in Utility::PearsonCorrelation!" << endl;
		return(0.0);
	}

	double qSum, qVar, tSum, tVar;

	Utility::GetDataInfo( vector1, qSum, qVar );
	Utility::GetDataInfo( vector2, tSum, tVar );

	double product = 0.0;
	for(unsigned i=0; i<vector1.size(); i++) 
	{
		product += vector1[i] * vector2[i];
	}

	if(qVar == 0.0 || tVar == 0.0)
	{
//		cout << "Utility::PearsonCorrelation has zero value in the denominator!" << endl;
		return(0.0);
	}
	else
	{
		return ( vector1.size() * product - qSum * tSum ) / (qVar * tVar);
	}
}

int Utility::Combination42(const int size)
{
	if( size <= 1)
	{
		return 0;
	}
	else
	{
		return size*(size-1)/2;
	}
}

// standardize the factorized matrix into a comparable scale, make cluster assignment based on that
int Utility::GetZscores(const array1d &orig, array1d &zscore)
{
	double sum1 = 0.0;
	double sum2 = 0.0;

	if(orig.size() <= 1) 
	{
		return(-1);
	}

	for (unsigned i=0; i<orig.size(); i++)
	{
		sum1 += orig[i];
		sum2 += orig[i] * orig[i];
	}

	double sd = (sum2 - (sum1 * sum1) / orig.size()) / (orig.size() - 1);

	// in case all samples in the data set are more or less the same, sd ~ 0.0
	if(fabs(sd) < EPS || sd < 0.0) 
	{
		return -1;
	}

	sd = sqrt(sd);
	double average = sum1 / orig.size();

	for (unsigned i=0; i<orig.size(); i++) 
	{
		zscore[i] = (orig[i] - average) / sd;
	}

	return(1);
}


///////////////////////////////////////////////////////////////////////////////////
// ran2 provided a better portable random number generators (0,1), which
// library function rand() might not.
// For using rand(), only the first 1/10 of RAND_MAX numbers are
// safely random. For most of the current machines, RAND_MAX = 2147483647,
// which means about 150 millions numbers generated from rand() would be good. 
///////////////////////////////////////////////////////////////////////////////////
double Utility::ran2(long *idum)
{
	int   j;
	long  k;
	static long idum2 = 123456789;
	static long iy = 0;
	static int    IM1 = 2147483563;
	static int    IMM1 = IM1 - 1;
	static int    IM2 = 2147483399;
	static double AM  = 1.0 / (double)IM1;
	static int    IA1 = 40014;
	static int    IA2 = 40692;
	static int    IQ1 = 53668;
	static int    IQ2 = 52774;
	static int    IR1 = 12211;
	static int    IR2 = 3791;
	const  int    NTAB = 32;
	static int    NDIV = 1 + IMM1 / NTAB;
	static double EPS = 1.2e-7;
	static double RNMX = 1.0 - EPS;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) 
	{
		if (-(*idum) < 1) *idum = 1; else *idum = -1 * (*idum);
		idum2 = (*idum);

		for (j=NTAB+7; j>=0; j--)
		{
			k = (*idum) / IQ1;
			*idum = IA1 * (*idum - k * IQ1) - k * IR1;

			if(*idum < 0) *idum += IM1;
			if(j < NTAB) iv[j] = *idum;
		}

		iy = iv[0];
	}

	k = (*idum) / IQ1;
	*idum = IA1 * (*idum - k * IQ1) - k * IR1;

	if(*idum < 0) *idum += IM1;

	k = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;

	if(idum2 < 0) idum2 += IM2;

	j = iy / NDIV;
	iy = iv[j] - idum2;
	iv[j] = *idum;

	if(iy < 1) iy += IMM1;
	if( (temp = AM * iy) > RNMX) return RNMX; else return temp;
}

// void Tokenize(string &str, vector<string> &tokens, const string &delimiters)
void Utility::Tokenize(const string &str, vector<string> &tokens, const string &delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

