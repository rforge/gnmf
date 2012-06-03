#include "DataMatrix.h"

using namespace std;

const int BUFFSIZE = 200000;

DataMatrix::DataMatrix(const int nR, const int nC)
{
	HD = false;
	nRows = nR;
	nColumns = nC;

	// reserve fixed amount memory for the data members
	columnNames.resize(nColumns);
	rowNames.resize(nRows);

	data.resize( nRows, array1d(nColumns) );
}

DataMatrix::DataMatrix(const DataMatrix &_data)
{
	nRows = _data.nRows;
	nColumns = _data.nColumns;
	data = _data.data;
	rowNames = _data.rowNames;
	columnNames = _data.columnNames;

	properties = _data.properties;

	HD = _data.HD;
	hdIndex = _data.hdIndex;
	hdValue = _data.hdValue;
}

DataMatrix::~DataMatrix()
{
	data.clear();
	rowNames.clear();
	columnNames.clear();
}

ostream& operator << (ostream &out, const DataMatrix &myData)
{
	SSMap::const_iterator ip;
	for (ip=myData.properties.begin(); ip!=myData.properties.end(); ++ip)
	{
		out << "# " << ip->first << " " << ip->second << endl;
	}

	for (int i=0; i<(int)myData.columnNames.size(); i++)
	{
		out << "\t" << myData.columnNames[i];
	}
	out << endl;

	for (int i=0; i<(int)myData.rowNames.size(); i++)
	{
		out << myData.rowNames[i];

		for (int j=0; j<(int)myData.columnNames.size(); j++)
		{
			out << "\t" << myData.data[i][j];
		}
		out << endl;
	}

	return (out);
}

istream& operator >> (istream &in, DataMatrix &myData)
{
	// clear all memories
	myData.data.clear();
	myData.columnNames.clear();
	myData.rowNames.clear();

	char buff[BUFFSIZE];
	char dataUnit[50];
	stringstream line;
	while( in.getline( buff, BUFFSIZE, '\n' ) )
	{
		line << buff;
		if(buff[0] == '#')
		{
			string tmp, name, value;
			line >> tmp >> name >> value;
			myData.properties[name] = value;
		}
		else if(buff[0] == '\t')
		// else if(myData.columnNames.size() < 1)
		{
			while( line.getline(dataUnit, 50, '\t') )
			{
				myData.columnNames.push_back( dataUnit );
			}
			myData.columnNames.erase( myData.columnNames.begin() );
		}
		else
		{
			string rowName;
			line >> rowName;
			myData.rowNames.push_back( rowName );
			if(myData.HD)
			{
				int nTerms;
				line >> nTerms;
				array1d rowIndexes( nTerms ), rowValues( nTerms );
				string term;
				for (int i=0; i<nTerms; i++)
				{
					line >> term;
					rowIndexes[i] = atoi( term.substr(0, term.find(":")).c_str() );
					rowValues[i] = atof( term.substr(term.find(":") + 1).c_str() );
				}
				myData.hdIndex.push_back( rowIndexes );
				myData.hdValue.push_back( rowValues );
			}
			else
			{
				array1d rowData( myData.columnNames.size() );
				for (unsigned i=0; i<myData.columnNames.size(); i++)
				{
					line >> rowData[i];
				}
				myData.data.push_back(rowData);
			}
		}
		line.clear();
	}

	myData.nRows = myData.rowNames.size();
	myData.nColumns = myData.columnNames.size();

	return (in);
}

int DataMatrix::GetStatistical()
{
	double sum1 = 0.0;
	double sum2 = 0.0;

	if(nRows <= 1) 
	{
		return(-1);
	}

	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			sum1 += data[i][j];
			sum2 += data[i][j] * data[i][j];
		}
	}

	int nSampleSize = nRows * nColumns;
	average2 = sum2 / nSampleSize; 
	stdd = (sum2 - (sum1 * sum1) / nSampleSize) / (nSampleSize - 1);

	if(stdd <= 0.0) stdd = EPS; 

	stdd = sqrt(stdd);
	average = sum1 / nSampleSize;

	return 1;
}

int DataMatrix::GetZscores()
{
	// GetStatistical();

	for (int i=0; i<nRows; i++) 
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] = (data[i][j] - average) / stdd;
		}
	}

	return(1);
}

int DataMatrix::TakePower(const double alpha)
{
	if(alpha == 1.0) return 1;

	if(alpha == 2.0)
	{
		for (int i=0; i<nRows; i++)
		{
			for (int j=0; j<nColumns; j++)
			{
				data[i][j] *= data[i][j];
			}
		}
	}
	else if(alpha == 0.5)
	{
		for (int i=0; i<nRows; i++)
		{
			for (int j=0; j<nColumns; j++)
			{
				data[i][j] = sqrt(data[i][j]);
			}
		}
	}
	else
	{
		for (int i=0; i<nRows; i++)
		{
			for (int j=0; j<nColumns; j++)
			{
				// try
				{
					data[i][j] = exp( alpha * log(data[i][j]) );
				}
				/* catch( std::runtime_error &e )
				{
					std::cerr << "DataMatrix::TakePower(const double alpha): " << e.what() << endl;
					return -1;
				}
				*/
			}
		}
	}

	return 1;
}


void DataMatrix::Transposing(array2d &outputMatrix)
{
	// allocate memory for the output matrix
	outputMatrix.resize( nColumns, array1d(nRows) );

	for (int i=0; i<nColumns; i++)
	{
		for (int j=0; j<nRows; j++)
		{
			outputMatrix[i][j] = data[j][i];
		}
	}
}

void DataMatrix::Transposing()
{
	int nTemp = nRows;
	nRows = nColumns;
	nColumns = nTemp;

	// allocate memory for the output matrix
	array2d tmpMatrix(nRows, array1d(nColumns));

	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			tmpMatrix[i][j] = data[j][i];
		}
	}

	data.resize( nRows, array1d(nColumns) );
	data = tmpMatrix;

	vector<string> names( rowNames );
	rowNames.resize(nRows);
	rowNames = columnNames;
	columnNames.resize(nColumns);
	columnNames = names;
}


void DataMatrix::Adding(const array2d &vec)
{
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] += vec[i][j];
		}
	}
}

void DataMatrix::Scaling(const double scaler)
{
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] *= scaler;
		}
	}
}



int DataMatrix::Scaling(const array2d &scaler)
{
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] *= scaler[i][j];
		}
	}

	return 1;
}


int DataMatrix::DScaling(const array2d &scaler)
{
	for (int i=0; i<nRows; i++)
	{
		for (int j=0; j<nColumns; j++)
		{
			data[i][j] /= scaler[i][j];
		}
	}

	return 1;
}

void DataMatrix::SetProperties( SSMap & maps)
{
	properties = maps;
}

SSMap & DataMatrix::GetProperties()
{
	SSMap *maps = 0;

	return *maps;
}
