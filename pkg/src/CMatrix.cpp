#include "CMatrix.h"

using namespace std;

void CMatrix::ReadTriangle( string fileName)
{
	nRows = 0;
	nColumns = 0;

	ifstream iData(fileName.c_str(), ios::in);
	while( ! iData.eof() )
	{
		string line;
		// string line = iData.getline();
		getline(iData, line);

		vector<string> tokens;
		Utility::Tokenize(line, tokens);
	
		rowNames.push_back( tokens[0] );
	
		array1d rowData;
		// lower triangle copied from upper triangle
		for (unsigned i=0; i<data.size(); i++)
		{
			rowData.push_back( data[i][ nRows ] );
		}
		// diagonal element
		rowData.push_back( 0.0 ); 
		// upper triangle read from file
		for (unsigned i=1; i<tokens.size(); i++)
		{
			rowData.push_back( atof( tokens[i].c_str() ) );
		}
		data.push_back( rowData );

		nRows = data.size();
	}

	columnNames = rowNames;
	nColumns = nRows;
}

