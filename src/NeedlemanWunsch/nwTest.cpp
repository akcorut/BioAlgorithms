#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <map>
#include <iomanip>

using namespace std;

void createSimilarityMatrix(string seqFirst, string seqSecond)
{
	seqFirst = seqFirst;
	seqSecond = seqSecond;

	lengthFirst = seqFirst.length();
	lengthSecond = seqSecond.length();

	similarity = 0;

	sMatrix = new int*[lengthFirst + 2];
	for (int i = 0; i < (lengthFirst + 1); i++)
	{
		sMatrix[i] = new int[lengthSecond + 2];
	}

	for (int i = 0; i < lengthFirst; i++)
	{
		for (int j = 0; j < lengthFirst; j++)
		{
			sMatrix[i][j] = 0;
		}
	}
	
}

int main(int argc, char **argv) 
{
	// File Parsing.
	if (argc > 1) {
        cout << "First fasta file: " << argv[1] << endl;
        cout << "Second fasta file: " << argv[2] << endl; 

    } else {
        cout << "No file detected. Please provide the proper paths of the input files.";
        return -1;
    }

	
	string seqFirst, seqSecond, line; // Declare string objects for reading the files

	ifstream fastaFirst(argv[1]); // Read the first file from command line
	for (int i = 0; getline(fastaFirst, line); i++)
	{
		if (i != 0) 
		{
			seqFirst += line;
		}
	}

    ifstream fastaSecond(argv[2]); // Read the second file from command line
	for (int i = 0; getline(fastaSecond, line); i++)
	{
		if (i != 0)
		{
			seqSecond += line;
		}
	}
    
	cout << seqFirst << endl;
    cout << seqFirst.length() << endl;
    cout << seqSecond << endl;
    cout << seqSecond.length() << endl;
	
	
	createSimilarityMatrix(seqFirst, seqSecond);
	return 0;
}