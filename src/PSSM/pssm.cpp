#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <map>
#include <iomanip>
#include <stdio.h>
#include <vector>

using namespace std;

void freqMatrix(vector <string> motifSeq, vector <string> dnaSeq)
{
	motifSeq = motifSeq;
	dnaSeq = dnaSeq;	
	//int fMatrix[motifSeq[0].size()+1][5];
	int** fMatrix = new int*[motifSeq[0].size()+1];
	for(int i = 0; i < motifSeq[0].size()+1; ++i)
    	fMatrix[i] = new int[5];

	for(int i=1; i<=motifSeq[0].size(); i++){
    	fMatrix[i][0] = i;
  	}
  	
	fMatrix[0][1] = 'A';
	fMatrix[0][2] = 'T';
	fMatrix[0][3] = 'G';
	fMatrix[0][4] = 'C';

	for(int j=0; j<motifSeq.size();j++){
		for (int i=0; i<motifSeq[0].size();i++){
			if(toupper(motifSeq[j][i])=='A'){
				fMatrix[i+1][1]+=1;
			}else if(toupper(motifSeq[j][i])=='T'){
				fMatrix[i+1][2]+=1;
			}else if(toupper(motifSeq[j][i])=='G'){
				fMatrix[i+1][3]+=1;
			}else if(toupper(motifSeq[j][i])=='C'){
				fMatrix[i+1][4]+=1;	
			}
		}
	}

	for(int j=0; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			cout << fMatrix[j][i] << "\t";
			if (i==4)
			{
				cout<< "\n";
			}
		}
	}
	
}

int main(int argc, char **argv) 
{	
	int cutoff;
  	
	cout<<"Enter the minimum score cutoff: ";
  	cin>> cutoff;

	// File Parsing.
	if (argc > 1) {
        cout << "Motif sequence file: " << argv[1] << endl;
        cout << "DNA fasta file: " << argv[2] << endl; 

    } else {
        cout << "No file detected. Please provide the proper paths of the input files.";
        return -1;
    }

	vector <string> motifSeq;
	vector <string> dnaSeq;
	string line; // Declare the required string objects for reading the files


	ifstream fastaFirst(argv[1]); // Read the first file from command line
	for (int i = 0; getline(fastaFirst, line); i++)
	{
			motifSeq.push_back(line);
	}

    ifstream fastaSecond(argv[2]); // Read the second file from command line
	for (int i = 0; getline(fastaSecond, line); i++)
	{
		if (i != 0)
		{
			dnaSeq.push_back(line);
		}
	}
    
    cout << dnaSeq.size() << endl;
	cout << motifSeq[0].size() << endl;

	cout<<"The dimension of the motif data is ";
    cout<<motifSeq.size()<<" x ";
    cout<<motifSeq[0].size()<<endl;

	freqMatrix(motifSeq, dnaSeq);
    return 0;
}