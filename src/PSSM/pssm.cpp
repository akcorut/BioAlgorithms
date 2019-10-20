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
#include <math.h>

using namespace std;

double getGC(vector <string> dnaSeq){
	dnaSeq = dnaSeq;

	double GCcount=0;
	for(int i=0;i<dnaSeq.size();i++){
		for(int j=0;j<dnaSeq[0].size();j++){
			if(dnaSeq[i][j]=='G' || dnaSeq[i][j]== 'C'){
				GCcount += 1;
			}
		}
	}
	double sizeGeno=dnaSeq.size()*dnaSeq[0].size();
	double GCfreq = (GCcount/sizeGeno);

	return GCfreq;
}

void freqMatrix(vector <string> motifSeq, vector <string> dnaSeq)
{
	motifSeq = motifSeq;
	dnaSeq = dnaSeq;	
	double fMatrix[motifSeq[0].size()+1][5];
	/*double** fMatrix = new double*[motifSeq[0].size()+1];
	for(int i = 0; i < motifSeq[0].size()+1; ++i)
    	fMatrix[i] = new double[5];*/

	// Frequency matrix
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

	/*for(int j=0; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			cout << fMatrix[j][i] << "\t";
			if (i==4)
			{
				cout<< "\n";
			}
		}
	}*/
	
	// Add pseudocounts
	double pseCount = 0.25;

	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=1; i<=4;i++){
			fMatrix[j][i] += pseCount;
		}
	}

	// Probability matrix
	double pMatrix[motifSeq[0].size()+1][5];
	for(int i=1; i<=motifSeq[0].size(); i++){
    	pMatrix[i][0] = i;
  	}
	pMatrix[0][1] = 'A';
	pMatrix[0][2] = 'T';
	pMatrix[0][3] = 'G';
	pMatrix[0][4] = 'C';

	//int numSeq = motifSeq[0].size()*motifSeq.size();
	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=1; i<=4;i++){
			pMatrix[j][i] = fMatrix[j][i]/motifSeq.size();
		}
	}	
	
	/*for(int j=0; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			cout << pMatrix[j][i] << "\t";
			if (i==4)
			{
				cout<< "\n";
			}
		}
	}*/

	/*double GCcount=0;
	for(int i=0;i<dnaSeq.size();i++){
		for(int j=0;j<dnaSeq[0].size();j++){
			if(dnaSeq[i][j]=='G' || dnaSeq[i][j]== 'C'){
				GCcount += 1;
			}
		}
	}
	double sizeGeno=dnaSeq.size()*dnaSeq[0].size();
	cout << endl << GCcount;
	cout << endl << sizeGeno;
	double GCfreq = (GCcount/sizeGeno);*/
	
	double GCfreq = getGC(dnaSeq);
	//cout << endl << GCfreq;
	double freqG, freqC, freqA, freqT;
	
	freqG=GCfreq/2;
	freqC=GCfreq/2;
	freqA=(1-GCfreq)/2;
	freqT=(1-GCfreq)/2;

	//cout << endl << freqG << endl;
	//cout << endl << freqC << endl;
	//cout << endl << freqA << endl;
	//cout << endl << freqT << endl;

	double PSSM[motifSeq[0].size()+1][5];
	
	for(int i=1; i<=motifSeq[0].size(); i++){
		PSSM[i][0] = i;
	}
	
	PSSM[0][1] = 'A';
	PSSM[0][2] = 'T';
	PSSM[0][3] = 'G';
	PSSM[0][4] = 'C';

	for(int j=1; j<=motifSeq[0].size();j++){
		PSSM[j][1] = log(pMatrix[j][1]/freqA);
		PSSM[j][2] = log(pMatrix[j][2]/freqT);
		PSSM[j][3] = log(pMatrix[j][3]/freqG);
		PSSM[j][4] = log(pMatrix[j][4]/freqC);
	}

	cout << endl;
	
	for(int j=0; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			cout << PSSM[j][i] << "\t";
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