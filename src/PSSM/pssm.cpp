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

double getGC(vector <string> dnaSeq)
{
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

void getReverse(char seq[])
{
    while (*seq) {
        switch(*seq){
        	case 'A':
            *seq = 'T';
            break;

        	case 'G':
            *seq = 'C';
            break;

        	case 'C':
            *seq = 'G';
            break;

        	case 'T':
            *seq = 'A';
            break;  
        }
        ++seq;
    }
    return;
}


void freqMatrix(vector <string> motifSeq, vector <string> dnaSeq, int cutoff)
{
	motifSeq = motifSeq;
	dnaSeq = dnaSeq;
	cutoff = cutoff;

	double fMatrix[motifSeq[0].size()+1][5];

	// Frequency matrix
	for(int i=1; i<=motifSeq[0].size(); i++){
    	fMatrix[i][0] = i;
  	}
  	
	fMatrix[0][0] = 0;
	fMatrix[0][1] = 'A';
	fMatrix[0][2] = 'T';
	fMatrix[0][3] = 'G';
	fMatrix[0][4] = 'C';

	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=1; i<=4;i++){
			fMatrix[j][i] = 0;
		}
	}

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
	ofstream fout1("frequencyMatrix.txt");
	for(int j=0; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			fout1 << fMatrix[j][i] << "\t";
			if (i==4)
			{
				fout1<< "\n";
			}
		}
	}
	fout1.close();
	
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
	
	pMatrix[0][0] = 0;
	pMatrix[0][1] = 'A';
	pMatrix[0][2] = 'T';
	pMatrix[0][3] = 'G';
	pMatrix[0][4] = 'C';

	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=1; i<=4;i++){
			pMatrix[j][i] = fMatrix[j][i]/motifSeq.size();
		}
	}	
	
	ofstream fout2("probabilityMatrix.txt");
	for(int j=0; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			fout2 << pMatrix[j][i] << "\t";
			if (i==4)
			{
				fout2<< "\n";
			}
		}
	}
	fout2.close();
	
	double GCfreq = getGC(dnaSeq);
	double freqG, freqC, freqA, freqT;
	
	freqG=GCfreq/2;
	freqC=GCfreq/2;
	freqA=(1-GCfreq)/2;
	freqT=(1-GCfreq)/2;

	double PSSM[motifSeq[0].size()+1][5];
	
	for(int i=1; i<=motifSeq[0].size(); i++){
		PSSM[i][0] = i;
	}
	
	PSSM[0][0] = 0;
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
	
	ofstream fout3("PSSM.txt");
	for(int j=0; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			fout3 << PSSM[j][i] << "\t";
			if (i==4)
			{
				fout3<< "\n";
			}
		}
	}
	fout3.close();
	
	int numMotifs = motifSeq.size();
	int numNucMotif = motifSeq[0].size();

	double score[numMotifs];
	
	for(int i=0;i<numMotifs;i++){
		for(int j=0;j<numNucMotif;j++){
			if(toupper(motifSeq[i][j])=='A'){
				score[i]+=PSSM[j+1][1];
			}else if(toupper(motifSeq[i][j])=='T'){
				score[i]+=PSSM[j+1][2];
			}else if(toupper(motifSeq[i][j])=='G'){
				score[i]+=PSSM[j+1][3];
			}else if(toupper(motifSeq[i][j])=='C'){
				score[i]+=PSSM[j+1][4];
			}
		}
	}
	
	ofstream fout4("motifScore.txt");
	for(int i=0;i<numMotifs;i++){
		fout4 << "Motif Sequence: " << motifSeq[i] << ",	" << "Score: " << score[i] << endl;
	}
	fout4.close();

	int k=0;
	char dnaNuc[dnaSeq[0].size()*dnaSeq.size()];
	for(int i=0;i<dnaSeq.size();i++){
		for(int j=0;j<dnaSeq[0].size();j++){
			dnaNuc[k++] = dnaSeq[i][j];
		}
	}

	int n = (dnaSeq[0].size()*dnaSeq.size());
	int m = 16;

	vector<pair<string, double> > dnaMotif(n - m + 1);
	vector<pair<int, int> > positions(n - m + 1);

	for (int i = 0; i < n - m + 1; i++) {
        for (int j = 0; j < m; j++){
			if(toupper(dnaNuc[i+j])=='A'){
				dnaMotif[i].second+=PSSM[j+1][1];
				dnaMotif[i].first+=dnaNuc[i+j];
				positions[i].first=i;
				positions[i].second=i+j;

			}else if(toupper(dnaNuc[i+j])=='T'){
				//dnaScore[i]+=PSSM[j+1][2];
				dnaMotif[i].second+=PSSM[j+1][2];
				dnaMotif[i].first+=dnaNuc[i+j];
				positions[i].first=i;
				positions[i].second=i+j;

			}else if(toupper(dnaNuc[i+j])=='G'){
				dnaMotif[i].second+=PSSM[j+1][3];
				dnaMotif[i].first+=dnaNuc[i+j];
				positions[i].first=i;
				positions[i].second=i+j;

			}else if(toupper(dnaNuc[i+j])=='C'){
				dnaMotif[i].second+=PSSM[j+1][4];
				dnaMotif[i].first+=dnaNuc[i+j];
				positions[i].first=i;
				positions[i].second=i+j;
			}
			
		}
	}

	ofstream fout5("results.txt");
	for (int i=0; i<dnaMotif.size(); i++){ 
        if(dnaMotif[i].second>cutoff){
			fout5 << "Start: " << positions[i].first << ", " << "End: " << positions[i].second << ", " << "Sequence: " << dnaMotif[i].first << ", " << "Score: " << dnaMotif[i].second << endl;     
    	}	
	}
	fout5.close();
	
	getReverse(dnaNuc);

	int nn = (dnaSeq[0].size()*dnaSeq.size());
	int mm = 16;

	vector<pair<string, double> > dnaMotifRev(nn - mm + 1);
	vector<pair<int, int> > positionsRev(nn - mm + 1);

	for (int i = 0; i < nn - mm + 1; i++) {
        for (int j = 0; j < mm; j++){
			if(toupper(dnaNuc[i+j])=='A'){
				dnaMotifRev[i].second+=PSSM[j+1][1];
				dnaMotifRev[i].first+=dnaNuc[i+j];
				positionsRev[i].first=i;
				positionsRev[i].second=i+j;

			}else if(toupper(dnaNuc[i+j])=='T'){
				dnaMotifRev[i].second+=PSSM[j+1][2];
				dnaMotifRev[i].first+=dnaNuc[i+j];
				positionsRev[i].first=i;
				positionsRev[i].second=i+j;

			}else if(toupper(dnaNuc[i+j])=='G'){
				dnaMotifRev[i].second+=PSSM[j+1][3];
				dnaMotifRev[i].first+=dnaNuc[i+j];
				positionsRev[i].first=i;
				positionsRev[i].second=i+j;

			}else if(toupper(dnaNuc[i+j])=='C'){
				dnaMotifRev[i].second+=PSSM[j+1][4];
				dnaMotifRev[i].first+=dnaNuc[i+j];
				positionsRev[i].first=i;
				positionsRev[i].second=i+j;
			}
			
		}
	}

	ofstream fout6("resultsReverse.txt");
	for (int i=0; i<dnaMotifRev.size(); i++){ 
        if(dnaMotifRev[i].second>cutoff){
			fout6 << "Start: " << positionsRev[i].first << ", " << "End: " << positionsRev[i].second << ", " << "Sequence: " << dnaMotifRev[i].first << ", " << "Score: " << dnaMotifRev[i].second << endl;     
    	}	
	}
	fout6.close();
}

int main(int argc, char **argv) 
{	
	int cutoff;
  	
	cout<<"Enter the minimum score cutoff: ";
  	cin>> cutoff;

	// File Parsing.
	if (argc > 1) {
        cout << "Motif sequences file: " << argv[1] << endl;
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

	freqMatrix(motifSeq, dnaSeq, cutoff);

    return 0;
}