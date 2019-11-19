/* 
Author: Adnan Kivanc Corut
Usage:
This program takes one input fasta file and returns
2 output files including matrices and final motifs with score and locations. 
User is also expected to define the motif length from command line.

./gibbs <input fasta>

Example:
./gibbs ../../data/GibbsSampler/H.pyloriRpoN-sequences-10-300nt.fasta
*/

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
#include <random>

using namespace std;

double getGC(vector <string> dnaSeq) // Function to get GC freq of given dna sequence
{

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

// Function to randomly select motifs with given size
void getRandomMotif(vector <string>& dnaSeq, vector <string>& motifSeq, int motifSize, vector <int>& randPos)
{
	random_device ran_dev;
	mt19937 gen(ran_dev());
	
	int randLimit = dnaSeq[0].size() - motifSize;
    uniform_int_distribution<> dist1(0, randLimit);

    for(int k=0; k < dnaSeq.size(); k++){
		int i = dist1(gen);
		cout << i << endl;
		randPos.push_back(i);
        motifSeq.push_back(dnaSeq[k].substr(i,motifSize));
    }
}

// Function to change motif size and/or start position of th motifs
void changeMotifSize(vector <string>& dnaSeq, vector <string>& newMotifs, int motifSize, vector <int>& startPos, int i, int j)
{
    
	for(int k=0; k < dnaSeq.size(); k++){
        newMotifs.push_back(dnaSeq[k].substr(startPos[k]+i,motifSize+j));
    }
	
}

void getPSSM(vector<string>& dnaSeq, vector<string>& motifSeq, int motifSize, int skipSeq, double PSSM[][5])
{
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

	ofstream fout2("matrices.txt");
	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			fout2 << fMatrix[j][i] << "\t";
			if (i==4)
			{
				fout2<< "\n";
			}
		}
	}
	fout2<< "\n";

	for(int j=0; j<motifSeq.size();j++){
		for (int i=0; i<motifSeq[0].size();i++){
			if(j !=skipSeq ){
				if(motifSeq[j][i]=='A' || motifSeq[j][i]=='a'){
					fMatrix[i+1][1]+=1;
				}else if(motifSeq[j][i]=='T' || motifSeq[j][i]=='t'){
					fMatrix[i+1][2]+=1;
				}else if(motifSeq[j][i]=='G' || motifSeq[j][i]=='g'){
					fMatrix[i+1][3]+=1;
				}else if(motifSeq[j][i]=='C' || motifSeq[j][i]=='c'){
					fMatrix[i+1][4]+=1;	
				}
			}
		}
	}
	string header[5] = {" ", "A", "T", "G", "C"};

	
	for (int i=0; i<=4;i++){
			fout2 << header[i] << "\t";
	}
	fout2 << "\n";
	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			fout2 << fMatrix[j][i] << "\t";
			if (i==4)
			{
				fout2<< "\n";
			}
		}
	}
	fout2<< "\n";

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
	
	for (int i=0; i<=4;i++){
			fout2 << header[i] << "\t";
	}
	fout2 << "\n";
	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			fout2 << pMatrix[j][i] << "\t";
			if (i==4)
			{
				fout2<< "\n";
			}
		}
	}
	fout2<< "\n";

	vector <string> seqToFreq;
	seqToFreq.push_back(dnaSeq[skipSeq]);

	double GCfreq = getGC(seqToFreq);
	double freqG, freqC, freqA, freqT;
	
	freqG=GCfreq/2;
	freqC=GCfreq/2;
	freqA=(1-GCfreq)/2;
	freqT=(1-GCfreq)/2;

	//double PSSM[motifSeq[0].size()+1][5];
	
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
	
	for (int i=0; i<=4;i++){
			fout2 << header[i] << "\t";
	}
	fout2 << "\n";
	for(int j=1; j<=motifSeq[0].size();j++){
		for (int i=0; i<=4;i++){
			fout2 << PSSM[j][i] << "\t";
			if (i==4)
			{
				fout2<< "\n";
			}
		}
	}
	fout2<< "\n";

}

void getMotifScore(vector <string>& motifSeq, int motifSize, double PSSM[][5], double& motifScore, int targetSeq)
{
	//int numMotifs = motifSeq.size();

	double score;
	
	
	for(int j=0;j<motifSeq[0].size();j++){
		if(motifSeq[targetSeq][j]=='A' || motifSeq[targetSeq][j]=='a'){
			score+=PSSM[j+1][1];
		}else if(motifSeq[targetSeq][j]=='T' || motifSeq[targetSeq][j]=='t'){
			score+=PSSM[j+1][2];
		}else if(motifSeq[targetSeq][j]=='G' || motifSeq[targetSeq][j]=='g'){
			score+=PSSM[j+1][3];
		}else if(motifSeq[targetSeq][j]=='C' || motifSeq[targetSeq][j]=='c'){
			score+=PSSM[j+1][4];
		}
	}
	
	
	ofstream fout5("motifScore.txt");
	
	fout5 << "Motif Sequence: " << motifSeq[targetSeq] << ",	" << "Score: " << score << "\n";
	
	fout5.close();

	
	motifScore = score;
	
}

void findMotif(vector<string>& dnaSeq, vector<string>& motifSeq, int motifSize, int targetIndex, double PSSM[][5], vector<pair<int, int> >& newPositions, int shift, vector<int>& startPos)
{	
	random_device ran_dev;
	mt19937 gen(ran_dev());
	
	char *targetSeq = new char[dnaSeq[targetIndex].size()];
	
	int k=0;
	for(int j=0;j<dnaSeq[targetIndex].size();j++){
			targetSeq[k++] = dnaSeq[targetIndex][j];
	}
	
	int n = dnaSeq[targetIndex].size();
	int m = motifSeq[0].size();

	vector<pair<string, double> > dnaMotif(n - m + 1);
	vector<pair<int, int> > positions(n - m + 1);

	for (int i = 0; i < n - m + 1; i++) {
        for (int j = 0; j < m; j++){
			if(targetSeq[i+j]=='A' || targetSeq[i+j]=='a'){
				dnaMotif[i].second+=PSSM[j+1][1];
				dnaMotif[i].first+=targetSeq[i+j];
				positions[i].first=i;
				positions[i].second=i+j;

			}else if(targetSeq[i+j]=='T' || targetSeq[i+j]=='t'){
				dnaMotif[i].second+=PSSM[j+1][2];
				dnaMotif[i].first+=targetSeq[i+j];
				positions[i].first=i;
				positions[i].second=i+j;

			}else if(targetSeq[i+j]=='G' || targetSeq[i+j]=='g'){
				dnaMotif[i].second+=PSSM[j+1][3];
				dnaMotif[i].first+=targetSeq[i+j];
				positions[i].first=i;
				positions[i].second=i+j;

			}else if(targetSeq[i+j]=='C' || targetSeq[i+j]=='c'){
				dnaMotif[i].second+=PSSM[j+1][4];
				dnaMotif[i].first+=targetSeq[i+j];
				positions[i].first=i;
				positions[i].second=i+j;
			}
			
		}
	}

	ofstream fout("motifs.txt");
	for (int i=0; i<dnaMotif.size(); i++){
		fout << "Start: " << positions[i].first << ", " << "End: " << positions[i].second << ", " << "Sequence: " << dnaMotif[i].first << ", " << "Score: " << dnaMotif[i].second << "\n";     
    	
	}

	double totalScore = 0.0;
	

	vector<double> expScore;

	for (int i=0; i<dnaMotif.size(); i++){
		expScore.push_back(exp2(dnaMotif[i].second));
		totalScore += expScore[i];
	}
	

	
	uniform_real_distribution<> dist2(0.0, 1.0);
	double randNum = dist2(gen);
	//cout << randNum << endl;
	
	vector<double> probs;
	double tempSum = 0.0;
	for (int i=0; i<expScore.size(); i++){
		probs.push_back(expScore[i]/totalScore+tempSum);
      	tempSum=tempSum+(expScore[i]/totalScore);
	}
	
	
	int window = 0;
	while (probs[window] < randNum){
		window++;
	}
	
	bool outOfindex = false;
	if (window+shift < 0 || window+shift >= n - m + 1)
		outOfindex = true; 
	else 
		outOfindex = false;  
	if(!outOfindex){
		motifSeq[targetIndex] = dnaMotif[window+shift].first;
	
		newPositions[targetIndex].first = positions[window+shift].first+1;
		newPositions[targetIndex].second = positions[window+shift].second+1;
	}
	
	startPos[targetIndex] = positions[window+shift].first;
	
	expScore.clear();
	probs.clear();
	
}

int main(int argc, char **argv) 
{	
	int motifSize; // Motif size
  	
	cout<<"Enter the motif size: ";
  	cin>> motifSize;

	// File Parsing.
	if (argc > 1) {
        cout << "Input fasta file: " << argv[1] << endl;

    } else {
        cout << "No file detected. Please provide the proper path of the input fasta file.";
        return -1;
    }

	vector<string> dnaSeq; // Declare a string vector for DNA sequences
    vector<string> header; // Declare another string vector for headers
    string line;
    string seq;

    ifstream fastaDNA(argv[1]); // Read the file from command line
    while (getline(fastaDNA, line)){
        if(line[0] == '>'){
            header.push_back(line);
            if(!seq.empty()){
                dnaSeq.push_back(seq);
                seq.clear();
            }
        }else{
            seq = seq + line;
        }
    }

    dnaSeq.push_back(seq);
    for (int i=0; i<dnaSeq.size(); i++){ 
        cout << "ID: " << header[i] << "\t" << "Sequence: " << dnaSeq[i]<< endl;
    }
	cout << endl;
    /*for(string x:dnaSeq){
        cout << x << endl;
    }*/

	for (auto &i : dnaSeq)
    	for (auto &j : i)      
        	j = toupper(j);

	vector <int> randPos;
	vector<string> motifSeq;

	getRandomMotif(dnaSeq, motifSeq, motifSize, randPos);
	
	/*vector<string> testM;
	changeMotifSize(dnaSeq, testM, motifSize, randPos, 0, 2);
	for(string h:testM){
		cout<<h<<endl;
	}*/
	
	ofstream fout3("finalResults.txt");
	double PSSM[motifSeq[0].size()+1][5];
	//vector<double> tempScore;
	
	double finalScore =0.0;
	vector<string> finalMotif;
	//vector<string> sMotif;
	//vector<string> bMotif;

	vector<int> startPos(motifSeq.size(), 0);

	vector<pair<int, int> > newPositions(motifSeq.size(), make_pair(0, 0));
	vector<pair<int, int> > lastPositions;
	double motifScore;
	int iter = 0;
	
	while(iter < 300){
		for (int j=0; j<=300; j++){
			double currentScore =0.0;
			//double tmpScore =0.0;
			//double tmpScore2 =0.0;
			for (int i=0; i<motifSeq.size(); i++){
				getPSSM(dnaSeq, motifSeq, motifSize, i, PSSM);
				findMotif(dnaSeq, motifSeq, motifSize, i, PSSM, newPositions, 0, startPos);
				//cout << motifSeq[i] << endl;
				//currentScore += tempScore[i];
				//fout3 << currentScore << endl;
			}
			if (j % 5==0) {
				for (int i=0; i<motifSeq.size(); i++){
					getPSSM(dnaSeq, motifSeq, motifSize, i, PSSM);
					findMotif(dnaSeq, motifSeq, motifSize, i, PSSM, newPositions, 1, startPos);
				}
			}
			if (j % 12==0) {
				for (int i=0; i<motifSeq.size(); i++){
					getPSSM(dnaSeq, motifSeq, motifSize, i, PSSM);
					findMotif(dnaSeq, motifSeq, motifSize, i, PSSM, newPositions, -1, startPos);
				}
			}
			/*if (j % 7==0) {
				changeMotifSize(dnaSeq, sMotif, motifSize, startPos, 0, -2);
				double tempPSSM[sMotif[0].size()+1][5];
				for (int k=0; k<motifSeq.size(); k++){
					getPSSM(dnaSeq, sMotif, motifSize, k, tempPSSM);
					getMotifScore(sMotif, motifSize, tempPSSM, motifScore, k);
					//cout << motifScore << endl;
					tmpScore += motifScore;
				}
			}*/
			/*if (j % 13==0) {
				changeMotifSize(dnaSeq, bMotif, motifSize, startPos, 0, 2);
				double tempPSSM[bMotif[0].size()+1][5];
				for (int k=0; k<motifSeq.size(); k++){
					getPSSM(dnaSeq, bMotif, motifSize, k, tempPSSM);
					getMotifScore(bMotif, motifSize, tempPSSM, motifScore, k);
					//cout << motifScore << endl;
					tmpScore2 += motifScore;
				}
			}*/
			

			double currentPSSM[motifSeq[0].size()+1][5];
			for (int k=0; k<motifSeq.size(); k++){
				getPSSM(dnaSeq, motifSeq, motifSize, k, currentPSSM);
				getMotifScore(motifSeq, motifSize, currentPSSM, motifScore, k);
				//cout << motifScore << endl;
				currentScore += motifScore;
			}
			//getMotifScore(motifSeq, motifSize, currentPSSM, currentScore);
			//fout3 << currentScore << endl;
			
			if(currentScore <= finalScore){
				iter++;
			}
			
			else if(currentScore > finalScore)
			{
				finalScore = currentScore;
				finalMotif = motifSeq;
				lastPositions = newPositions;
				iter=0;
			}
			fout3 << finalScore << endl;

			/*ofstream fout6("alternativeMotifs.txt");
			if(tmpScore > finalScore){
				fout6 << tmpScore << endl;
				for(string y:sMotif){
					fout6<<y<<endl;
				}
			}*/

			/*if(tmpScore2 > finalScore){
				fout6 << tmpScore2 << endl;
				for(string z:bMotif){
					fout6<<z<<endl;
				}
			}*/
		}
	}
	/*for(string y:finalMotif){
		fout3<<y<<endl;
	}
	fout3 << finalScore << endl;*/

	for (int i=0; i<finalMotif.size(); i++){
        fout3 << "Start: " << lastPositions[i].first << "\t" << "End: " << lastPositions[i].second<< "\t" << "Motif: " << finalMotif[i]<< endl;
    }
	fout3 << "Final Score: " << finalScore << endl;
    return 0;
}