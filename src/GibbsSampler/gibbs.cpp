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



double getGC(vector <string> dnaSeq)
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

void getRandomMotif(vector <string>& dnaSeq, vector <string>& motifSeq, int motifSize)
{
	random_device ran_dev;
	mt19937 gen(ran_dev());
	
	int randLimit = dnaSeq[0].size() - motifSize;
    uniform_int_distribution<> dist1(0, randLimit);

    for(int k=0; k < dnaSeq.size(); k++){
		int i = dist1(gen);
		cout << i << endl;
        motifSeq.push_back(dnaSeq[k].substr(i,motifSize));
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

	double GCfreq = getGC(dnaSeq);
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

void getMotifScore(vector <string>& motifSeq, int motifSize, double PSSM[][5], double& currentScore)
{
	int numMotifs = motifSeq.size();

	double score[numMotifs];
	
	for(int i=0;i<numMotifs;i++){
		for(int j=0;j<motifSize;j++){
			if(motifSeq[i][j]=='A' || motifSeq[i][j]=='a'){
				score[i]+=PSSM[j+1][1];
			}else if(motifSeq[i][j]=='T' || motifSeq[i][j]=='t'){
				score[i]+=PSSM[j+1][2];
			}else if(motifSeq[i][j]=='G' || motifSeq[i][j]=='g'){
				score[i]+=PSSM[j+1][3];
			}else if(motifSeq[i][j]=='C' || motifSeq[i][j]=='c'){
				score[i]+=PSSM[j+1][4];
			}
		}
	}
	
	ofstream fout5("motifScore.txt");
	for(int i=0;i<numMotifs;i++){
		fout5 << "Motif Sequence: " << motifSeq[i] << ",	" << "Score: " << score[i] << "\n";
	}
	fout5.close();

	for(int i=0;i<numMotifs;i++){
		currentScore += score[i];
	}
}

void findMotif(vector<string>& dnaSeq, vector<string>& motifSeq, int motifSize, int targetIndex, double PSSM[][5], vector<pair<int, int> >& newPositions)
{	
	random_device ran_dev;
	mt19937 gen(ran_dev());
	
	char *targetSeq = new char[dnaSeq[targetIndex].size()];
	
	int k=0;
	for(int j=0;j<dnaSeq[targetIndex].size();j++){
			targetSeq[k++] = dnaSeq[targetIndex][j];
	}
	
	int n = dnaSeq[targetIndex].size();
	int m = motifSize;

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
	/*for (int i=0; i<dnaMotif.size(); i++){
		totalScore +=  exp(dnaMotif[i].second);
	}
	cout<< totalScore;*/

	vector<double> expScore;

	for (int i=0; i<dnaMotif.size(); i++){
		expScore.push_back(exp(dnaMotif[i].second));
		totalScore += expScore[i];
	}
	//cout<< totalScore;


	/*for (int i=0; i<dnaMotif.size(); i++){
		cout << "Start: " << positions[i].first << ", " << "End: " << positions[i].second << ", " << "Sequence: " << dnaMotif[i].first << ", " << "Score: " << dnaMotif[i].second << "\n";     
    	
	}*/

	
	uniform_real_distribution<> dist2(0.0, 1.0);
	double randNum = dist2(gen);
	//cout << randNum << endl;
	
	vector<double> probs;
	double tempSum = 0.0;
	for (int i=0; i<expScore.size(); i++){
		probs.push_back(expScore[i]/totalScore+tempSum);
      	tempSum=tempSum+(expScore[i]/totalScore);
	}
	
	/*for (int i=0; i<probs.size(); i++){
		cout << probs[i] << "\t";
	}
	cout << endl;*/
	//fout4.close();
	
	int window = 0;
	while (probs[window] < randNum){
		window++;
	}
	//cout << "\t" << positions[window].first << "\t" << dnaMotif[window].first << endl;
	//cout << window << endl;
	//vector<double> tempScore;
	motifSeq[targetIndex] = dnaMotif[window].first;
	//tempScore.push_back(dnaMotif[window].second);

	newPositions[targetIndex].first = positions[window].first+1;
	newPositions[targetIndex].second = positions[window].second+1;

	expScore.clear();
	probs.clear();
	//tempScore.clear();
	/*cout << endl;
	for (string x:motifSeq){ 
        cout << x << endl;
    }*/
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

	vector<string> motifSeq;
	getRandomMotif(dnaSeq, motifSeq, motifSize);
	
	/*for (string x:motifSeq){ 
        cout << x << endl;
    }*/
	ofstream fout3("finalResults.txt");
	double PSSM[motifSeq[0].size()+1][5];
	vector<double> tempScore;
	//double currentScore;
	
	double finalScore =0.0;
	vector<string> finalMotif;

	vector<pair<int, int> > newPositions(motifSeq.size(), make_pair(0, 0));
	vector<pair<int, int> > lastPositions;

	for (int j=0; j<=500; j++){
		double currentScore =0.0;
		for (int i=0; i<motifSeq.size(); i++){
			getPSSM(dnaSeq, motifSeq, motifSize, i, PSSM);
			findMotif(dnaSeq, motifSeq, motifSize, i, PSSM, newPositions);
			//currentScore += tempScore[i];
			//fout3 << currentScore << endl;
		}
		double currentPSSM[motifSeq[0].size()+1][5];
		getPSSM(dnaSeq, motifSeq, motifSize, motifSeq.size()+10, currentPSSM);
		getMotifScore(motifSeq, motifSize, currentPSSM, currentScore);
		//fout3 << currentScore << endl;
		if (currentScore > finalScore)
		{
			finalScore = currentScore;
			finalMotif = motifSeq;
			lastPositions = newPositions;
		}
		fout3 << finalScore << endl;
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