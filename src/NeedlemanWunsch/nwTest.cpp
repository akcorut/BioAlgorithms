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


void initMatrices(int ** sMatrix, int ** tbMatrix, string seqFirst, string seqSecond, int gap){
	
	seqFirst = seqFirst;
	seqSecond = seqSecond;

	int lengthFirst = seqFirst.length();
	int lengthSecond = seqSecond.length();
	
	sMatrix[0][0]=0;
	tbMatrix[0][0]=0;
  	
	for(int i=1; i<=lengthFirst; i++){
    	sMatrix[i][0] = sMatrix[i-1][0] + gap;
    	tbMatrix[i][0]=0;
  	}

  	for(int i=1; i<=lengthSecond; i++){
    	sMatrix[0][i] = sMatrix[0][i-1] + gap;
    	tbMatrix[0][i]=2;
  	}
}

void nwAlign(int ** sMatrix, int ** tbMatrix, string seqFirst, string seqSecond, string& alignmentFirst, string& alignmentSecond, int match, int mismatch, int gap)
{
	seqFirst = seqFirst;
	seqSecond = seqSecond;
	tbMatrix = tbMatrix;
	sMatrix = sMatrix;	
	int lengthFirst = seqFirst.length();
	int lengthSecond = seqSecond.length();

	int TBx[]={-1,-1,0};
  	int TBy[]={0,-1,-1};

  	for(int i=1; i<=lengthFirst; i++){
		  
    	for(int j=1; j<=lengthSecond; j++){
			int diagScore = 0;
      		int verScore = sMatrix[i][j-1] + gap;
      		int horScore = sMatrix[i-1][j] + gap;
      		
			if(seqFirst.at(i-1) == seqSecond.at(j-1))
			{
        		diagScore = sMatrix[i-1][j-1] + match;
      		
			}else
			{
        		diagScore = sMatrix[i-1][j-1] + mismatch;
      		}
      		
			if(diagScore >= verScore && diagScore >= horScore)
			{
        		sMatrix[i][j] = diagScore;
        		tbMatrix[i][j]=1;
      		
			}else if(verScore >= diagScore && verScore >= horScore)
			{
        		sMatrix[i][j] = verScore;
        		tbMatrix[i][j]=2;
      		
			}else
			{
        		sMatrix[i][j] = horScore;
        		tbMatrix[i][j]=0;
      		}
    	}	
  	}

	while(!(lengthFirst==0&&lengthSecond==0)){
		switch (tbMatrix[lengthFirst][lengthSecond])
		{
			case 2:      	alignmentFirst += "-";
                        	alignmentSecond += seqSecond.at(lengthSecond-1);;
                        	break;

        	case 1:      	alignmentFirst += seqFirst.at(lengthFirst-1);
                        	alignmentSecond += seqSecond.at(lengthSecond-1);
                        	break;

        	case 0:      	alignmentFirst += seqFirst.at(lengthFirst-1);
							alignmentSecond += "-";
		}
		int x=lengthFirst+TBx[tbMatrix[lengthFirst][lengthSecond]];
    	int y=lengthSecond+TBy[tbMatrix[lengthFirst][lengthSecond]];
    	lengthFirst=x;
    	lengthSecond=y;
	}
	reverse(alignmentFirst.begin(), alignmentFirst.end());
    reverse(alignmentSecond.begin(), alignmentSecond.end());
}

int needlemanWunsch(string seqFirst, string seqSecond, string& alignmentFirst, string& alignmentSecond, int match, int mismatch, int gap){

	seqFirst = seqFirst;
	seqSecond = seqSecond;
	int lengthFirst = seqFirst.length();
    int lengthSecond = seqSecond.length();

    int **sMatrix = new int *[lengthSecond+1];
    for(int i = 0; i <= lengthSecond; i++)  
        sMatrix[i] = new int [lengthFirst];
 
    int **tbMatrix = new int *[lengthSecond+1];
    for(int i = 0; i <= lengthSecond; i++) 
        tbMatrix[i] = new int [lengthFirst];
        
    initMatrices(sMatrix, tbMatrix, seqFirst, seqSecond, gap);
    nwAlign(sMatrix, tbMatrix, seqFirst, seqSecond, alignmentFirst, alignmentSecond, match, mismatch, gap);

    return 0;
}

void printAlignment(string& alignmentFirst, string& alignmentSecond){
        ofstream fout("outputAlignment.txt");
        fout << alignmentFirst << endl;
        fout << alignmentSecond << endl;
        fout.close();
}

int main(int argc, char **argv) 
{	
	int match,mismatch,gap;
  	cout<<"Enter the score for matches: ";
  	cin>>match;
  	cout<<"Enter the score for mismatches: ";
  	cin>>mismatch;
  	cout<<"Enter the gap penalty: ";
  	cin>>gap;

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
    
	//cout << seqFirst << endl;
    //cout << seqFirst.length() << endl;
    //cout << seqSecond << endl;
    //cout << seqSecond.length() << endl;
	
	string alignmentFirst;
    string alignmentSecond;
    
    needlemanWunsch(seqFirst, seqSecond, alignmentFirst, alignmentSecond, match, mismatch, gap);
    printAlignment(alignmentFirst, alignmentSecond);
	return 0;
}