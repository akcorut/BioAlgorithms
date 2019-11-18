/* 
Usage:
This program takes two files in fasta format as inputs and returns
an alignment file as output. User is also expected to define the 
required scores for mismatch, match and gap penalty from command line
by typing them.

./nwAlign <input fasta sequence 1> <input fasta sequence 2>

Example:
./nwAlign /path/to/RpoB-B.subtilis.fasta /path/to/RpoB-E.coli.fasta
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

using namespace std;

void reverseAlignment(string& alignment) 
{ 
	int lastIndex = alignment.length()-1;
    for (int i=0; i<alignment.length()/2; i++){
		swap(alignment[i], alignment[lastIndex-i]);
	}
}

void initMatrices(int ** sMatrix, int ** tbMatrix, string seqFirst, string seqSecond, int gap)
{	

	int lengthFirst = seqFirst.length(); // Get length of the first sequence
	int lengthSecond = seqSecond.length(); // Get length of the second sequence
	
	sMatrix[0][0]=0; // Scoring matrix
	tbMatrix[0][0]=0; // Traceback matrix
  	
	// Initialize the matrices with first row and column
	for(int i=1; i<=lengthFirst; i++){
    	sMatrix[i][0] = sMatrix[i-1][0] + gap;
    	tbMatrix[i][0]=0;
  	}

  	for(int i=1; i<=lengthSecond; i++){
    	sMatrix[0][i] = sMatrix[0][i-1] + gap;
    	tbMatrix[0][i]=2;
  	}

	ofstream fout3("InitScoreMatrix.txt");
	for(int i = 0; i <= lengthFirst; i++){
      	for(int j = 0; j < lengthSecond; j++) {
        	fout3 << sMatrix[i][j] << '\t';
    	}
    	fout3 << sMatrix[i][lengthSecond] << '\n';
    }

	ofstream fout5("inittbmatrix.csv");
	for(int i = 0; i <= lengthFirst; i++){
      	for(int j = 0; j < lengthSecond; j++) {
        	fout5 << tbMatrix[i][j] << ',';
    	}
    	fout5 << tbMatrix[i][lengthSecond] << '\n';
    }
}

void nwAlign(int ** sMatrix, int ** tbMatrix, string& seqFirst, string& seqSecond, string& alignmentFirst, string& alignmentSecond, int match, int mismatch, int gap)
{	
	int diagScore, verScore, horScore;
	int lengthFirst = seqFirst.length(); // Get length of the first sequence
	int lengthSecond = seqSecond.length(); // Get length of the second sequence

  	for(int i=1; i<=lengthFirst; i++){
		  
    	for(int j=1; j<=lengthSecond; j++){
			int diagScore = 0;
      		int verScore = sMatrix[i][j-1] + gap; // Deletion cost
      		int horScore = sMatrix[i-1][j] + gap; // Insertion cost
      		
			if(seqFirst.at(i-1) == seqSecond.at(j-1))
			{
        		diagScore = sMatrix[i-1][j-1] + match; // If match
      		
			}else
			{
        		diagScore = sMatrix[i-1][j-1] + mismatch; // If mismatch
      		}
      		
			// Calculate maximum alignment score / best path
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

	ofstream fout2("scorematrix.csv");
	for(int i = 0; i <= lengthFirst; i++){
      	for(int j = 0; j < lengthSecond; j++) {
        	fout2 << sMatrix[i][j] << ',';
    	}
    	fout2 << sMatrix[i][lengthSecond] << '\n';
    }

	ofstream fout4("tbmatrix.csv");
	for(int i = 0; i <= lengthFirst; i++){
      	for(int j = 0; j < lengthSecond; j++) {
        	fout4 << tbMatrix[i][j] << ',';
    	}
    	fout4 << tbMatrix[i][lengthSecond] << '\n';
    }
	
	int k=0;
	while(!(lengthFirst == 0 && lengthSecond == 0)){
		switch (tbMatrix[lengthFirst][lengthSecond])
		{
			// If the second seqeunce aligned with a gap
			case 2:     alignmentFirst += "-";
                        alignmentSecond += seqSecond[lengthSecond-1];
						lengthSecond--;
                        break;

			// If both seqeunces are aligned
        	case 1:     alignmentFirst += seqFirst[lengthFirst-1];
                        alignmentSecond += seqSecond[lengthSecond-1];
						lengthFirst--; lengthSecond--;
                        break;

			// If the first seqeunce aligned with a gap
        	case 0:     alignmentSecond += "-";
						alignmentFirst += seqFirst[lengthFirst-1];			
						lengthFirst--;
		}
		k++;
	}

	// Trace aligned sequences by reversing them
	reverseAlignment(alignmentFirst);
    reverseAlignment(alignmentSecond);
}

int needlemanWunsch(string seqFirst, string seqSecond, string& alignmentFirst, string& alignmentSecond, int match, int mismatch, int gap)
{
	
	
	int lengthFirst = seqFirst.length(); // Get length of the first sequence
    int lengthSecond = seqSecond.length(); // Get length of the second sequence

	// Dynamically allocate the scoring matrix
	int **sMatrix = new int *[lengthFirst+1];
    for(int i = 0; i <= lengthFirst; i++)  
        sMatrix[i] = new int [lengthSecond+1];
	
	// Dynamically allocate the traceback matrix
    int **tbMatrix = new int *[lengthFirst+1];
    for(int i = 0; i <= lengthFirst; i++) 
        tbMatrix[i] = new int [lengthSecond+1];
    
	initMatrices(sMatrix, tbMatrix, seqFirst, seqSecond, gap); // Call the function for initializing the matrices
    nwAlign(sMatrix, tbMatrix, seqFirst, seqSecond, alignmentFirst, alignmentSecond, match, mismatch, gap); // Call the alignment function
	
	return 0;
}

void printAlignment(string& alignmentFirst, string& alignmentSecond)
{
		char *signs = new char[alignmentFirst.length()];
		// Generate the signs for the alignment output
		for(int i=0; i<=alignmentFirst.length(); i++){
			if(toupper(alignmentFirst[i]) == toupper(alignmentSecond[i])){
				signs[i]='|'; // If match
			}else if(alignmentFirst[i] == '-' || alignmentSecond[i] == '-'){
				signs[i]=' '; // If gap
			}else{
				signs[i]='*'; // If mismatch
			}
		}
		
		ofstream fout("outAlignment.txt"); // Write the final alignment to an output file
		for (int i = 0; i < alignmentFirst.length(); i++)
		{
			fout << alignmentFirst[i]; // Print out the first alignment in the first line
		}
		fout << "\n";

		for (int i = 0; i < alignmentFirst.length(); i++)
		{
			fout << signs[i]; // Print out the signs in the second line 
		}
		fout << "\n";

		for (int i = 0; i < alignmentSecond.length(); i++)
		{
			fout << alignmentSecond[i]; // Print out the second alignment in the third line 
		}
		fout << "\n";
		fout.close();
}

int main(int argc, char **argv) 
{	
	int match,mismatch,gap;
  	
	/*cout<<"Enter the score for matches: ";
  	cin>>match;
  	
	cout<<"Enter the score for mismatches: ";
  	cin>>mismatch;
  	
	cout<<"Enter the gap penalty: ";
  	cin>>gap;*/

	// File Parsing.
	if (argc > 1) {
        cout << "First fasta file: " << argv[1] << endl;
        cout << "Second fasta file: " << argv[2] << endl; 

    } else {
        cout << "No file detected. Please provide the proper paths of the input files.";
        return -1;
    }

	string seqFirst, seqSecond, line; // Declare the required string objects for reading the files

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
    
	match = atoi(argv[3]);  
    mismatch = atoi(argv[4]); 
    gap = atoi(argv[5]); 

    cout << "Length of the first sequence: " << seqFirst.length() << endl;
    cout << "Length of the second sequence: " << seqSecond.length() << endl;
	
	// Aligned sequences
	string alignmentFirst;
    string alignmentSecond;
    
	for (int i=0; i < seqFirst.length(); i++)
		seqFirst[i] = toupper(seqFirst[i]);
	
	for (int i=0; i < seqSecond.length(); i++)
		seqSecond[i] = toupper(seqSecond[i]);

    needlemanWunsch(seqFirst, seqSecond, alignmentFirst, alignmentSecond, match, mismatch, gap);
    printAlignment(alignmentFirst, alignmentSecond);
	
	cout << gap;
	return 0;
}