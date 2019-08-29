#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main()
{
    string line1, line2, line3, line4;
    ifstream f( "/Users/kivanccorut/Desktop/courses_uga/fall_2019/BINF_8500/data/sample1k.fastq" );
    string *name = new string[1000]; 
    string *seq = new string[1000];
    string *pls = new string[1000];
    string *qual= new string[1000];
    int i = 0;
    while(getline(f,line1) && getline(f,line2) && getline(f,line3) && getline(f,line4)){
        name[i] = line1;
        seq[i] = line2;
        pls[i] = line3;
        qual[i] = line4;
        i++;
        for (int y = 0; y <= 4; y++) 
            cout<< name[y] << endl << seq[y] <<endl;
    };
    return 0;
} 