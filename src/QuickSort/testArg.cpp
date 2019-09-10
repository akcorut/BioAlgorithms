#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{   
    if (argc > 1) {
        cout << "File = " << argv[1] << endl; 
    } else {
        cout << "No file name entered. Exiting...";
        return -1;
    }
    ifstream infile(argv[1]);

    string line1, line2, line3, line4;
    string *id = new string[100000]; 
    string *seq = new string[100000];
    string *pls = new string[100000];
    string *qual= new string[100000];
    int i = 0;
    while(getline(infile,line1) && getline(infile,line2) && getline(infile,line3) && getline(infile,line4)){
        id[i] = line1;
        seq[i] = line2;
        pls[i] = line3;
        qual[i] = line4;
        i++;
    };
    //cout << id[1] << endl;
    //return 0;

    ofstream outfile;
    outfile.open("sorted.fastq");
    outfile << id[1];
    outfile.close();
    return 0;
}