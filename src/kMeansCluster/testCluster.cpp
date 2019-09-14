#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{  
    if (argc > 1) {
        cout << "k-means clustering input file: " << argv[1] << endl; 
    } else {
        cout << "No file detected. Please provide a proper path for the input.";
        return -1;
    }
    ifstream infile(argv[1]);

    int colCount = 20;
    int rowCount= 92;

    double** Point = new double*[rowCount];
    for(int i = 0; i < rowCount; ++i)
        Point[i] = new double[colCount];

    int len=92;
    double Aa[len];
    double Rr[len];
    double Nn[len];
    double Dd[len];
    double Cc[len];	
    double Qq[len];	
    double Ee[len];
    double Gg[len];
    double Hh[len];
    double Ii[len];
    double Ll[len];
    double Kk[len];
    double Mm[len];
    double Ff[len];
    double Pp[len];
    double Ss[len];
    double Tt[len];
    double Ww[len];
    double Yy[len];
    double Vv[len];
    string line;
    //double am[92];
    while(getline(infile, line)) {
        stringstream linestream(line);
        string genome;
        getline(linestream, genome, '\t');
        int kk;
        linestream >> Aa[kk];
        linestream >> Rr[kk];
        linestream >> Nn[kk];
        linestream >> Dd[kk];
        linestream >> Cc[kk];	
        linestream >> Qq[kk];	
        linestream >> Ee[kk];
        linestream >> Gg[kk];
        linestream >> Hh[kk];
        linestream >> Ii[kk];
        linestream >> Ll[kk];
        linestream >> Kk[kk];
        linestream >> Mm[kk];
        linestream >> Ff[kk];
        linestream >> Pp[kk];
        linestream >> Ss[kk];
        linestream >> Tt[kk];
        linestream >> Ww[kk];
        linestream >> Yy[kk];
        linestream >> Vv[kk];
        //linestream >> aas[j][0];
        //linestream >> aas[j][1];
        ++kk;
            //linestream>> amino.Nn >> amino.Dd >> amino.Cc >> amino.Qq 
            //>> amino.Ee >> amino.Gg >> amino.Hh >> amino.Ii >> amino.Ll >> amino.Kk >> amino.Mm 
            //>> amino.Ff >> amino.Pp >> amino.Ss >> amino.Tt >> amino.Ww >> amino.Yy >> amino.Vv;
    }
    int index;
    for (int k=0; k < rowCount; k++) {
        Point[k][1] = Aa[index];
        Point[k][2] = Rr[index];
        Point[k][3] = Nn[index];
        Point[k][4] = Dd[index];
        Point[k][5] = Cc[index]; 
        Point[k][6] = Qq[index];
        Point[k][7] = Ee[index];
        Point[k][8] = Gg[index];
        Point[k][9] = Hh[index];
        Point[k][10] = Ii[index];
        Point[k][11] = Ll[index];
        Point[k][12] = Kk[index];
        Point[k][13] = Mm[index];
        Point[k][14] = Ff[index];
        Point[k][15] = Pp[index];
        Point[k][16] = Ss[index];
        Point[k][17] = Tt[index];
        Point[k][18] = Ww[index]; 
        Point[k][19] = Yy[index];
        Point[k][20] = Vv[index];
        ++index;
    }
    //cout << aas[0][6];
    /*for (int y = 0; y < rowCount; ++y) {
        for (int z = 0; z < colCount; ++z) {
            cout << Point[y][z] << "\t";
        }
    }*/
    cout << Point[91][1];
    return 0;
}