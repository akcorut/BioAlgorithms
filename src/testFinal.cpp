#include <fstream>
#include <iostream>
#include <string>

using namespace std;

void qSort (string a[], string b[], string c[], string d[], int left, int right){
    int j=left; 
    int k=right;
    string tmp;
    string pivot= a[(left+right)/2];
    while(j<=k){
        while((a[j]).compare(pivot) < 0)
            j++;
        while((a[k]).compare(pivot) > 0)
            k--;
        if(j<=k){
            tmp=a[j];
            a[j]=a[k];
            a[k]=tmp;
            tmp=b[j];
            b[j]=b[k];
            b[k]=tmp;
            tmp=c[j];
            c[j]=c[k];
            c[k]=tmp;
            tmp=d[j];
            d[j]=d[k];
            d[k]=tmp;
            j++;
            k--;
        }
    }
    if(left<k){
        qSort(a, b, c, d, left, k);
    }
    if(j<right){
        qSort(a, b, c, d, j, right);
    }
}

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

    qSort(seq, id, pls, qual, 0, 999);
    ofstream outfile;
    outfile.open("sorted.fastq");
    for (int y = 0; y < 10000; y++)
    {
        outfile << id[y] << endl << seq[y] << pls[y] << endl << qual[y] << endl;
    }
    outfile.close();
    //ofstream outfile;
    //outfile.open("sorted.fastq");
    //outfile << id[1];
    //outfile.close();
    return 0;
}