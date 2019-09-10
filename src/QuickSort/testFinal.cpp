#include <fstream>
#include <iostream>
#include <string>

using namespace std;

void FastqSort (string a[], string b[], string c[], string d[], int start, int end){
    int j=start; 
    int k=end;
    string tmp;
    string pivot= a[(start+end)/2];
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
    if(start<k){
        FastqSort(a, b, c, d, start, k);
    }
    if(j<end){
        FastqSort(a, b, c, d, j, end);
    }
}

int main(int argc, char* argv[])
{   
    if (argc > 1) {
        cout << "Fastq file to be sorted = " << argv[1] << endl; 
    } else {
        cout << "No fastq file detected. Please provide a proper path for the input.";
        return -1;
    }
    ifstream fastq_in(argv[1]);

    string line1, line2, line3, line4;
    string *id = new string[1000000]; 
    string *seq = new string[1000000];
    string *pls = new string[1000000];
    string *qual= new string[1000000];
    int i = 0;
    while(getline(fastq_in,line1) && getline(fastq_in,line2) && getline(fastq_in,line3) && getline(fastq_in,line4)){
        id[i] = line1;
        seq[i] = line2;
        pls[i] = line3;
        qual[i] = line4;
        i++;
    };
    //cout << id[1] << endl;
    //return 0;

    FastqSort(seq, id, pls, qual, 0, 999999);
    ofstream fastq_out;
    fastq_out.open("sorted.fastq");
    for (int y = 0; y < 1000000; y++)
    {
        fastq_out << id[y] << endl << seq[y] << endl << pls[y] << endl << qual[y] << endl;
    }
    fastq_out.close();
    return 0;
}