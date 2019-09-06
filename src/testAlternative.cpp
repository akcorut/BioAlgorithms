#include <fstream>
#include <iostream>
#include <string>

using namespace std;

/* FastqSort function gets 4 string array and two integers as arguments. 
First string array must be sequence array since it is the array of comparable items. 
The other three string arrays correspond to id, plus sign and quality. 
These arrays represent 4 distinct lines (each entry in a FASTQ files consists of 4 lines) 
in input fastq file. Integer arguments represents left-most index of the subarray and 
the right-most index of the subarray.
 */

int median3(string arr[], int left, int right) {
    
    int center = (left + right)/2;
    if(arr[center].compare(arr[ left ]) < 0 )
        swap( arr[ left ], arr[ center ] );
    if( arr[ right ].compare(arr[ left ]) < 0 )
        swap( arr[ left ], arr[ right ] );
    if( arr[ right ].compare(arr[ center ]) < 0 )
        swap( arr[ center ], arr[ right ] );
    
    swap( arr[ center ], arr[ right-1]);
    int pivot=right-1;
    return pivot;
}

void FastqSort (string a[], string b[], string c[], string d[], int start, int end){
    int j=start; 
    int k=end;
    string tmp; // Define a temporary string objet to use it during swapping step
    string pivot= a[median3(a, start, end)]; // Median-of-Three Partitioning: Choosing the pivot as the median of the left, right, and center elements. 
    /* Partition: Compare only sequences (string a[] should be sequence array) but swap everything */
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
    /* Recursion: Imply quick sort recursively on both sides (left and right of the pivot) */
    if(start<k){
        FastqSort(a, b, c, d, start, k); //Sort elements on the left
    }
    if(j<end){
        FastqSort(a, b, c, d, j, end); //Sort elements on the right
    }
}

int main(int argc, char* argv[])
{   
    /* Read file from command line */
    if (argc > 1) {
        cout << "Fastq file to be sorted = " << argv[1] << endl; 
    } else {
        cout << "No fastq file detected. Please provide a proper path for the input.";
        return -1;
    }
    ifstream fastq_in(argv[1]);

    /* Get each line and store them into string arrays 
    by their type (id, sequence, plus sign and quality) */
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

    /* Run quick sort function on input file and write it to an output file */
    FastqSort(seq, id, pls, qual, 0, 999999);
    ofstream fastq_out;
    fastq_out.open("sorted_Alt.fastq");
    for (int y = 0; y < 1000000; y++)
    {
        fastq_out << id[y] << endl << seq[y] << endl << pls[y] << endl << qual[y] << endl;
    }
    fastq_out.close();
    return 0;
}