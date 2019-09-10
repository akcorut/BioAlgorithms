/* Usage:
./FastQSort <input fastq>
Example:
./FastQSort /path/to/sample1M.fastq
*/

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

using namespace std;

/* FastqSort function gets 4 string array and two integers as arguments. 
First string array must be sequence array since it is the array of comparable items. 
The other three string arrays correspond to id, plus sign and quality. 
These arrays represent 4 distinct lines (each entry in a FASTQ files consists of 4 lines) 
in input fastq file. Integer arguments represents left-most index of the subarray and 
the right-most index of the subarray.
 */

void FastqSort (string a[], string b[], string c[], string d[], int left_index, int right_index){
    int left=left_index; 
    int right=right_index;
    string tmp; // Define a temporary string objet to use it during swapping step
    string pivot= a[left_index + ((right_index - left_index) / 2)]; // Choosing the pivot as the median.
    /* Partition: Compare only sequences (string a[] should be sequence array) but swap everything */
    while(left <= right){
        while((a[left]).compare(pivot) < 0) // Compare left-most item with pivot
            left++; // when left-most element is smaller than pivot increment index
        while((a[right]).compare(pivot) > 0) // Compare right-most item with pivot
            right--; // when right-most element is bigger than pivot decrement index
        if(left <= right){ // Swapping step
            tmp=a[left]; // swap sequences
            a[left]=a[right];
            a[right]=tmp; 
            tmp=b[left]; // swap everything else (other three lines)
            b[left]=b[right];
            b[right]=tmp;
            tmp=c[left];
            c[left]=c[right];
            c[right]=tmp;
            tmp=d[left];
            d[left]=d[right];
            d[right]=tmp;
            left++;
            right--;
        }
    }
    /* Recursion: Imply quick sort recursively on both sides (left and right of the pivot) */
    if(left_index < right){
        FastqSort(a, b, c, d, left_index, right); //Sort elements on the left
    }
    if(left < right_index){
        FastqSort(a, b, c, d, left, right_index); //Sort elements on the right
    }
}

int main(int argc, char* argv[])
{   
    /* Read file from command line */
    if (argc > 1) {
        cout << "Fastq file to be sorted = " << argv[1] << endl << "Sorting... " << endl; 
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
    
    FILE * pFile;

    pFile = fopen ("out.fastq","w");
    for (int y = 0; y < 1000000; y++)
    {
        fprintf (pFile, "%s\n%s\n%s\n%s\n", id[y].c_str(), seq[y].c_str(), pls[y].c_str(), qual[y].c_str());
    }
    fclose (pFile);
    return 0;
}