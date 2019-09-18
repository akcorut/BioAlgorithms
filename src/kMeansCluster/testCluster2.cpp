#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <math.h>

using namespace std;

double sum;
double mean;
double stdv;
double var;
int i;
int j;
void norm(double **a, int row, int col)
{
    for (int j=1; j < col+1; j++) {
        for (int i=1; i < row; i++) {
            sum += a[i][j];
            mean = sum/row;
        }
        for (int i=1; i < row; i++) {
            var += ((a[i][j] - mean)*(a[i][j] - mean))/row;
            stdv = sqrt(var);
        }
        //cout << stdv << "\t" << "\n";
        for (int i=1; i < row; i++) {
            a[i][j] = a[i][j] - mean;
        }
        for (int i=1; i < row; i++) {
            a[i][j] = a[i][j]/stdv;
        } 
        sum = 0;
        mean = 0;
        var= 0;
        stdv = 0;
    }
}

double** new2DArray(int m, int n)
{
	double** resultSet;
 
	resultSet = new double*[m];
	for (int i = 0; i < m; i++)
	{
		resultSet[i] = new double[n];
	}
 
	return resultSet;
}
 
/*
	create two dimensional array, allocate the memory
*/
template<class T>
void new2DArray(T** resultSet, int m, int n)
{
	resultSet = new T*[m];
	for (int i = 0; i < m; i++)
	{
		resultSet[i] = new T[n];
	}
}


double** randCent(double **a, int row, int col, int k)
{
	srand (time(NULL));
 
	double** centroids;
	centroids = new2DArray(k, col+1);
 
	for (int j = 1; j < col+1; j++)
	{
		// get max and min value of column j
		double max = a[0][j];
		double min = a[0][j];
 
		for(int i = 1; i < row; i++)
		{
			if (a[i][j] > max)
			{
				max = a[i][j];
			}
			else if (a[i][j] < min)
			{
				min = a[i][j];
			}
		}
 
		// get the range of column j
		double range = max - min;
 
		// get k random centroids
		for(int i = 1; i < k; i++)
		{
			centroids[i][j] = min + range * rand() / ((double) RAND_MAX);
			cout << centroids[i][j] << " ";
		}
		cout << endl;
	}
 
	return centroids;
}

double disteuclid(double *vecA, double *vecB, int size)
{
	double dresult = 0.0;
 
	for (int i = 1; i < size; i++)
	{
		dresult += pow(vecA[i] - vecB[i], 2);
	}
 
	return sqrt(dresult);
}

void getNewCentroids(double **a, double** centroids, double** clustersRecord, int row, int col, int k)
{
	int* kcount = new int[k];
	memset(kcount, 0, k * sizeof(int));
	/*for (int i = 0; i < k; i++)
	{
		kcount[i] = 0;
	}*/
 
	// it seems that memset can't initialize this 2d array
	//memset(centroids, 0.0, k * n * sizeof(centroids[0][0]));
	for (int i = 1; i < k; i++)
	{
		for (int j = 1; j < col+1; j++)
		{
			centroids[i][j] = 0.0;
		}
	}
 
	for (int i = 1; i < row; i++)
	{
		int kk = (int)clustersRecord[i][1];
		kcount[kk] += 1;
 
		for (int j = 1; j < col+1; j++)
		{
			centroids[kk][j] += a[i][j];
		}
	}
 
	for (int i = 1; i < k; i++)
	{
		for (int j = 1; j < col+1; j++)
		{
			centroids[i][j] = centroids[i][j] / kcount[i];
		}
	}
}
void kMeans(double **a, int row, int col, int k)
{
	double** clustersRecord;	// clustersRecord is a two dimensional array with 2 columns, 
								// first the cluster index, second the euclid distance
	clustersRecord = new2DArray(row, 2);
 
	// get k random centroids
	double** centroids = randCent(a, row, col, k);
	bool clusterChanged = true;
 
	int count = 0;	// count is the number when to converge
 
	while (clusterChanged)
	{
		clusterChanged = false;
 
		for (int i = 1; i < row; i++)
		{
			double minDist = 100000;
			int minIndex = -1;
 
			for (int j = 1; j < k; j++)
			{
				double distJ = disteuclid(centroids[j], a[i], col);
				if (distJ < minDist)
				{
					minDist = distJ;
					minIndex = j;
				}
			}
 
			if (clustersRecord[i][1] != minIndex)
			{
				clusterChanged = true;
			}
 
			clustersRecord[i][1] = minIndex;
			clustersRecord[i][2] = minDist;
		}
 
		count++;
		cout << "The change time is : " << count << endl;
 
		// recalculate controids method 1
		getNewCentroids(a, centroids, clustersRecord, row, col, k);
 
		//// recalculate controids method 2
		//for (int cent = 0; cent < k; cent++)
		//{
		//	// get all the points in this cluster and assign centroid to mean
		//	centroids[cent] = getdataSetMean(dataSet, clustersRecord, m, n, cent);
		//}
	}
 
	// print the cluster records, the first colum is the point belong to which cluster, 
	// the second column is distance to that centroid
	for (int i = 1; i < row; i++)
	{
		cout << "clustersRecord[" << i << "][0] : " << clustersRecord[i][1] << " clustersRecord[" << i << "][1] : " << clustersRecord[i][2] << endl;
	}
}

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
        ++kk;
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

    /*for (int y = 0; y < rowCount; ++y) {
        for (int z = 0; z < colCount+1; ++z) {
            cout << Point[y][z] << "\t";
        }
    }*/
    //cout << Point[91][1];
    
    
    
    norm(Point, rowCount, colCount);
    for (int y = 0; y < rowCount; ++y) {
        for (int z = 0; z < colCount; ++z) {
            cout << Point[y][z] << "\t";
        }
    }
    //cout << std;

    //randCent(Point, rowCount, colCount, 3);

    //kMeans(Point, rowCount, colCount, 3);
    return 0;
}