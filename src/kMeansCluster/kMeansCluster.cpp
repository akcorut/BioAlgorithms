/* Usage:
./kMeansCluster <input data>
Example:
./kMeansCluster /path/to/Bacteria.txt
*/

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <random>

using namespace std;
vector<vector<double> > v;

double sum, mean, stdv, var;
int i,j;
// Function for the normalization step
void norm(vector<vector<double> >&a, int row, int col) 
{   
    
    for (int j=0; j < col; j++) {
        // Calculate mean of each column
        for (int i=0; i < row; i++) {
            sum += a[i][j];
            mean = sum/row; 
        }
        // Calculate standard deviation of each column 
        for (int i=0; i < row; i++) {
            var += ((a[i][j] - mean)*(a[i][j] - mean))/row;
            stdv = sqrt(var);
        }
        // Normalize the column
        for (int i=0; i < row; i++) {
            a[i][j] = a[i][j] - mean;
        }
        for (int i=0; i < row; i++) {
            a[i][j] = a[i][j]/stdv;
        } 
        sum = 0;
        mean = 0;
        var= 0;
        stdv = 0;
    }
}

//Function for initializing centroids randomly
void initCentro(vector<vector<double> >&centroids, int k) 
{
    mt19937 mt_rand(time(0));
	int matrixSize = v.size();
	vector<int> cen;
	for(int i=0; i<k; i++)
		cen.push_back(mt_rand() % matrixSize);
	for(int i=0; i<k; i++)
		centroids[i] = v[cen[i]];
}

//Function for calculating euclidian distances
double eucDist(vector<double> a, vector<double> b)
{
	double temp = 0.0;
	for(int i=0; i < a.size(); i++)
	{
		temp += (a[i] - b[i])*(a[i] - b[i]);
	}
	return sqrt(temp);
}

// Function for finding new centroids
void updateCentro(vector<vector<double> >&a, vector<vector<double> >&centroids, vector<pair<double, double> > distVec, int m, int n, int k)
{
	int* kcount = new int[k];
	memset(kcount, 0, k * sizeof(int));

	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < n; j++)
		{
			centroids[i][j] = 0.0;
		}
	}

	for (int i = 0; i < m; i++)
	{
		int kk = distVec[i].second;
		kcount[kk] += 1;

		for (int j = 0; j < n; j++)
		{
			centroids[kk][j] += a[i][j];
		}
	}

	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < n; j++)
		{
			centroids[i][j] = centroids[i][j] / kcount[i];
		}
	}
}

//Function for calculating within cluster sum of squares (WCSS)
double getWCSS(vector< pair<string, pair<double, double> > > myvec)
{
    double wcss=0.0;
    for (int i=0; i<myvec.size(); i++)
    {
        wcss += (myvec[i].second.first)*(myvec[i].second.first);
    }
    return wcss;
}

double getMean(vector< pair<string, pair<double, double> > > myvec){
    double sum, mean;
    for (int i=0; i < myvec.size(); i++) 
    {
        sum += myvec[i].second.first;
        mean = sum/myvec.size();
    } 
    return mean;
}

// Function for calculating The Akaike information criterion (AIC)
double getAIC(vector<vector<double> >&a, double wcss, int k){
    double d = a[0].size();
    double aic=((2*k*d) + wcss);
    return aic;
}

// Function for calculating the Bayesian information criterion (BIC)
double getBIC(vector<vector<double> >&a, double wcss, int k){
    double d = a[0].size();
    double n= (a.size()*a[0].size());
    double bic = ((log(n)*k*d) + wcss);
    return bic;
}

/* Function for doing actual k-means clustering. This function returns a vector with 3 columns 
which first column indicates the name of the species, second colunm indicates the distance between 
that species and its cluster and the third column indicates its cluster index */
vector< pair<string, pair<double, double> > > kMeans(vector<vector<double> >&centroids, int k, vector <string> gnm)
{

	int m = v.size(); // no of rows in input file
	int n = v[0].size(); // no of columns in input file
	vector<pair<double, double> > distVec(m);
    vector<size_t> assignments(v.size());
	for(int id=0; id<=100; id++)
	{
		for(int i=0; i<m; i++)
		{
			double minDist = INT_MAX;
			int minIndex = -1;
			for(int j=0; j<k; j++){
				double disTemp = eucDist(centroids[j], v[i]);

				if(disTemp < minDist)
				{
					minDist = disTemp;
					minIndex = j;
				}
			}
			distVec[i] = (make_pair(minDist, minIndex));
            assignments[i] = minIndex;
		}
        updateCentro(v, centroids, distVec , m, n, k);
    }
    vector< pair<string, pair<double, double> > > myvec;
    for(int i = 0; i < gnm.size(); i++) 
        myvec.push_back(make_pair(gnm[i], make_pair(distVec[i].first, distVec[i].second)));
    
    return myvec;
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
    vector <string> gnm;
    string line;
    string unused;
    getline(infile, unused);
    while(getline(infile, line)) {
        stringstream linestream(line);
        string genome;
        std::vector<double> aminoData;
        double aminoValue;
        getline(linestream, genome, '\t');
        gnm.push_back(genome);
        while(linestream >> aminoValue)
        {
            aminoData.push_back(aminoValue);
        }
        v.push_back(aminoData);
    }
    cout<<"The dimension of the data set is ";
    cout<<v.size()<<" x ";
    cout<<v[0].size()<<endl;

    norm(v, v.size(), v[0].size()); // Normalize the dataset for each column
    
    ofstream fout2("outputStats.txt");
    fout2 << "k" << "\t" << "Avg. Mean" << "\t" << "WCSS" << "\t" << "AIC" << setw(10) << "BIC" << endl; 
    for(int k=1; k<=10; k++){
        vector< pair<string, pair<double, double> > > myvec;
        vector<vector<double> >centroids(k);
        initCentro(centroids, k);
        double aic, bic, wcss, means;
        for(int i=1; i<=100; i++){
            
            myvec=kMeans(centroids, k, gnm);
            double minWCSS = INT_MAX;
            double mean = getMean(myvec);
            double tempWCSS = getWCSS(myvec);
            //double aic = getAIC(v, wcss, k);
            //double bic = getBIC(v, wcss, k);
            if (tempWCSS < minWCSS){
                minWCSS = tempWCSS;
            }
        means = mean;
        wcss = minWCSS; 
        aic = getAIC(v, wcss, k);
        bic = getBIC(v, wcss, k);
        }
        fout2 << "k= " << k << "\t" << means << "\t" << wcss << "\t" << aic << "\t" << bic << endl;
    }
    fout2.close();
    
    vector<vector<double> >centroids(3);
    initCentro(centroids, 3);
    vector< pair<string, pair<double, double> > > myvec = kMeans(centroids, 3, gnm);
    ofstream fout("outputClusters.txt");
    fout << "Species" << " -------- " <<  "Minimum Distance" << " -------- " << "Cluster" << endl;  
    for (int i=0; i<myvec.size(); i++) 
    { 
        fout << myvec[i].first << " -------- " << myvec[i].second.first 
            << " -------- " << myvec[i].second.second << endl;
    }
    fout.close();
    return 0;
}