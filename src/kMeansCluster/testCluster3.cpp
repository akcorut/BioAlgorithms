#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <vector>
#include <cstdlib>

using namespace std;
vector<vector<double> > v;

double sum, mean, stdv, var;
int i,j;
void norm(vector<vector<double> >&a, int row, int col)
{
    for (int j=0; j < col; j++) {
        for (int i=0; i < row; i++) {
            sum += a[i][j];
            mean = sum/row;
        }
        for (int i=0; i < row; i++) {
            var += ((a[i][j] - mean)*(a[i][j] - mean))/row;
            stdv = sqrt(var);
        }
        //cout << stdv << "\t" << "\n";
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

void initCentro(vector<vector<double> >&centroids, int k)
{
    srand(time(NULL));
	int matrixSize = v.size();
	vector<int> cen;
	for(int i=0; i<k; i++)
		cen.push_back(rand() % matrixSize);
	for(int i=0; i<k; i++)
		centroids[i] = v[cen[i]];
}

double eucDist(vector<double> a, vector<double> b)
{
	double temp = 0.0;
	for(int i=0; i < a.size(); i++)
	{
		temp += (a[i] - b[i])*(a[i] - b[i]);
	}
	return sqrt(temp);
}

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

double getWCSS(vector< pair<string, pair<double, double> > > myvec)
{
    double wcss=0.0;
    for (int i=0; i<myvec.size(); i++)
    {
        wcss += myvec[i].second.first;
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

double getAIC(vector<vector<double> >&a, double wcss, int k){
    double d = a[0].size();
    double aic=((2*k*d) + wcss);
    return aic;
}

double getBIC(vector<vector<double> >&a, double wcss, int k){
    double d = a[0].size();
    double n= (a.size()*a[0].size());
    double bic = ((log(n)*k*d) + wcss);
    return bic;
}

void kMeans(vector<vector<double> >&centroids, int k, vector <string> gnm)
{

	int m = v.size(); // no of rows in input file
	int n = v[0].size(); // no of columns in input file
    //vector<double> tempV;
	vector<pair<double, double> > distVec(m);
    vector<size_t> assignments(v.size());
    //vector<pair<string, double> >clusters(m);
	for(int id=0; id<=300; id++)
	{
		for(int i=0; i<m; i++)
		{
			/*for(int j=0; j<n; j++)
			{
				tempV.push_back(v[i][j]);
			}*/
			double minDist = INT_MAX;
			int minIndex = 0;
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
			//tempV.clear();
		}
        updateCentro(v, centroids, distVec , m, n, k);
    }
    /*
    for(int i = 0; i < clusters.size(); i++)
        cout << clusters[i].first << ", " << clusters[i].second << endl;
    */
    /*for(int i=0; i<m; i++)
        cout << assignments[i] << endl;*/


    /*cout << "Species" << ", " << "Distance" << ", " << "Cluster" << endl; 
    for(int i = 0; i < distVec.size(); i++)
    {
        cout << gnm[i] << ", " <<  distVec[i].first << ", " << distVec[i].second << endl;
    }*/
    vector< pair<string, pair<double, double> > > myvec;
    for(int i = 0; i < gnm.size(); i++) 
        myvec.push_back(make_pair(gnm[i], make_pair(distVec[i].first, distVec[i].second)));
    
    /*for (int i=0; i<myvec.size(); i++) 
    { 
        cout << myvec[i].first << ", " << myvec[i].second.first 
            << ", " << myvec[i].second.second << endl; 
    }*/
    /*double wcss=0.0;
    for (int i=0; i<myvec.size(); i++) 
        wcss += myvec[i].second.first;
    cout << wcss;*/
    cout << "k= " << k << "\t" << getMean(myvec) << "\t" << getWCSS(myvec) << "\t" << getAIC(v, getWCSS(myvec), k) << "\t" << getBIC(v, getWCSS(myvec), k) << endl;
    
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

    norm(v, v.size(), v[0].size());
    /*for(int i=0; i<v.size(); i++)
    {
		for(int j=0; j<v[0].size(); j++)
			cout<<v[i][j]<<" ";
		cout<<endl;
    }cout<<endl;*/
    cout << "k" << "\t" << "Avg. Mean" << "\t" << "WCSS" << "\t" << "AIC" << "\t" << "BIC" << endl; 
    for(int j=1; j<11; j++){
        vector<vector<double> >centroids(j);
        initCentro(centroids, j);
    /*for(int i=0; i<centroids.size(); i++)
    {
		for(int j=0; j<centroids[0].size(); j++)
			cout<<centroids[i][j]<<" ";
        cout<<endl;
    }cout<<endl;*/
    //for (int i=0; i<gnm.size(); i++)     
    //    cout << gnm[i] << "\n";
        kMeans(centroids, j, gnm);
    }
    return 0;
}