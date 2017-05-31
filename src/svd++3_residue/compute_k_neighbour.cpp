#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <algorithm>
#define NMOVIE 17770


using namespace std;
vector<size_t> sort_indexes(const vector<float> &v)
{
	vector<size_t> idx(v.size());
	iota(idx.begin(),idx.end(),0);
	cout <<idx[0] << idx[1] << idx[2] << endl;
	sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2){return v[i1] > v[i2];});
	return idx;
}

int main()
{
	vector<vector<float> > corre(NMOVIE, vector<float>(NMOVIE,0));
	vector<size_t> sort_index;


	vector<float> u(4,0);

	vector<int> ratings;
	vector<int> imovies;
	
	int curr_user = 0;
	int ind1,ind2,i1,j1;
	
	
	ifstream infile("./correlationmatrixsparse.dta");
	ofstream outfile("./k_neighbour_index.dta");
	string filename;
	while(!infile.eof())
	{
//		cout << u[0] <<"\t"<< u[1]<<"\t" << u[2]<<"\t" << u[3]<<"\t"<< u[4] << endl;
		infile >> u[0] >> u[1] >> u[2] >> u[3] ;
		if(curr_user < u[0])
		{	
			cout << u[0] << endl;
			curr_user ++;
		}

		corre[u[0]-1][u[1]-1] = (u[3] * u[2]) / (100.0 + u[2]); 
		corre[u[1]-1][u[0]-1] = corre[u[0]-1][u[1]-1];

// calculate the total number of movie one rates in a specific day
	}
	for(int i = 0 ; i <NMOVIE ; i ++)
	{ 
		sort_index = sort_indexes(corre[i]);
		sort(corre[i].begin(), corre[i].end(),std::greater<float>());
		for(int j = 0; j < 100; j++)
		{
			outfile << sort_index[j] + 1 << "\t" ;
			outfile << corre[i][j] << "\t";
			cout << corre[i][j] << "\t";
		}
		outfile << endl;
		cout << endl;
	}	


	
}


