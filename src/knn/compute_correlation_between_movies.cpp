#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <algorithm>
#define NMOVIE 17770
#define NB 5000
#define alpha 500


using namespace std;
template<typename t1>
t1 get_matrix(vector<vector<t1>>&a, size_t x, size_t y)
{
	if(x >y) return a[x][y];
	else return a[y][x];
}

int main()
{
	vector<vector<float> > correX(NMOVIE, vector<float>(NMOVIE,0));
	vector<vector<float> > correY(NMOVIE, vector<float>(NMOVIE,0));
	vector<vector<float> > correXY(NMOVIE, vector<float>(NMOVIE,0));
	vector<vector<float> > correXX(NMOVIE, vector<float>(NMOVIE,0));
	vector<vector<float> > correYY(NMOVIE, vector<float>(NMOVIE,0));
	vector<vector<int> > correN(NMOVIE, vector<int>(NMOVIE,0));
	vector<vector<float> > correR(NMOVIE, vector<float>(NMOVIE,0));


	int u[6];
	vector<float> ratings;
	vector<int> imovies;
	
	int curr_user = 1;
	int ind1,ind2,i1,j1;
	vector<float> residue;
	vector<int> id;
	float tmpresidue;
	int currentid=5;

	
	ifstream infile("../../data/1map.dta");
	ifstream residf("../../data/residue1.txt");
	ifstream idfile("../../data/all.idx");

//	ofstream outfile("correlationmatrixsparse.dta");
	string filename;
	while(!infile.eof() && (curr_user < 10000000000) )
	{
//		cout << u[0] <<"\t"<< u[1]<<"\t" << u[2]<<"\t" << u[3]<<"\t"<< u[4] << endl;
		infile >>u[0]>>u[1] >> u[2] >> u[3] >> u[4]>> u[5] ;
		idfile >> currentid;
		residf >> tmpresidue;
		while(currentid ==5) idfile >> currentid;
		//cout << u[0] << "\t" << u[1] <<"\t" << currentid << "\t" << tmpresidue << endl;
		if(u[0] > curr_user)
		{
			cout <<curr_user << endl;
			for( int i = 0 ; i < ratings.size(); i++)
			{
				for(int j = 0; j < i ; j++)
				{
					ind1 = imovies[i] - 1;
					i1 = i;
					ind2 = imovies[j] - 1;
					j1 = j;
//					cout << ind1 <<"\t" << ind2 << endl;
					if(ind1 < ind2) 
					{
						ind1 = ind1 + ind2;
						ind2 = ind1 - ind2;
						ind1 = ind1 - ind2;
						i1 = j;
						j1 = i;
						//switch the two values;
					}
					correX[ind1][ind2] += ratings[i1];
					correY[ind1][ind2] += ratings[j1];
					correXY[ind1][ind2] += ratings[i1] * ratings[j1];
					correXX[ind1][ind2] += ratings[i1] * ratings[i1];
					correYY[ind1][ind2] += ratings[j1] * ratings[j1];
					correN[ind1][ind2]  ++;
//					cout << ind1  << "\t " << ind2 << endl; 
				}
			}
			ratings.clear();
			imovies.clear();
			curr_user = u[0];
		}
		imovies.push_back(u[1]);
		ratings.push_back(tmpresidue);
		//cout << tmpresidue<< endl;

// calculate the total number of movie one rates in a specific day
	}
/*	cout << ratings.size() <<"aaa" << imovies[0] << imovies[1] << imovies[2] << imovies[3] << endl;
	for( int i = 0 ; i < ratings.size(); i++)
	{
		for(int j = 0; j < i ; j++)
		{
			ind1 = imovies[i] - 1;
			i1 = i;
			ind2 = imovies[j] - 1;
			j1 = j;
			if(ind1 < ind2) 
			{
				ind1 = ind1 + ind2;
				ind2 = ind1 - ind2;
				ind1 = ind1 - ind2;
				i1 = j;
				j1 = i;
				//switch the two values;
			}
			
			correX[ind1][ind2] += ratings[i1];
			correY[ind1][ind2] += ratings[j1];
			correXY[ind1][ind2] += ratings[i1] * ratings[j1];
			correXX[ind1][ind2] += ratings[i1] * ratings[i1];
			correYY[ind1][ind2] += ratings[j1] * ratings[j1];
			correN[ind1][ind2]  ++;
			cout << ind1 << "\t" << ind2 << endl;
		
		}
	}
*/
	for(int i = 0 ; i <NMOVIE ; i ++)
	{
		for(int j = 0; j < NMOVIE; j++)
		{
			if(correN[i][j] == 0) 
			{
				correR[i][j] = -2;
			}
			else
			{
				correR[i][j] = (correN[i][j]*correXY[i][j]  - correX[i][j] * correY[i][j])/sqrt(correN[i][j] * correXX[i][j] - correX[i][j] * correX[i][j]+ 1e-6)/ sqrt(correN[i][j] * correYY[i][j] - correY[i][j] * correY[i][j] + 1e-6);
				correR[i][j] *= (1.0 * correN[i][j])/(correN[i][j]+alpha);
				//penalize sparse data
//				cout <<i+1<<"\t"<<j+1<<"\t"<<"\t" <<correN[i][j] <<"\t"<<correR[i][j] << endl;
//				outfile <<i<<"\t"<<j<<"\t"<< correX[i][j] << "\t" << correY[i][j] << "\t" << correXY[i][j] << "\t"<<correXX[i][j] << "\t"<< correYY[i][j] << "\t"<<correN[i][j] <<"\t"<<correR[i][j] << endl;
			}

		}
	}	
	vector<pair<int,float>> indecies(NMOVIE,std::make_pair(0,0.0));
	ofstream opt("../../data/movie_neighbour_200_only500.txt");
	for(size_t i = 0 ; i < NMOVIE; i++)
	{
		for(size_t j=0; j < NMOVIE; j++)
		{
			indecies[j].first = j+1;
			indecies[j].second = get_matrix(correR, i, j);
		}
		std::sort(indecies.begin(), indecies.end(),[](pair<int,float>&i1, pair<int,float>&i2){return i1.second > i2.second;});
		opt << i+1 << "\t" ;
		for(int k = 0 ; k < NB; k++) opt << indecies[k].first << "\t" << indecies[k].second<<"\t" << get_matrix(correN, (size_t)indecies[k].first-1, i) << "\t" ;
		opt << endl;
	}
	opt.close();
}




