#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#define NMOVIE 17770


using namespace std;


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
	vector<int> ratings;
	vector<int> imovies;
	
	int curr_user = 0;
	int ind1,ind2,i1,j1;
	
	
	ifstream infile("./1.dta");
	ofstream outfile("correlationmatrixsparse.dta");
	string filename;
	while(!infile.eof())
	{
//		cout << u[0] <<"\t"<< u[1]<<"\t" << u[2]<<"\t" << u[3]<<"\t"<< u[4] << endl;
		infile >>u[0]>>u[1] >> u[2] >> u[3] >> u[4] ;
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
			curr_user ++ ;
		}
		imovies.push_back(u[1]);
		ratings.push_back(u[4]);

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
				outfile <<i+1<<"\t"<<j+1<<"\t"<<"\t" <<correN[i][j] <<"\t"<<correR[i][j] << endl;
//				outfile <<i<<"\t"<<j<<"\t"<< correX[i][j] << "\t" << correY[i][j] << "\t" << correXY[i][j] << "\t"<<correXX[i][j] << "\t"<< correYY[i][j] << "\t"<<correN[i][j] <<"\t"<<correR[i][j] << endl;


			}

				}
	}	


	
}




