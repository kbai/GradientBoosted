#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <algorithm>
#define NMOVIE 17770
#define NB 200

using namespace std;
void intersect(vector<int> a, vector<int> b, vector<int>& inter)
{
	sort(a.begin(), a.end());
	sort(b.begin(), b.end());
	auto i = a.begin();
	auto j = b.begin();
	while((i != a.end()) && ( j != b.end()))
	{
		if((*i) < (*j)) i++;
		if((*i) > (*j)) j++;
		if((*i)==(*j))
		{
			inter.push_back(*i);
			i++;j++;
		}
	}
}



int main()
{
	int u[6];
	vector<float> ratings;
	vector<int> imovies;
	
	int curr_user = 1;
	int ind1,ind2,i1,j1;
	vector<float> residue;
	vector<int> id;
	vector<int> cat;
	float tmpresidue;
	int currentid=5;
	vector<vector<float>> corr(NMOVIE,vector<float>(NB,0.0));
	vector<vector<int>> neig(NMOVIE,vector<int>(NB,0));
	vector<vector<int>> count(NMOVIE,vector<int>(NB,0));


	
	ifstream infile("../../data/1234map.dta");
	ifstream residf("../../data/residue.txt");
	ifstream idfile("../../data/all.idx");
	ifstream neighbour("../../data/movie_neighbour.txt");
	float rmse = 0.0;
	int ttct = 0.0;
	float rmse1 = 0.0;

	string filename;

	for(int i = 0; i < NMOVIE; i++)
	{
		int id;
		neighbour >> id;
//		cout << id << endl;
		for(int j = 0; j < NB; j++)
		{
			neighbour >> neig[i][j] >> corr[i][j] >> count[i][j];
		}	
	}
	neighbour.close();


	while(!infile.eof())
	{
		infile >>u[0]>>u[1] >> u[2] >> u[3] >> u[4]>> u[5] ;
		idfile >> currentid;
		residf >> tmpresidue;
		while(currentid ==5) idfile >> currentid;
		if(u[0] > curr_user)
		{
			if(curr_user%10000 == 0)
			{
				//	cout <<"rmse:"<< sqrt(rmse/ttct) << endl;
				//	cout <<"original rmse:" << sqrt(rmse1/ttct) << endl;
				cout << "error reduction" << sqrt(rmse/ttct) - sqrt(rmse1/ttct) << endl; 
			}

			//cout <<curr_user << endl;
			for( int i = 0 ; i < ratings.size(); i++)
			{
				if(cat[i]==4)// if sample belongs to testset
				{
					float neiresidue = 0.0;
					float weight = 0.0;
					int currentmovie = imovies[i];
					vector<int> inter;
					intersect(imovies, neig[currentmovie-1],inter);
					auto p1 = imovies.begin();
					auto p2 = inter.begin();
					auto p3 = ratings.begin();
					while((p1 != imovies.end())&&(p2 != inter.end()))
					{
						if(*p2 == *p1)
						{
							if(*p1 != currentmovie)
							{
								neiresidue += (*p3) * exp(2*corr[currentmovie-1][(*p1)-1]);
								weight += exp(2*corr[currentmovie-1][(*p1)-1]);
							}
							p1++;
							p2++;
							p3++;
						}
						else
						{
							p1++;
							p3++;
						}

					}
					if(weight <  0.5) neiresidue = 0;
					else neiresidue /= (weight + 5.0);
					//cout <<ratings[i] <<"\t" <<neiresidue <<"\t" <<count << endl; 
					cout << "rrr"<< neiresidue << endl;
					rmse += (neiresidue - ratings[i])*(neiresidue - ratings[i]);
					rmse1 += (ratings[i])*(ratings[i]);
					ttct ++;

				}
			}
			ratings.clear();
			imovies.clear();
			cat.clear();
			curr_user = u[0];
					}
		imovies.push_back(u[1]);
		ratings.push_back(tmpresidue);
		cat.push_back(currentid);

	}
	for( int i = 0 ; i < ratings.size(); i++)
	{
		if(cat[i]==4)// if sample belongs to testset
		{
			float neiresidue = 0.0;
			int count = 0;
			int currentmovie = imovies[i];
			vector<int> inter;
			intersect(imovies, neig[currentmovie-1],inter);
			auto p1 = imovies.begin();
			auto p2 = inter.begin();
			auto p3 = ratings.begin();
			while((p1 != imovies.end())&&(p2 != inter.end()))
			{
				if(*p2 == *p1)
				{
					if(*p1 != currentmovie)
					{
						neiresidue += *p3;
						count++;
					}
					p1++;
					p2++;
					p3++;
				}
				else
				{
					p1++;
					p3++;
				}

			}
			if(count == 0) neiresidue = 0;
			else neiresidue /= count;
			//cout <<ratings[i] <<"\t" <<neiresidue <<"\t" <<count << endl; 
			rmse += (neiresidue - ratings[i])*(neiresidue - ratings[i]);
			rmse1 += (ratings[i])*(ratings[i]);
			ttct ++;

		}
	}
	cout << "nratings" << ttct << endl;



}


