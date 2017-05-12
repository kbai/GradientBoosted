#include <iostream>

#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <assert.h>
#include <algorithm>
#define NMOVIE 17770
#define NB 5000
#define UNB 100

using namespace std;
void intersect(vector<int> b, const vector<int> &a,const vector<float> &c ,vector<pair<int,float>>& inter)
{
	vector<pair<int,float>> d;
	for(int i = 0; i < a.size(); i++)
	{
		d.push_back(std::make_pair(a[i],c[i]));
	}
	sort(d.begin(), d.end(), [](pair<int,float>& m, pair<int,float>&n){return m.first < n.first;});
//	sort(b.begin(), b.end());
	auto i = d.begin();
	auto j = b.begin();
//	cout << "check if order is correct!" << endl;
//	for(auto o: d) cout << o.first <<"\t" << o.second <<  endl;
	while((i != d.end()) && ( j != b.end()))
	{
//		cout << (*i).first <<"\t" << *j << endl;
		if(((*i).first) < (*j)) i++;
		else if(((*i).first) > (*j)) j++;
		else if(((*i).first)==(*j))
		{
			inter.push_back(*i);

			i++;j++;
		}
	}
//	cout<< inter.size() << endl;
//	cout << "done!" << endl;
//
}
float sigmoid(float y)
{
	return 1./(1.+exp(-y));
}


float weightfunc(float x)
{
	float p =  x*x/(1.-x*x);
	return sigmoid(100.0*p -3); 
}

int main()
{
	int u[6];
	vector<float> ratings;
	vector<int> imovies;
	int t[7];
	
	float t6;
	int curr_user = 1;
	int ind1,ind2,i1,j1;
	vector<float> residue;
	vector<int> id;
	vector<int> cat;
	float tmpresidue;
	int currentid=5;
	vector<vector<float>> corr(NMOVIE,vector<float>(UNB,0.0));
	vector<vector<int>> neig(NMOVIE,vector<int>(UNB,0));
	vector<vector<int>> count(NMOVIE,vector<int>(UNB,0));
	vector<float> predictions;


	
	ifstream infile("../../data/1map.dta");
	ifstream residf("../../data/residue1.txt");
	ifstream idfile("../../data/all.idx");
	ifstream neighbour("../../data/movie_neighbour_200_only500.txt");
	ifstream testset("../../data/5map.dta");
	ifstream testr("../../data/residue1.txt");
	ofstream output("knnpredict_withprobe.txt");
	float rmse = 0.0;
	int ttct = 0;
	float rmse1 = 0.0;
	vector<int> tuser;
	vector<int> tmovie;
	vector<float> tresidue;
	vector<int> times;
	vector<int> ttime;
	while( !testset.eof() ) 
	{
		testset >> t[0] >> t[1] >> t[2] >> t[3] >> t[4] >> t[5] ;
		testr >>t6; 
		tuser.push_back(t[0]);
		tmovie.push_back(t[1]);
		ttime.push_back(t[2]);
		tresidue.push_back(t6);

	}
	int iii = 0;
	cout << tuser.size() << " test data loaded!" << endl;

	string filename;
	float trash;

	for(int i = 0; i < NMOVIE; i++)
	{
		int id;
		neighbour >> id;
		assert(id == i+1);
//		cout << id << endl;
		for(int j = 0; j < UNB; j++)
		{
			if(j < UNB) neighbour >> neig[i][j] >> corr[i][j] >> count[i][j];
//			else neighbour >> trash >> trash >> trash;// discard the rest
		}	
		neighbour.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	neighbour.close();


	while(!infile.eof())
	{
		infile >>u[0]>>u[1] >> u[2] >> u[3] >> u[4]>> u[5] ;
		idfile >> currentid;
		residf >> tmpresidue;
		while(currentid > 3) idfile >> currentid;
		if(u[0] > curr_user)
		{
			if(curr_user%10000 == 0)
			{
				cout << curr_user << endl;
					cout <<"rmse:"<< sqrt(rmse/ttct) << endl;
					cout <<"original rmse:" << sqrt(rmse1/ttct) << endl;
				cout << "error reduction" << sqrt(rmse/ttct) - sqrt(rmse1/ttct) << endl; 
			}
			while(tuser[iii] == curr_user)
			{
					float neiresidue = 0.0;
					float weight = 0.0;
					int currentmovie = tmovie[iii];
					vector<pair<int,float>> inter;
					intersect(imovies, neig[currentmovie-1],corr[currentmovie-1],inter);
					//cout << inter.size() << endl;
					auto p1 = imovies.begin();
					auto p2 = inter.begin();
					auto p3 = ratings.begin();
					auto p4 = times.begin();
					while((p1 != imovies.end())&&(p2 != inter.end()))
					{
						if((*p2).first == *p1)
						{
							if(*p1 != currentmovie)
							{
		//						if((*p2).second > 0.4)
		//						{
									float tf = (1./(1.0+ 0.00*(abs((*p4) - ttime[iii]))));
//									cout << tf << endl;
									neiresidue += (*p3) * weightfunc(((*p2).second*tf));//(*p3) * exp(2*corr[currentmovie-1][(*p1)-1]);
								//cout << *p3 << endl;
									weight +=  weightfunc((*p2).second*tf);//exp(2*corr[currentmovie-1][(*p1)-1]);
		//						}
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
//					cout << "1111" << endl;
					assert(p2 == inter.end());// all movies in the intersection should be found
					//cout<< "size of neighbours:" << inter.size() << endl;

					if(weight < 0.1) neiresidue = 0;
					else neiresidue /= (weight + 2.0);
					//cout <<ratings[i] <<"\t" <<neiresidue <<"\t" <<count << endl; 
					rmse += (neiresidue - tresidue[iii])*(neiresidue - tresidue[iii]);
					rmse1 += (tresidue[iii])*(tresidue[iii]);
					ttct ++;
					iii++;
					predictions.push_back(neiresidue);
		//		output<< neiresidue << endl;
			}

			
			ratings.clear();
			imovies.clear();
			cat.clear();
			times.clear();
			curr_user = u[0];
		}
		imovies.push_back(u[1]);
		ratings.push_back(tmpresidue);
		cat.push_back(currentid);
		times.push_back(u[2]);

	}
			while(tuser[iii] == curr_user)
			{
					float neiresidue = 0.0;
					float weight = 0.0;
					int currentmovie = tmovie[iii];
					vector<pair<int,float>> inter;
					intersect(imovies, neig[currentmovie-1],corr[currentmovie-1],inter);
					auto p1 = imovies.begin();
					auto p2 = inter.begin();
					auto p3 = ratings.begin();
					while((p1 != imovies.end())&&(p2 != inter.end()))
					{
						if((*p2).first == *p1)
						{
							if(*p1 != currentmovie)
							{
								neiresidue += *p3;//(*p3) * exp(2*corr[currentmovie-1][(*p1)-1]);
								//cout << *p3 << endl;
								weight += 1.0;//exp(2*corr[currentmovie-1][(*p1)-1]);
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
					assert(p2 == inter.end());// all movies in the intersection should be found
					//cout<< "size of neighbours:" << inter.size() << endl;

					if(weight < 0.5) neiresidue = 0;
					else neiresidue /= (weight);
					//cout <<ratings[i] <<"\t" <<neiresidue <<"\t" <<count << endl; 
					rmse += (neiresidue - tresidue[iii])*(neiresidue - tresidue[iii]);
					rmse1 += (tresidue[iii])*(tresidue[iii]);
					ttct ++;
					iii++;
					predictions.push_back(neiresidue);
	
								}
			cout << predictions.size() <<"done!" << endl;
			for(auto & a: predictions) output << a << endl;


}


