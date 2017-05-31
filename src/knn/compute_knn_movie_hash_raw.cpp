#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <unordered_map>
#define NMOVIE 17770
#define Prior 3.6095162
#define NB 5000
#define UNB 5000


using namespace std;
struct dat
{
	int iuse;
	int imov;
	int icat;
	int itim;	
	float irat;

};
struct corrtabl
{
	int imov;
	float corr;
	int count;
	int itim;
	float irat;

};
typedef struct dat data;
typedef struct corrtabl corrtable;

float sigmoid(float y)
{
	return 1./(1.+exp(-y));
}


float weightfunc(float x, float tf)
{
	float p =  x*x/(1.-x*x) * tf;
//	cout << sigmoid(100.0*p-3)<< endl;
	return sigmoid(1000.0*p -10.0); 
}

int main()
{
	int u[6];
	vector<float> ratings;
	vector<int> imovies;
	int t[7];
	
	float t6;
	int curr_user = 1;
	vector<float> residue;
	vector<int> id;
	vector<int> cat;
	float tmpresidue;
	int currentid=5;
	vector<float> predictions;
	vector<vector<corrtable>> cot(NMOVIE, vector<corrtable>(UNB,corrtable()));


	
	ifstream infile("../../data/1map.dta");
	ifstream residf("../../data/1only.txt");
	ifstream idfile("../../data/all.idx");
	ifstream neighbour("../../data/1rawcorrelation.txt");
	ifstream testsetf("../../data/4map.dta");
	ifstream testr("../../data/4only.txt");
	ofstream output("knnpredict_raw_movie.txt");
	float rmse = 0.0;
	int ttct = 0;
	float rmse1 = 0.0;
	unordered_map<int,data> map; 
	vector<data> testset;

	while( !testsetf.eof()  && ! testr.eof()) 
	{
		data tmp;
		//cout << sizeof(tmp)<< endl;
		testsetf >> t[1] >> tmp.imov >> tmp.itim >> t[3] >> t[4] >> t[5] ;
		tmp.iuse = t[1];
		testr >>tmp.irat; 
		testset.push_back(tmp);
	}
	testsetf.close();
	int iii = 0;
	cout << testset[iii].iuse << "user" << endl;
	cout << testset.size() << " test data loaded!" << endl;

	string filename;

	for(int i = 0; i < NMOVIE; i++)
	{
		int id;
		neighbour >> id;
		assert(id == i+1);
		for(int j = 0; j < UNB; j++)
		{
			if(j < UNB) neighbour >> cot[i][j].imov >> cot[i][j].corr >> cot[i][j].count;

		}	
		neighbour.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//		sort(cot[i].begin(), cot[i].end(), [](corrtable &x, corrtable &y){return x.corr > y.corr;});
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
			while(testset[iii].iuse == curr_user)
			{
					float neiresidue = 0.0;
					float weight = 0.0;
					int currentmovie = testset[iii].imov;
					vector<corrtable> inter;
					int ic  = 0;
					auto p = cot[currentmovie-1].begin();
					while(ic < 30 && p != cot[currentmovie-1].end())
					{
						auto f = map.find((*p).imov);
						if(f!= map.end())	
						{
							(*p).itim = (*f).second.itim;
							(*p).irat = (*f).second.irat;
							inter.push_back(*p); // get 30 neighbour from the hash table
							ic ++;
						}	
						p++;
					}
					for(auto & a: inter)
					{
						float tf = (1./(1.0+ 0.002*(abs(a.itim - testset[iii].itim))));
						//cout << a.itim - testset[iii].itim << endl;
						neiresidue += (a.irat) * weightfunc((a.corr),tf);
						//cout <<a.corr<<"\t"<< weightfunc(a.corr,tf)<<endl;
						weight += weightfunc(a.corr,tf);

					}
					if(weight < -1.0) neiresidue = Prior;
					else neiresidue = (neiresidue + 1.0*Prior)/(weight+1.0); 
					rmse += (neiresidue - testset[iii].irat)*(neiresidue - testset[iii].irat);
					rmse1 += testset[iii].irat*testset[iii].irat;
					ttct ++;
					iii++;
					predictions.push_back(neiresidue);
			}

			
			map.clear();
			curr_user = u[0];
		}
		data tmp{0};
		tmp.irat = tmpresidue;
		tmp.icat = currentid;
		tmp.itim = u[2];
		map.insert(std::make_pair(u[1],tmp));

	}
	while(testset[iii].iuse == curr_user && iii < testset.size())
	{
			float neiresidue = 0.0;
			float weight = 0.0;
			int currentmovie = testset[iii].imov;
			vector<corrtable> inter;
			int ic  = 0;
			auto p = cot[currentmovie-1].begin();
			while(ic < 30 && p != cot[currentmovie-1].end())
			{
				auto f = map.find((*p).imov);
				if(f!= map.end())	
				{
					(*p).itim = (*f).second.itim;
					(*p).irat = (*f).second.irat;
					inter.push_back(*p); // get 30 neighbour from the hash table
					ic ++;
				}	
				p++;
			}
			for(auto & a: inter)
			{
				float tf = (1./(1.0+ 0.002*(abs(a.itim - testset[iii].itim))));
				neiresidue += (a.irat) * weightfunc((a.corr),tf);
				weight += weightfunc(a.corr,tf);

			}
			if(weight < -1.0) neiresidue = Prior;
			else neiresidue /= (weight+1.0); 
			rmse += (neiresidue - testset[iii].irat)*(neiresidue - testset[iii].irat);
			rmse1 += testset[iii].irat*testset[iii].irat;
			ttct ++;
			iii++;
			predictions.push_back(neiresidue);
	}

	for(auto &a: predictions)  output << a << endl;

}


