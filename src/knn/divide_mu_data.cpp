#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <cmath>
#define NUSER 458293
#define NMOVIE 17770
#define DATAPATH "../../data/mu/"

// Modified script to split the mu data

using namespace std;

void writedataintofile(int curr_movie, vector<int>& iuser, vector<int> & itime, vector<int> & irating, vector<int>& nrating,vector<ofstream*> & ostreams)
{
	int user;
	int movie;
	int rating;
	int avetime;
	avetime = std::accumulate(itime.begin(),itime.end(),0.0)/itime.size();


	for(int i = 0 ; i < iuser.size(); i++)
	{
		(*ostreams[i]) << iuser[i] << "\t" << curr_movie << "\t" << itime[i]/75<<"\t" << nrating[itime[i] - 1] <<"\t" << irating[i] << endl;
	}
}


int main()
{
	vector<int> iuser;
	vector<int> itime;
	vector<int> irating;
	vector<int> nrating(2243,0);
	vector<int> u(4,0);
	vector<ofstream*> ostreams;
	int curr_movie = 1;
	int counter = 0;
	int averagetime = 0;
	int idx = -1;
	bool training;
	ifstream infile(std::string(DATAPATH)+"all.dta");
	ifstream idfile(std::string(DATAPATH)+"all.idx");
	ofstream outfile[5];
	string filename;
	for(int i = 1; i < 6; i++)
	{
		filename = std::string(DATAPATH) + std::to_string(i) + ".dta";
		
		outfile[i-1].open(filename);
	}
	while(!infile.eof())
	{
//		cout << "qqq" << endl;

//		cout << u[0] << u[1] << u[2] << u[3] << endl;

		infile >>u[0]>>u[1] >> u[2] >> u[3];
//		cout << "jjj" << endl;
		idfile >> idx;
//		cout << "iii" << endl;
		if(u[1] > curr_movie)
		{
			cout << curr_movie << endl;


			writedataintofile(curr_movie,iuser,itime,irating,nrating,ostreams);
			std::fill(nrating.begin(),nrating.end(),0);
			counter = 0;
			iuser.clear();
			itime.clear();
			irating.clear();
			ostreams.clear();
			curr_movie++;
		//	if(curr_movie == 3) exit(0);
		}
//		cout << "ppp" << endl;
		iuser.push_back(u[1]);
		itime.push_back(u[2]);
		irating.push_back(u[3]);
		ostreams.push_back(&(outfile[idx-1]));
//		cout << u[2] << endl;
		nrating[u[2] - 1] ++; // accumulating rating of a specific user at day u[2]
		
		counter ++;
// calculate the total number of movie one rates in a specific day
	}
	writedataintofile(curr_movie,iuser,itime,irating,nrating,ostreams);



	
}



