#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include "commonheader.h"


using namespace std;

void writedataintofile(int curr_user, vector<int>& imovie, vector<int> & itime, vector<int> & irating, vector<int>& nrating,vector<ofstream*> & ostreams)
{
	int user;
	int movie;
	int rating;
	int avetime;
	avetime = std::accumulate(itime.begin(),itime.end(),0.0)/itime.size();


	for(int i = 0 ; i < imovie.size(); i++)
	{
		(*ostreams[i]) << curr_user << "\t" << imovie[i] << "\t" << itime[i] -avetime<<"\t" << nrating[itime[i] - 1] <<"\t" << irating[i] << endl;
	}
}


int main()
{
	vector<int> imovie;
	vector<int> itime;
	vector<int> irating;
	vector<int> earlymovie(NMOVIE,3000);
	vector<int> earlyuser(NUSER, 3000);
	vector<int> moviecount(NMOVIE,0);
	vector<int> usercount(NUSER,0);
	vector<int> u(4,0);
	bool training;
	ifstream infile("../../data/all.dta");
	ifstream idfile("../../data/all.idx");
	int idx;
	int count = 0;
	while(!infile.eof())
	{
		count ++;
		infile >>u[0]>>u[1] >> u[2] >> u[3];
		if(count %10000 ==0) cout << count <<"   "<< u[0] <<  endl;

		idfile >> idx;
		earlymovie[u[1]-1] = (u[2] < earlymovie[u[1]-1]) ? u[2]:earlymovie[u[1]-1];
		earlyuser[u[0]-1] = (u[2] < earlyuser[u[0]-1]) ? u[2]:earlyuser[u[0]-1];
		usercount[u[0]-1]++;
		moviecount[u[1]-1]++;
	}

	vector<int> tmp(6,0);
	ifstream file4("../../data/4.dta");
	ofstream file4out("../../data/4addinfo.txt");
	ifstream file5("../../data/5.dta");
	ofstream file5out("../../data/5addinfo.txt");

	while(!file4.eof())
	{
		file4 >> tmp[0]>> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] ; 
		file4out << moviecount[tmp[1]-1] <<"\t" << usercount[tmp[0]-1]  <<"\t"<< tmp[2] - earlymovie[tmp[1]-1] <<"\t" << tmp[2] - earlyuser[tmp[0]-1] << endl;
	}

	while(!file5.eof())
	{
		file5 >> tmp[0]>> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] ; 
		file5out << moviecount[tmp[1]-1] <<"\t" << usercount[tmp[0]-1]  <<"\t"<< tmp[2] - earlymovie[tmp[1]-1] <<"\t" << tmp[2] - earlyuser[tmp[0]-1] << endl;
	}




}



