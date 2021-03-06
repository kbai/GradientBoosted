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
		(*ostreams[i]) << curr_user << "\t" << imovie[i] << "\t" << itime[i]<<"\t" << nrating[itime[i] - 1] <<"\t" << irating[i] << endl;
	}
}


int main()
{
	vector<int> imovie;
	vector<int> itime;
	vector<int> irating;
	vector<int> nrating(2243,0);
	vector<int> u(4,0);
	vector<ofstream*> ostreams;
	srand(0);
	int curr_user = 1;
	int counter = 0;
	int averagetime = 0;
	int idx = -1;
	bool training;
	ifstream infile(std::string(DATAPATH)+"all.dta");
	ifstream idfile(std::string(DATAPATH)+"all.idx");
	ofstream outfile[2];
	string filename;
	outfile[0].open( std::string(DATAPATH)  + "1234.dta");
	outfile[1].open( std::string(DATAPATH)  + "55.dta");



	while(!infile.eof())
	{
//		cout << "qqq" << endl;

//		cout << u[0] << u[1] << u[2] << u[3] << endl;

		infile >>u[0]>>u[1] >> u[2] >> u[3];
//		cout << "jjj" << endl;
		idfile >> idx;
//		cout << "iii" << endl;
		if(u[0] > curr_user)
		{
			cout << curr_user << endl;


			writedataintofile(curr_user,imovie,itime,irating,nrating,ostreams);
			std::fill(nrating.begin(),nrating.end(),0);
			counter = 0;
			imovie.clear();
			itime.clear();
			irating.clear();
			ostreams.clear();
			curr_user++;
		//	if(curr_user == 3) exit(0);
		}
//		cout << "ppp" << endl;
		imovie.push_back(u[1]);
		itime.push_back(u[2]);
		irating.push_back(u[3]);
		if(idx < 5) ostreams.push_back(&(outfile[0]));
		else ostreams.push_back(&(outfile[1]));
//		cout << u[2] << endl;
		nrating[u[2] - 1] ++; // accumulating rating of a specific user at day u[2]
		
		counter ++;
// calculate the total number of movie one rates in a specific day
	}
	writedataintofile(curr_user,imovie,itime,irating,nrating,ostreams);



	
}



