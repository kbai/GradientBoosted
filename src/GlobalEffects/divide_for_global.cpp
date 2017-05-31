#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <cmath>
#define NUSER 458293
#define NMOVIE 17770
#define DATAPATH "../../data/"

// Modified script to split the mu data

using namespace std;

int main()
{
	vector<int> u(4,0);
	int idx = -1;
//	bool training;
	ifstream infile(std::string(DATAPATH)+"all.dta");
	ifstream idfile(std::string(DATAPATH)+"all.idx");
	ofstream outfile1("../../data/1bethere.dta");
	ofstream outfile2("../../data/2bethere.dta");
	ofstream outfile3("../../data/3bethere.dta");
	ofstream outfile4("../../data/4bethere.dta");



	cout << "starting writing now" << endl;
	while(!infile.eof())
	{
//		cout << "qqq" << endl;

//		cout << u[0] << u[1] << u[2] << u[3] << endl;

		infile >>u[0]>>u[1] >> u[2] >> u[3];
//		cout << "jjj" << endl;
		idfile >> idx;
		if (idx < 4)
		{
			outfile1 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
			outfile2 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
		}
		if (idx == 4)
		{
			outfile2 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
			outfile3 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
		}
		if (idx == 5)
		{
			outfile4 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
			
		}

// calculate the total number of movie one rates in a specific day
	}

	infile.close();
	idfile.close();
	outfile2.close();
	outfile1.close();
	outfile4.close();
	outfile3.close();

	cout << "starting writing MOVIE data now" << endl;

	ifstream infile1(std::string(DATAPATH)+"mu/all.dta");
	ifstream idfile1(std::string(DATAPATH)+"mu/all.idx");
	ofstream outfile11("../../data/mu/1bethere.dta");
	ofstream outfile21("../../data/mu/2bethere.dta");
	ofstream outfile31("../../data/mu/3bethere.dta");
	ofstream outfile41("../../data/mu/4bethere.dta");

	
	while(!infile1.eof())
	{
//		cout << "qqq" << endl;

//		cout << u[0] << u[1] << u[2] << u[3] << endl;

		infile1 >>u[0]>>u[1] >> u[2] >> u[3];
//		cout << "jjj" << endl;
		idfile1 >> idx;
		if (idx < 4)
		{
			outfile11 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
			outfile21 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
		}
		if (idx == 4)
		{
			outfile21 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
			outfile31 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
		}
		if (idx == 5)
		{
			outfile41 << u[0] << "\t" << u[1] << "\t" << u[2]<<"\t" << 1 <<"\t" << u[3] << endl;
			
		}

// calculate the total number of movie one rates in a specific day
	}


	infile1.close();
	idfile1.close();
	outfile21.close();
	outfile11.close();
	outfile41.close();
	outfile31.close();










	return 0;


	
}



