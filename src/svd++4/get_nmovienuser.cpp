#include "commonheader.h"
using namespace std;
int main()
{
	vector<int> uuser(NUSER,0);
	vector<int> umovie(NMOVIE,0);
	vector<float> avetime(NUSER,0);
	vector<int> u(4,0);
	ifstream infile(std::string(DATAPATH)+"/all.dta");
	ofstream userfile(std::string(DATAPATH)+"/userfile.dta");
	ofstream moviefile(std::string(DATAPATH)+"/moviefile.dta");
	ofstream tufile(std::string(DATAPATH)+"/averagetime.dta");
	while(!infile.eof())
	{
		infile >> u[0] >> u[1] >> u[2] >>
		   u[3];
		uuser[u[0] -1] ++;
		umovie[u[1] -1] ++;	
		avetime[u[0] -1] += u[2];
			
	}
	for(int i = 0 ; i < avetime.size(); i++) avetime[i] /= uuser[i];
	for(int& i: uuser) userfile << i << endl;
	for(int& i: umovie) moviefile << i << endl;
	for(auto &i : avetime) tufile << i << endl;
	infile.close();
	userfile.close();
	moviefile.close();
	tufile.close();
	
	
}
