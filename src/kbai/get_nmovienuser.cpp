#include "commonheader.h"
using namespace std;
int main()
{
	vector<int> uuser(NUSER,0);
	vector<int> umovie(NMOVIE,0);
	vector<int> u(4,0);
	ifstream infile(std::string(DATAPATH)+"/all.dta");
	ofstream userfile(std::string(DATAPATH)+"/userfile.dta");
	ofstream moviefile(std::string(DATAPATH)+"/moviefile.dta");
	while(!infile.eof())
	{
		infile >> u[0] >> u[1] >> u[2] >>
		   u[3];
		uuser[u[0] -1] ++;
		umovie[u[1] -1] ++;	
	}
	for(int& i: uuser) userfile << i << endl;
	for(int& i: umovie) moviefile << i << endl;

	infile.close();
	userfile.close();
	moviefile.close();
	
	
	
}
