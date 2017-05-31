#include "commonheader.h"
#include <stdlib.h>
#include <cstdlib>
using namespace std;
int main()
{
	ifstream ifile(std::string(DATAPATH)+"1.dta");
	ofstream ofile(std::string(DATAPATH)+"1RAND.dta");
	vector<vector<int>> dataset;
	vector<int> ent(5,0);
	int i = 0;
	cout << RAND_MAX << endl;
	
	while((ifile >> ent[0] >> ent[1] >> ent[2] >> ent[3] >> ent[4]))
	{
		i++;
		dataset.push_back(ent);
		if(i%10000 == 0) cout<< i/10000<<endl;
	}
	srand(0);
	cout <<"nelement" << dataset.size() << endl;
	int n = dataset.size();
	int a = 0;
	int start = 0;
	vector<int> tmp(5,0);
	while(n > 0)
	{
		int a = rand()%n; //a random int from 0 to n-1
		tmp = dataset[a+start];
		dataset[a+start] = dataset[start];
		dataset[start] = tmp;
		n--;
		if(n%10000 == 0) cout << n/10000 << endl;
		start ++;
	}
	for(auto &u : dataset) ofile << u[0] <<"\t" <<u[1] << "\t" 
		<<u[2] <<"\t" << u[3] << "\t" <<u[4] << endl;
	ofile.close();
	ifile.close();




	return 0;
}
