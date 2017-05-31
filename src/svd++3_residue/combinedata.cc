#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
int main()
{
	int id;
	float rat;
	ifstream idfile("../../data/all.idx");
	ifstream file4("./testtestsvd++3123.dta");
	ifstream file5("./testsvd++3123.dta");
	ifstream file1234("./residuetestsvd++3123.dta");
	ofstream outfile("./svd200residue_train.dta");
	vector<float> allrat;

    int count = 0;
	while(!(idfile.eof()))
	{
		idfile >> id;
		if(count %10000==0) cout << count << endl;
		count ++;
		if(id < 4) file1234 >> rat;
		if(id == 4) file4 >> rat;
		if(id == 5) file5 >> rat;
		allrat.push_back(rat);
	}
	for(float &y : allrat) outfile << y <<"\n";
	outfile.close();


}
