#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
int main()
{
	int id;
	float rat;
	ifstream idfile("../../data/all.idx");
//	ifstream file4("./train4_sum.dta");
	ifstream file5("./predict5_sum.dta");
	ifstream file1234("./trbm1234.dta");
	ofstream outfile("./tfrbmoutputpredict.dta");
	vector<float> allrat;

    int count = 0;
	while(!(idfile.eof()))
	{
		idfile >> id;
		if(count %10000==0) cout << count << endl;
		count ++;
		if(id < 5) file1234 >> rat;
//		if(id == 4) file4 >> rat;
		if(id == 5) file5 >> rat;
		allrat.push_back(rat);
	}
	for(float &y : allrat) outfile << y <<"\n";
	outfile.close();


}
