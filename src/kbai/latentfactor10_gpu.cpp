#include "commonheader.h"
#include <cstdlib>

using namespace std;
#include "feature.hh"
#include "bkmodel.hh"
#include "bkmodel1.hh"
#include "bkmodel_gpu.hh"
#include "utils.hh"







void usage()
{
	cout << "./latentfactor10 outputfilename" << endl;
}
int main(int argc, char** argv)
{

	if(argc != 2) 
	{
		usage();
		exit(1);
	}

	string outfilename = std::string(argv[1]);
	float testscore;
	float mintestscore = 0.94;
	float moving_ave = 3.0;
	float new_moving;
	
	srand(0);

	ifstream infile(std::string(DATAPATH)+"1RAND.dta");
	ifstream testfile(std::string(DATAPATH)+"2.dta");
	ifstream testfile2(std::string(DATAPATH)+"4.dta");
	ifstream testfile3(std::string(DATAPATH)+"5.dta");
	bkmodel_gpu abk;
	float bestscore = 0.91;
	float currentscore=1.5;
	cout << "111"<< endl;
	feature fset(infile,true);  //training set , donot load data, use streaming mode
	feature tset1(testfile,true); // testsets do load data
	feature tset2(testfile2,true);
	feature tset3(testfile3,true);
	fset.load_gpu();
	tset1.load_gpu();
	tset2.load_gpu();

	cout << tset1.compute_RMSE(abk) <<endl;


	abk.loaddata(fset,tset1,tset2);
	float lr = 0.01;
	for(int i = 0 ; i < 300; i++)
	{
		abk.test(lr);
		currentscore = abk.compute_error();
		lr *=0.95;
		cout << "learning rate:" << lr << endl;
		if(currentscore <= bestscore)
		{
			abk.retrieve_gpu();
			cout <<tset1.compute_RMSE(abk) << endl;
			cout <<tset2.compute_RMSE(abk) << endl;
			tset3.compute_QUAL(abk,outfilename);
			bestscore = currentscore;
		}

	}


}
