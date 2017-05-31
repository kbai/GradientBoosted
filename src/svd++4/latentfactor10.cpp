#include "commonheader.h"
#include <cstdlib>

using namespace std;
#include "feature.hh"
#include "bkmodel.hh"
#include "bkmodel1.hh"
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
	float mintestscore = 0.92;
	float moving_ave = 3.0;
	float new_moving;
	
	srand(0);

	ifstream infile(std::string(DATAPATH)+"1RAND.dta");
	ifstream testfile(std::string(DATAPATH)+"2.dta");
	ifstream testfile2(std::string(DATAPATH)+"4.dta");
	ifstream testfile3(std::string(DATAPATH)+"5.dta");
	bkmodel abk;
	feature fset(infile);  //training set , donot load data, use streaming mode
	feature tset1(testfile,true); // testsets do load data
	feature tset2(testfile2,true);
	feature tset3(testfile3,true);
//	cout << abk.bst[2][2][2] << endl;
	
	
	for(int j = 0; j < 300; j++)	
	{

		for(int i = 0; i < 10000000; i ++)
		{
			abk.alternating = !abk.alternating;
			abk.update_param_sgd(fset);

			cout << i << endl;
		}
		cout << tset1.compute_RMSE(abk) <<endl;

		testscore = tset2.compute_RMSE(abk); 

		new_moving = moving_ave* 0.9 + testscore * 0.1;
		if(new_moving >= moving_ave)
		{
			cout << "half lr" << endl;
			abk.half_lr();
		}
		moving_ave = new_moving;
		cout <<testscore << "minscore"<< mintestscore<<endl;
		if(testscore < mintestscore) 
		{
			mintestscore = testscore;
			tset3.compute_QUAL(abk,outfilename);
			cout << "score recorded!"<< endl;
		}
	}
}
