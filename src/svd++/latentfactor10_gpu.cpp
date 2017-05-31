#include "commonheader.h"
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

using namespace std;
#include "feature.hh"
#include "bkmodel.hh"
#include <thread>
#include <mutex>
#include "bkmodel1.hh"
#include "bkmodel_gpu_old.hh"
#include "utils.hh"







void usage()
{
	cout << "./latentfactor10 outputfilename" << endl;
}
void output(string& outputfilename, 
			bkmodel& abk,
			feature& f1,
			feature& f2,
			feature& f3)
{
			cout <<f1.compute_RMSE(abk) << endl;
			cout <<f2.compute_RMSE(abk) << endl;
			f3.compute_QUAL(abk,outputfilename);
			f2.compute_QUAL(abk,std::string("test")+outputfilename);
			return;

}
int main(int argc, char** argv)
{

	if(argc != 2) 
	{
		usage();
		exit(1);
	}
	char hostname[100];
	gethostname(hostname,100);
	puts(hostname);

	string outfilename = std::string(argv[1]);
	float testscore;
	float mintestscore = 0.94;
	float moving_ave = 3.0;
	float new_moving;
	thread th1;
	
	srand(0);

	ifstream infile(std::string(DATAPATH)+"1RANDmap.dta");
	ifstream testfile(std::string(DATAPATH)+"2map.dta");
	ifstream testfile2(std::string(DATAPATH)+"4map.dta");
	ifstream testfile3(std::string(DATAPATH)+"5map.dta");
	bkmodel_gpu abk;
	float bestscore = 0.500;
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
	float lr = 0.014;
	clock_t begin = clock(),end;
	float previous_score = 0.0;

	for(int i = 0 ; i < 100; i++)
	{
		begin = clock();

		abk.test(lr);
		currentscore = abk.compute_error();
//		cout <<"maximum tba:"<<*max_element(abk.bta.begin(),abk.bta.end())<<endl;
		lr*=0.91;
//
		//if(currentscore > previous_score) lr *=0.90;
		cout << "learning rate:" << lr << endl;
		if(currentscore <= bestscore)
		{
			if(th1.joinable()) th1.join();
			abk.retrieve_gpu();
			th1 = std::thread(output,std::ref(outfilename),std::ref(abk),std::ref(tset1),std::ref(tset2),std::ref(tset3));
			bestscore = currentscore;
		}
		if(i%20 == 0)
		{
			if(th1.joinable()) th1.join();
			abk.retrieve_gpu();
			th1 = std::thread(output,std::ref(outfilename),std::ref(abk),std::ref(tset1),std::ref(tset2),std::ref(tset3));
		}
		end = clock();
		cout << "Elapsed time this iteration:" << (end-begin)/CLOCKS_PER_SEC <<" seconds"<< endl; 
		previous_score =currentscore;

	}
}

