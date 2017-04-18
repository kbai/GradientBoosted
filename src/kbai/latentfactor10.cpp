#include "commonheader.h"
#include <cstdlib>
#define NLAT 10 // number of latent factors

using namespace std;

float dotprod(const vector<float>& v1, const vector<float>& v2)
{
	float sum = 0.0;
	for(int i = 0; i < v1.size(); i++)
		sum += v1[i] * v2[i];
	return sum;
}


struct feature
{
	int iu;
	int im;
	int it;
	int rate;
	ifstream & _ifile;
	feature(ifstream& a):_ifile(a){}
	void retrieve_feature();
};


void feature::retrieve_feature()
{
	int a;
	if(!(_ifile >> iu  >> im  >> it >> a>> rate))
	{
		cout << "here is the end" << endl;
		_ifile.clear();
		_ifile.seekg(0,ios::beg);// if read fails go to beginning
		_ifile >> iu  >> im  >> it >> a>> rate;
	}
}

struct bkmodel
{
	vector<int> uuser;
	vector<int> umovie;
	vector<vector<float>> bm; // 50 latent factors
	vector<float> bu;
	vector<vector<float>> bs; // 50 latent factors
	vector<float> btu;
	vector<vector<float>> bt;
	float mean;
	float _lr;
	float _alpha;
	float _lambda;

	bkmodel();
	void half_lr(); // setting learning rate
	float g(int iu, int im, int it);
	void update_param_sgd(feature &a);
};


bkmodel::bkmodel(): uuser(NUSER,0), umovie(NMOVIE,0), bm(NMOVIE,vector<float>(NLAT,0)) // 50 latent factors
 ,bu(NUSER,0), bs(NUSER,vector<float>(NLAT,0)), btu(NUSER,0), bt(NMOVIE,vector<float>(30,0))
{
	for(auto &i : bu)  i = 1.0*rand()/RAND_MAX;  //random initialization
	for(auto &i : bs) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT; // initialization of latent factors
	for(auto &i : bm) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT;
	for(auto &i : btu) i = 1.0*rand()/RAND_MAX;
	for(auto &i : bt) for(auto &j : i)  j = 1.0*rand()/RAND_MAX;
	mean = rand()*1.0/RAND_MAX;
	ifstream userfile(std::string(DATAPATH)+"/userfile.dta");
	ifstream moviefile(std::string(DATAPATH)+"/moviefile.dta");
	for(int &i : uuser) userfile >>i;
	for(int &i : umovie) moviefile >>i;
	userfile.close();
	moviefile.close();
	_lr = 0.01;//set basic learning rate to be 0.01
	_alpha = 4.91299;
	_lambda = 0.0903385;

}
void bkmodel::half_lr()
{
	_lr *= 0.5;
}

float bkmodel::g(int iu, int im, int it)
{
	float gv = 0.0;
	gv = mean + bu[iu] + dotprod(bm[im],bs[iu])  + bt[im][it]*btu[iu];
	return gv;
}

void bkmodel::update_param_sgd(feature &a)
{
	a.retrieve_feature();
	int iu = a.iu -1;
	int im = a.im -1;
	int it = a.it;
	int rate = a.rate;
	float error = -g(iu,im,it) + rate ;
	bu[iu] -= _lr*(-2.0*error + 2.0*bu[iu]*(_lambda + _alpha/uuser[iu]));
	for(int i = 0; i < NLAT; i++)
	{
		bm[im][i] -= _lr*(-2.0*error*bs[iu][i] + 2.0*bm[im][i]*(0.00234879 + 0.0610418/umovie[im])) ;
		bs[iu][i] -= _lr*(-2.0*error*bm[im][i] + 2.0*bs[iu][i] *(0.0753211 + 3.11288/uuser[iu]));
	}
	bt[im][it] -= _lr*(-2.0*error*btu[iu])+ bt[im][it] * 0.0001;
	btu[iu] -= _lr*(-2.0*error*bt[im][it]) + btu[iu] * 0.001;
	mean -= _lr*(-2.0*error)/100;
}


double compute_RMSE(bkmodel &bk, const vector<vector<int>> &testset)
{
	int im,iu,rate,it;
	int n = testset.size();
	double rmse = 0.0;
	for(auto & p: testset)
	{
		iu = p[0]-1;
		im = p[1]-1;
		it = p[2];
		rate = p[4];
		rmse +=pow( (bk.g(iu,im,it) - rate),2.0)/n;
	}
	rmse = sqrt(rmse);
	return rmse;
}
double compute_QUAL(bkmodel &bk, const vector<vector<int>> &testset, string ofilename)
{
	int im,it,iu,rate;
	int n = testset.size();
	double rmse = 0.0;
	ofstream qual(ofilename);
	for(auto & p: testset)
	{
		iu = p[0]-1;
		im = p[1]-1;
		it = p[2];
		rate = p[4];
		qual << bk.g(iu,im,it) << endl;
	}
	rmse = sqrt(rmse);
	qual.close();
	return rmse;
}
int loadtestset(ifstream &validfile, vector<vector<int>> &p)
{
	int u[5];
	vector<int> data(5,0);
	while(!validfile.eof())
	{
		if(validfile >> data[0] >> data[1] >> data[2] >> data[3] >> data[4])
		p.push_back(data);
	}
	cout << p.size() << endl;
	return 0;
}

int loadtestset2(ifstream &validfile, vector<vector<int>> &p)
{
	int u[5];
	vector<int> data(5,0);
	while(!validfile.eof())
	{
		if(validfile >> data[0] >> data[1] >> data[2])
		p.push_back(data);
	}
	cout << p.size() << endl;
	return 0;
}

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
	vector<vector<int>> testset;
	vector<vector<int>> testset2;
	vector<vector<int>> testset3;
	float testscore;
	float mintestscore = 0.95;
	float moving_ave = 3.0;
	float new_moving;
	
	srand(0);

	ifstream infile(std::string(DATAPATH)+"1.dta");
	ifstream testfile(std::string(DATAPATH)+"2.dta");
	ifstream testfile2(std::string(DATAPATH)+"4.dta");
	ifstream testfile3(std::string(DATAPATH)+"5.dta");
	bkmodel abk;
	feature fset(infile);
	loadtestset(testfile,testset);
	loadtestset(testfile2,testset2);
	loadtestset(testfile3,testset3);
	
	for(int j = 0; j < 300; j++)	
	{
		for(int i = 0; i < 10000000; i ++)
		{
			abk.update_param_sgd(fset);
		}
		cout << compute_RMSE(abk,testset)<<endl;
		testscore = compute_RMSE(abk,testset2); 
		new_moving = moving_ave* 0.9 + testscore * 0.1;
		if(new_moving >= moving_ave)
		{
			cout << "half lr" << endl;
			abk.half_lr();
		}
		moving_ave = new_moving;
		cout <<testscore << endl;
		if(testscore < mintestscore) 
		{
			mintestscore = testscore;
			compute_QUAL(abk,testset3,outfilename);
			cout << "score recorded!"<< endl;
		}
	}
}
