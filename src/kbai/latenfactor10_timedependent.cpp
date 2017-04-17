#include "commonheader.h"
#include <cstdlib>
#define NLAT 20

using namespace std;

float dotprod(const vector<float>& v1, const vector<float>& v2)
{
	float sum = 0.0;
	for(int i = 0 ; i < v1.size(); i++)
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
	vector<vector<vector<float>>> bst;
	vector<float> btu;
	vector<vector<float>> bt;
	float mean;
	float lr;
	float alpha;
	float lambda;

	bkmodel();
	float g(int iu, int im, int it);
	void update_param_sgd(feature &a);
};

bkmodel::bkmodel(): uuser(NUSER,0), umovie(NMOVIE,0), bm(NMOVIE,vector<float>(NLAT,0)) // 50 latent factors
 ,bu(NUSER,0), bs(NUSER,vector<float>(NLAT,0)), btu(NUSER,0), bt(NMOVIE,vector<float>(30,0)), bst(NMOVIE,vector<vector<float>>(30,vector<float>(NLAT,0)))
{
	for(auto &i : bu)  i = 1.0*rand()/RAND_MAX;  //random initialization
	for(auto &i : bs) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT; // initialization of latent factors
	for(auto &i : bm) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT;
	for(auto &i : btu) i = 1.0*rand()/RAND_MAX;
	for(auto &i : bt) for(auto &j : i)  j = 1.0*rand()/RAND_MAX;
	for(auto &i : bst) for(auto &j : i) for(auto &k : j) k = 0.0;
	mean = rand()*1.0/RAND_MAX;
	ifstream userfile("../data/userfile.dta");
	ifstream moviefile("../data/moviefile.dta");
	for(int &i : uuser) userfile >>i;
	for(int &i : umovie) moviefile >>i;
	userfile.close();
	moviefile.close();
	lr = 0.01;//set basic learning rate to be 0.01
	alpha = 4.91299;
	lambda = 0.0903385;

}

float bkmodel::g(int iu, int im, int it)
{
	float gv = 0.0;
	gv = mean + bu[iu] + dotprod(bm[im],bs[iu]) + dotprod(bst[im][it],bs[iu])  + bt[im][it]*btu[iu];
	return gv;
}

void bkmodel::update_param_sgd(feature &a)
{
//	cout << "aaa"<< endl;
	a.retrieve_feature();
//	cout << "bbb" << endl;
	int iu = a.iu -1;
	int im = a.im -1;
	int it = a.it;
	int rate = a.rate;
//	cout <<"11a" << endl;
	float error = -g(iu,im,it) + rate ;
//	cout <<"11b" <<endl;
	bu[iu] -= lr*(-2.0*error + 2.0*bu[iu]*(lambda + alpha/uuser[iu]));
//	cout << bu[17772] << endl;
//	cout <<"ccc" << endl;
	for(int i = 0; i < NLAT; i++)
	{
		bm[im][i] -= lr*(-2.0*error*bs[iu][i] + 2.0*bm[im][i]*(0.00234879 + 0.0610418/umovie[im])) ;
		bs[iu][i] -= lr*(-2.0*error*(bm[im][i]+bst[im][it][i]) + 2.0*bs[iu][i] *(0.0753211 + 3.11288/uuser[iu]));
		bst[im][it][i]  -= 0.1*lr*(-2.0*error*bs[iu][i] + 2.0*bst[im][it][i]*(0.0234879 + 0.610418/umovie[im])) ;

	}
//	cout << "ddd" << endl;
	bt[im][it] -= lr*(-2.0*error*btu[iu])+ bt[im][it] * 0.0001;
	btu[iu] -= lr*(-2.0*error*bt[im][it]) + btu[iu] * 0.001;
	mean -= lr*(-2.0*error)/100;
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
/*		if(rmse> 1e6)
		{
			cout <<iu << "\t"<< bu[iu] <<"\t" << bm[im] << endl;
			exit(12);
		}*/
	}
	rmse = sqrt(rmse);
	return rmse;
}
double compute_QUAL(bkmodel &bk, const vector<vector<int>> &testset)
{
	int im,it,iu,rate;
	int n = testset.size();
	double rmse = 0.0;
	ofstream qual("timedependent_latenfactormodel10.dta");
	for(auto & p: testset)
	{
		iu = p[0]-1;
		im = p[1]-1;
		it = p[2];
		rate = p[4];
/*		if(rmse> 1e6)
		{
			cout <<iu << "\t"<< bu[iu] <<"\t" << bm[im] << endl;
			exit(12);
		}*/
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
int main()
{

	vector<vector<int>> testset;
	vector<vector<int>> testset2;
	vector<vector<int>> testset3;
	float testscore;
	float mintestscore = 0.95;
	
	srand(0);

	ifstream infile("../data/1.dta");
	ifstream testfile("../data/2.dta");
	ifstream testfile2("../data/4.dta");
	ifstream testfile3("../data/5.dta");
	bkmodel abk;
	feature fset(infile);

	cout << "1111"<< endl;
	loadtestset(testfile,testset);
	loadtestset(testfile2,testset2);
	loadtestset(testfile3,testset3);
	cout << "2222" << endl;
	for(int j = 0; j < 300; j++)	
	{
		for(int i = 0; i < 10000000; i ++)
		{
			abk.update_param_sgd(fset);
			//cout << i << endl;
		}
		cout << compute_RMSE(abk,testset)<<endl;
		testscore = compute_RMSE(abk,testset2); 
		cout <<testscore << endl;
		if(testscore < mintestscore) 
		{
			mintestscore = testscore;
			compute_QUAL(abk,testset3);
			cout << "score recorded!"<< endl;
		}
//		cout << mean << endl;
	}

	//compute_QUAL(abk,testset3);



	
}
