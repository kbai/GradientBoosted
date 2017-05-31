#include "commonheader.h"
#include <cstdlib>
using namespace std;
int sgd(vector<int> & uuser, vector<int> & umovie,float &mean, vector<float> &bu, vector<float> &bm,vector<float>& bs,vector<vector<float>>& bt,vector<float>& btu, ifstream &ofile,float lambda, float alpha ,float lr)
{
	int iu;
	int im;
	int it;
	int a;
	int rate;
	float error;
	if(!(ofile >> iu  >> im  >> it >> a>> rate))
	{
		cout << "here is the end" << endl;
		ofile.clear();
		ofile.seekg(0,ios::beg);// if read fails go to beginning
		ofile >> iu  >> im  >> it >> a>> rate;

	}
	//cout << iu <<"\t"<< im <<"\t" << rate << endl;
	iu--;
	im--;
	error = -(mean + bu[iu] + bm[im]*(1+bs[iu]) + bt[im][it]*btu[iu]) + rate ;
	bu[iu] -= lr*(-2.0*error + 2.0*bu[iu]*(lambda + alpha/uuser[iu]));
//	cout << bu[17772] << endl;
	bm[im] -= lr*(-2.0*error*(1+bs[iu]) + 2.0*bm[im]*(0.00234879 + 0.0610418/umovie[im])) ;
	bs[iu] -= lr*(-2.0*error*bm[im] + 2.0*bs[iu] *(0.0753211 + 3.11288/uuser[iu]));
	bt[im][it] -= lr*(-2.0*error*btu[iu])+ bt[im][it] * 0.0001;
	btu[iu] -= lr*(-2.0*error*bt[im][it]) + btu[iu] * 0.001;
	mean -= lr*(-2.0*error)/100;
	return 0;
}
double compute_RMSE(vector<vector<int>> &testset,vector<float> &bu,vector<float> &bm,vector<float> &bs,vector<vector<float>> & bt,vector<float> & btu, float mean)
{
	int im,iu,rate,it;
	int n = testset.size();
	double rmse = 0.0;
	for(vector<int> & p: testset)
	{
		iu = p[0]-1;
		im = p[1]-1;
		it = p[2];
		rate = p[4];
//		cout << iu << "\t" << im  << "\t" << rate << endl;
//
//		cout << rmse << endl;
		rmse +=pow( (mean + bu[iu] + bm[im]*(1+bs[iu])+bt[im][it]*btu[iu] -rate),2.0)/n;
		if(rmse> 1e6)
		{
			cout <<iu << "\t"<< bu[iu] <<"\t" << bm[im] << endl;
			exit(12);
		}
//		rmse = sqrt(rmse);
//		cout << testset.size() << endl;
//		exit(1);
	}
	rmse = sqrt(rmse);
	return rmse;
}
double compute_QUAL(vector<vector<int>> &testset,vector<float> &bu,vector<float> &bm,vector<float> &bs,vector<vector<float>> &bt,vector<float> &btu, float mean)
{
	int im,it,iu,rate;
	int n = testset.size();
	double rmse = 0.0;
	ofstream qual("qualresults_baseline2.dta");
	for(vector<int> & p: testset)
	{
		iu = p[0]-1;
		im = p[1]-1;
		it = p[2];
		rate = p[4];
//		cout << iu << "\t" << im  << "\t" << rate << endl;
//
//		cout << rmse << endl;
		rmse +=pow( (mean + bu[iu] + bm[im]*(1+bs[iu])+bt[im][it]*btu[iu] -rate),2.0)/n;
		if(rmse> 1e6)
		{
			cout <<iu << "\t"<< bu[iu] <<"\t" << bm[im] << endl;
			exit(12);
		}
		qual << mean + bu[iu] + bm[im]*(1+bs[iu])+bt[im][it]*btu[iu] << endl;
//		rmse = sqrt(rmse);
//		cout << testset.size() << endl;
//		exit(1);
	}
	rmse = sqrt(rmse);
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
}
int main()
{

	float mean = 0.0;
	vector<int> uuser(NUSER,0);
	vector<int> umovie(NMOVIE,0);
	vector<vector<int>> testset;
	vector<vector<int>> testset2;
	vector<vector<int>> testset3;
	
	vector<float> bm(NMOVIE,0);
	vector<float> bu(NUSER,0);
	vector<float> bs(NUSER,0);
	vector<float> btu(NUSER,0);
	vector<vector<float>> bt(NMOVIE,vector<float>(30,0));
	srand(0);

	for(float &i : bu) i = 1.0*rand()/RAND_MAX;  //random initialization
	for(float &i : bs) i = 1.0*rand()/RAND_MAX;
	for(float &i : bm) i = 1.0*rand()/RAND_MAX;
	for(float &i : btu) i = 1.0*rand()/RAND_MAX;
	mean = rand()*1.0/RAND_MAX;
	float alpha = 4.91299;
	float lbd = 0.0903385;
	ifstream infile("../data/1.dta");
	ifstream userfile("../data/userfile.dta");
	ifstream moviefile("../data/moviefile.dta");
	ifstream testfile("../data/2.dta");
	ifstream testfile2("../data/4.dta");
	ifstream testfile3("../data/5.dta");
	cout << "1111"<< endl;
	loadtestset(testfile,testset);
	loadtestset(testfile2,testset2);
	loadtestset(testfile3,testset3);
	cout << "2222" << endl;
	for(int &i : uuser) userfile >>i;
	for(int &i : umovie) moviefile >>i;
	for(int j = 0; j < 300; j++)	
	{
		for(int i = 0; i < 10000000; i ++)
		sgd(uuser,umovie,mean,bu,bm,bs,bt,btu,infile,lbd,alpha,0.01);
		cout << compute_RMSE(testset,bu,bm,bs,bt,btu, mean)<<endl;
		cout << compute_RMSE(testset2,bu,bm,bs,bt,btu,  mean) << endl;
//		cout << mean << endl;
	}

	compute_QUAL(testset3, bu, bm,bs,bt,btu, mean);



	
}
