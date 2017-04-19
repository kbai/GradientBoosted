#include "commonheader.h"
using namespace std;
#include "feature.hh"
#include "bkmodel.hh"
#include "utils.hh"
#include <cuda.h>

bkmodel::bkmodel(): 
	uuser(NUSER,0),
   	umovie(NMOVIE,0), 
	bm(NMOVIE,0) ,
	bu(NUSER,0), 
	pu(NUSER,vector<float>(NLAT,0)), 
	pm(NMOVIE,vector<float>(NLAT,0)), // 50 latent factors
	btm(NMOVIE,vector<float>(30,0)),
	btu(NUSER,0),
	bt(30,0)
{
	for(auto &i : bu)  i = 1.0*rand()/RAND_MAX;  //random initialization
	for(auto &i : bm)  i = 1.0*rand()/RAND_MAX;  //random initialization
	for(auto &i : pu) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT; // initialization of latent factors
	for(auto &i : pm) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT;
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
	gv = mean + bu[iu] + bm[im] + dotprod(pm[im],pu[iu])  + btm[im][it]*btu[iu] + bt[it] ;
	return gv;
}

void bkmodel::update_param_sgd(feature &a)
{
	a.retrieve_feature();
	int iu = a._iu -1;
	int im = a._im -1;
	int it = a._it;
	int rate = a._rate;
	float error = -g(iu,im,it) + rate ;

	bu[iu] -= _lr*(-2.0*error + 2.0*bu[iu]*(0.0753211 + 3.11288/uuser[iu])); //no need to regulating thses terms
	bm[im] -= _lr*(-2.0*error + 2.0*bm[im]*(0.00234879 + 0.0610418/umovie[im]));

	for(int i = 0; i < NLAT; i++)
	{
//		if(alternating)
//		{
		pm[im][i] -= _lr*(-2.0*error*pu[iu][i] + 2.0*pm[im][i]*(0.00234879 + 0.0610418/umovie[im])) ;
//		}
//		else
//		{
		pu[iu][i] -= _lr*(-2.0*error*pm[im][i] + 2.0*pu[iu][i] *(0.0753211 + 3.11288/uuser[iu]));
//		}
	}
	btm[im][it] -= _lr*(-2.0*error*btu[iu]+ btm[im][it] * 0.0001);
	btu[iu] -= _lr*(-2.0*error*btm[im][it] + btu[iu] * 0.0001);
	bt[it] -= _lr*(-2.0*error)/100;
	mean -= _lr*(-2.0*error)/1000;
}




