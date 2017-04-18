#include "commonheader.h"
using namespace std;
#include "feature.hh"
#include "bkmodel.hh"
#include "utils.hh"
#define NLAT 10

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





