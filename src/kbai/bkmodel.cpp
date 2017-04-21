#include "commonheader.h"
using namespace std;
#include <assert.h>
#include "feature.hh"
#include "bkmodel.hh"
#include "utils.hh"
#include <cuda.h>

float u(float a)
{
	if(a>0) return pow(a,0.4);
	else return -pow(abs(a),0.4);
}

bkmodel::bkmodel(): 
	uuser(NUSER,0),
   	umovie(NMOVIE,0), 
	bm(NMOVIE,0) ,
	bu(NUSER,0), 
	pu(NUSER,vector<float>(NLAT,0)), 
	pm(NMOVIE,vector<float>(NLAT,0)), // 50 latent factors
	pu1(NUSER,vector<float>(NLAT,0.0)),

	btm(NMOVIE,vector<float>(30,0)),
	btu(NUSER,0),
	bt(30,0),
	bf(8,0),
//	pfm(NMOVIE,vector<vector<float>>(8,vector<float>(NLAT,0)))
	bfm(NMOVIE,vector<float>(8,0))
{
	cout << "Number of latent factors" << NLAT << endl;
	for(auto &i : bu)  i = 1.0*rand()/RAND_MAX;  //random initialization
	for(auto &i : bm)  i = 1.0*rand()/RAND_MAX;  //random initialization
	for(auto &i : pu) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT; // initialization of latent factors
	for(auto &i : pm) for(auto &j : i)  j = 1.0*rand()/RAND_MAX/NLAT;
	//mean = rand()*1.0/RAND_MAX;
	mean = 3.6086089;
	ifstream userfile(std::string(DATAPATH)+"/userfile.dta");
	ifstream moviefile(std::string(DATAPATH)+"/moviefile.dta");
	for(int &i : uuser) userfile >>i;
	for(int &i : umovie) moviefile >>i;
	userfile.close();
	moviefile.close();
	_lr = 0.007;//set basic learning rate to be 0.01
	_alpha = 4.91299;
	_lambda = 0.0903385;

}
void bkmodel::half_lr()
{
	_lr *= 0.5;
}

float bkmodel::g(int iu, int im, int it, int ife,float tt)
{
	float gv = 0.0;
	vector<float> tmp(NLAT,0);
	for(int i = 0 ; i < NLAT;i++) tmp[i] = pu[iu][i]+pu1[iu][i]*u(tt);
	gv = mean + bu[iu] + bm[im] + dotprod(tmp,pm[im])  + btm[im][it]*btu[iu] + bt[it] + bf[ife]+bfm[im][ife];
	return gv;
}

void bkmodel::update_param_sgd(feature &a)
{
	a.retrieve_feature();
	int iu = a._iu -1;
	int im = a._im -1;
	int it = a._it;
	int ife = a._if;
	int rate = a._rate;
	float tt = a._tb;
	float gamma = 0.0785714;
	float error = -g(iu,im,it,ife,tt) + rate ;
	float tmp1,tmp2,tmp3;

	bu[iu] -= _lr*(-2.0*error); //no need to regulating thses terms
	bm[im] -= _lr*(-2.0*error);

	for(int i = 0; i < NLAT; i++)
	{
		
		tmp1= pm[im][i]-_lr*(-2.0*error*(pu[iu][i]+pu1[iu][i]*u(tt))) ;
		tmp2= pu[iu][i]-_lr*(-2.0*error*(pm[im][i]));
	    tmp3= pu1[iu][i]-0.001*_lr*(-2.0*error*u(tt)*pm[im][i]);//		}
	
		tmp1 *=(1-0.00001);//weight decay
		tmp2 *=(1-0.00001);
		tmp3 *=(1-0.00001);
		pm[im][i] = tmp1;
		pu[iu][i] = tmp2;
		pu1[iu][i] = tmp3;
	}
	tmp1 = btm[im][it] - _lr*(-2.0*error*btu[iu]);
	tmp2 = btu[iu] - _lr*(-2.0*error*btm[im][it]);
	tmp1 *= (1-0.00001);
	tmp2 *= (1-0.00001);
	btu[iu] = tmp2;
	btm[im][it] = tmp1;
	bt[it] -= _lr*(-2.0*error);
	bf[ife] -= _lr*(-2.0*error);
	bfm[im][ife] -= _lr*(-2.0*error);
	mean -= 0.001*_lr*(-2.0*error);
}




