#include "commonheader.h"
using namespace std;
#include "feature.hh"
#include "bkmodel.hh"
#include "utils.hh"
#include "bkmodel1.hh"


bkmodel1::bkmodel1():bkmodel(),bst(NUSER,vector<vector<float>>(30,vector<float>(NLAT,0)))// call superclass constructor
{


}

float bkmodel1::g(int iu, int im, int it)
{
//	cout << "bkmodel1" << endl;
	float gv = 0.0;
	vector<float> tmp(NLAT,0);
	for(int a = 0 ; a < NLAT; a++) tmp[a] = pu[iu][a] + bst[iu][it][a];	
	gv = mean + bu[iu] + bm[im] + dotprod(pm[im],tmp)  + btm[im][it]*btu[iu] + bt[it] ;
	return gv;
}


void bkmodel1::update_param_sgd(feature &a)
{
	a.retrieve_feature();
	int iu = a._iu -1;
	int im = a._im -1;
	int it = a._it;
	int rate = a._rate;
	float error = -g(iu,im,it) + rate ;

	bu[iu] -= _lr*(-2.0*error);
	bm[im] -= _lr*(-2.0*error);

	for(int i = 0; i < NLAT; i++)
	{
//		if(alternating)
//		{
		pm[im][i] -= _lr*(-2.0*error*(pu[iu][i] + bst[iu][it][i]) + 2.0*pm[im][i]*(0.00234879 + 0.0610418/umovie[im])) ;
//		}
//		else
//		{
		pu[iu][i] -= _lr*(-2.0*error*(pm[im][i]) + 2.0*pu[iu][i] *(0.0753211 + 3.11288/uuser[iu]));

		bst[iu][it][i]  -= 0.1*_lr*(-2.0*error*pm[im][i] + 2.0*bst[iu][it][i]*( 0.00753211 + 0.311288/uuser[iu])) ;
//		}


	}
	btm[im][it] -= _lr*(-2.0*error*btu[iu]+ btm[im][it] * 0.0001);
	btu[iu] -= _lr*(-2.0*error*btm[im][it] + btu[iu] * 0.001);
	bt[it] -= _lr*(-2.0*error);
	mean -= _lr*(-2.0*error)/100;
}


