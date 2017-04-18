#include "commonheader.h"
using namespace std;
#include "feature.hh"
#include "bkmodel.hh"
#include "utils.hh"
#define NLAT 10
#include "bkmodel1.hh"



float bkmodel1::g(int iu, int im, int it)
{
	cout << "bkmodel1" << endl;
	float gv = 0.0;
	gv = mean + bu[iu] + dotprod(bm[im],bs[iu])  + bt[im][it]*btu[iu];
	return gv;
}
