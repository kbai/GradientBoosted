#include "commonheader.h"
using namespace std;
#include "utils.hh"


float dotprod(const vector<float>& v1, const vector<float>& v2)
{
	float sum = 0.0;
	for(int i = 0; i < v1.size(); i++)
		sum += v1[i] * v2[i];
	return sum;
}

