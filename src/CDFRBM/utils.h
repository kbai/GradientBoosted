#pragma once

#include <iostream>
#include <math.h>
#include <assert.h>
using namespace std;


namespace utils {
  
  double uniform(double min, double max) {
    return rand() / (RAND_MAX + 1.0) * (max - min) + min;
  }

  int binomial(int n, double p) {
    if(p < 0 || p > 1) return 0;
  
    int c = 0;
    double r;
  
    for(int i=0; i<n; i++) {
      r = rand() / (RAND_MAX + 1.0);
      if (r < p) c++;
    }

    return c;
  }

  void binomial5(const double*p, int* s) {
  
    double u = 0.0;
	double r;
	s[0] = 0;
	s[1] = 0;
	s[2] = 0;
	s[3] = 0;
	s[4] = 0;
    r = rand() / (RAND_MAX + 1.0);
 
    for(int i=0; i<5; i++) 
	{
		u += p[i];
		if(r < u) 
		{
			s[i] = 1;
			break;
		}	
    }
//	cout << p[0] << p[1] <<p[2]<<p[3]<< p[4]<<u << r << endl;
//	cout << s[0] << s[1] << s[2] <<s[3] << s[4] << endl;
    return ;
  }


  double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
  }

  void sigmoid5(const double* x, double* y) {
	double sum = 0.0;
//	assert((x[0] > -10)&(x[0] < 10));
	for(int ns = 0; ns < 5; ns++) sum += exp(x[ns]);
	for(int ns = 0; ns < 5; ns++) y[ns] = exp(x[ns])/sum;
//	cout << y[0] <<"\t "<< y[1]<<"\t " << y[2]<<"\t " << y[3]<<"\t " << y[4] <<"hehe"<< endl;
//	cout << x[0] <<"\t"<< x[1] <<"\t"<< x[2]<<"\t" << x[3]<<"\t" << x[4] <<"hehe"<< endl;
  }

}
