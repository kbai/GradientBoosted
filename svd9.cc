#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>
#define NUSER 458293
#define NMOVIE 17770
#define NF 200 // number of latent factors
using namespace std;

/*
#iteration = 68
rmse on valid set = 0.798261
rmse on probe set = 0.902711
output output_svd(f=200).dta generated.
*/


default_random_engine e(time(0));

double operator*(vector<double>& v1, vector<double>&v2){
	if(v1.size() != v2.size()){
		cout << "vector dimension error." << endl;
		exit(0);
	}
	double res = 0;
	for(int i = 0; i < v1.size(); i++){
		res += v1[i]*v2[i];
	}
	return res;
}

struct Feature{
    int u;
    int m;
    int t;
    int r;
    Feature(int uu, int mm, int tt, int rr) : u(uu), m(mm), t(tt), r(rr) {}
};

void load_all_data(vector<Feature> &v1, vector<Feature> &v2, vector<Feature> &v3, vector<Feature> &v4, vector<Feature> &v5){
    ifstream all_dta("all.dta");
    ifstream all_idx("all.idx");
    int index;
    int uu, mm, tt, rr;
    while(all_idx >> index){
        all_dta >> uu >> mm >> tt >> rr;
        if(index == 1) v1.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 2) v2.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 3) v3.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 4) v4.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 5) v5.push_back(Feature(uu-1, mm-1, tt-1, rr));
    }
	cout << v1.size() << endl;
    cout << v2.size() << endl;
    cout << v3.size() << endl;
    cout << v4.size() << endl;
    cout << v5.size() << endl;

	all_dta.close();
	all_idx.close();
}

class SVD{
private:
	double mu = 3.6095162;
	vector<double> bu = vector<double>(NUSER, 0), bm = vector<double>(NMOVIE, 0);
	vector<vector<double>> p= vector<vector<double>>(NUSER, vector<double>(NF,0)), q= vector<vector<double>>(NMOVIE, vector<double>(NF,0));
	double eta_pq = 0.002, lambda_pq = 0.015, eta_b = 0.002, lambda_b = 0.008;

public:
    SVD(){
        uniform_real_distribution<double> u(0,1);
        //for(double &x: bu) x = u(e);
		//for(double &x: bm) x = u(e);
        for(vector<double> &v : p) for(double &x : v) x = u(e) / NF;
        for(vector<double> &v : q) for(double &x : v) x = u(e) / NF;
    }

	double r_hat(int u, int m){
		return mu + bu[u] + bm[m] + p[u]*q[m];
	}
	
	void sgd(Feature &data){	
		int u = data.u;
		int m = data.m;
		int r = data.r;
		double eum = r - r_hat(u, m);
		bu[u] = bu[u] + eta_b * (eum - lambda_b * bu[u]);
		bm[m] = bm[m] + eta_b * (eum - lambda_b * bm[m]);
		for(int l = 0; l < NF; l++){
			p[u][l] = p[u][l] + eta_pq * (eum * q[m][l] - lambda_pq * p[u][l]);
			q[m][l] = q[m][l] + eta_pq * (eum * p[u][l] - lambda_pq * q[m][l]);
		}
	}

    void show_parameters(){
        uniform_int_distribution<int> u_random(0, NUSER-1), m_random(0, NMOVIE - 1), l_random(0, NF-1);
        for(int i = 0; i < 10; i++){
            cout << "bu=" << bu[u_random(e)] << " ,bm=" << bm[m_random(e)] <<" ,p=" << p[u_random(e)][l_random(e)] << " ,q=" << q[m_random(e)][l_random(e)] << endl;
        }
    }

	void learning_rate_decay(){
		eta_b *= 0.90;
		eta_pq *= 0.90;
	}

	double rmse(vector<Feature> &set){
    	int N = set.size();
    	double res = 0;
    	for(Feature data : set){
        	res += pow((r_hat(data.u, data.m) - data.r), 2.0) / N;
		}
		return sqrt(res);
    }
	
	void predict(vector<Feature> &set, string s){
		int N = set.size();
		ofstream file(s);
		for(Feature data : set){
			file << r_hat(data.u, data.m) << endl;
		}
		file.close();
		cout << "output " << s << " generated." << endl;
	}

};


int main(){

	vector<Feature> train, valid, hidden, probe, qual;
    train.reserve(94362233);
    valid.reserve(1965045);
    hidden.reserve(1964391);
    probe.reserve(1374739);
    qual.reserve(2749898);
    load_all_data(train, valid, hidden, probe, qual);

	float testscore;
	float mintestscore = 0.95;
	
	SVD model;

	for(int i = 0; i < 250; i++){
		cout << "#iteration = " << i << endl;
		for(int j = 0; j < 94362233; j++){ 
			model.sgd(train[j]);
		}
		model.show_parameters();
		if(i > 55) model.learning_rate_decay();
		cout << "rmse on valid set = " << model.rmse(valid) << endl;
		testscore = model.rmse(probe);
		cout << "rmse on probe set = " << testscore << endl;
		if(testscore  < mintestscore){
			mintestscore = testscore;
			model.predict(qual, "output_svd(f=200).dta");
		}
	}
}

