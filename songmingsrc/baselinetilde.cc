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
#include <unordered_map>
#define NUSER 458293
#define NMOVIE 17770
#define NTIME 2243
#define NF 200 // number of latent factors
using namespace std;


/***************************
This program trains a simple baseline predictor without SGD, which combines Eq.(1), (2) and (3) of [Y. Koren, 2009].
It is also a component of several neighborhood models.
****************************/

/*
rmse on train set = 0.917583
rmse on valid set = 0.921833
rmse on hidden set = 0.921225
rmse on probe set = 0.981902
*/

default_random_engine e(time(0));

struct Feature{
    int u;
    int m;
    int t;
    int r;
    Feature(int uu, int mm, int tt, int rr) : u(uu), m(mm), t(tt), r(rr) {}
};

void load_all_data(vector<Feature> &v1, vector<Feature> &v2, vector<Feature> &v3, vector<Feature> &v4, vector<Feature> &v5, vector<double> &Tu, vector<unordered_map<int, int>> &Fut, vector<vector<int>> &Ru, vector<vector<int>> &Rm){
    ifstream all_dta("all.dta");
    ifstream all_idx("all.idx");
    int index;
    int uu, mm, tt, rr;
	int last_uu = 1, n_uu = 0, sum_t;
	while(all_idx >> index){
        all_dta >> uu >> mm >> tt >> rr;
        if(index == 1) v1.push_back(Feature(uu-1, mm-1, tt-1, rr)), Ru[uu-1].push_back(mm-1), Rm[mm-1].push_back(uu-1);
        if(index == 2) v2.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 3) v3.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 4) v4.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 5) v5.push_back(Feature(uu-1, mm-1, tt-1, rr));
		Fut[uu-1][tt-1]++;
        if(uu == last_uu){
            n_uu++;
            sum_t += tt-1;
        }
        else{
            Tu[last_uu-1] = (double)sum_t / n_uu;
            sum_t = tt-1;
            n_uu = 1;
            last_uu = uu;
        }
    }
	Tu[last_uu-1] = (double)sum_t / n_uu;
	cout << v1.size() << endl;
    cout << v2.size() << endl;
    cout << v3.size() << endl;
    cout << v4.size() << endl;
    cout << v5.size() << endl;
	all_dta.close();
	all_idx.close();
}

class BaselineTilde{
private:
	double mu = 3.6095162;
	vector<double> bu_tilde = vector<double>(NUSER, 0);
	vector<double> bm_tilde = vector<double>(NMOVIE, 0);
	vector<vector<int>> &Ru;
	vector<vector<int>> &Rm;

	int lambda1 = 25, lambda2 = 10;
public:
	BaselineTilde(vector<vector<int>> &Ruu, vector<vector<int>> &Rmm) : Ru(Ruu), Rm(Rmm) {}

	void initialize_btilde(vector<Feature> &set){
		for(Feature data : set){
			bm_tilde[data.m] += (double)(data.r - mu) / (lambda1 + Rm[data.m].size());
		}
		for(Feature data : set){
			bu_tilde[data.u] += (double)(data.r - mu - bm_tilde[data.m]) / (lambda2 + Ru[data.u].size());
		}
	}

	double r_hat(int u, int m, int t){
		return mu + bu_tilde[u] + bm_tilde[m];
	}
	
	double rmse(vector<Feature> &set){
    	int n = set.size();
    	double res = 0;
    	for(Feature data : set){
        	res += pow((r_hat(data.u, data.m, data.t) - data.r), 2.0) / n;
		}
		return sqrt(res);
    }
	
	void predict(vector<Feature> &set, string s){
		int n = set.size();
		ofstream file(s);
		for(Feature data : set){
			file << r_hat(data.u, data.m, data.t) << endl;
		}
		file.close();
		cout << "output " << s << " generated." << endl;
	}
    void show_parameters(){
        uniform_int_distribution<int> u_random(0, NUSER-1), m_random(0, NMOVIE-1);
        for(int i = 0; i < 50; i++){
			int u = u_random(e), m = m_random(e);
            cout << "u=" << u << " ,bu_tilde=" << bu_tilde[u] << " ,m=" << m << " ,bm_tilde=" << bm_tilde[m] << endl;;
       	}
    }
};



int main(){

	vector<Feature> train, valid, hidden, probe, qual;
	vector<double> Tu(NUSER,0);
	vector<unordered_map<int, int>> Fut(NUSER);
	vector<vector<int>> Ru(NUSER);
	vector<vector<int>> Rm(NMOVIE);
    train.reserve(94362233);
    valid.reserve(1965045);
    hidden.reserve(1964391);
    probe.reserve(1374739);
    qual.reserve(2749898);
    load_all_data(train, valid, hidden, probe, qual, Tu, Fut, Ru, Rm);

	BaselineTilde model(Ru, Rm);
	model.initialize_btilde(train);
	model.show_parameters();
	cout << "rmse on train set = " << model.rmse(train) << endl;
	cout << "rmse on valid set = " << model.rmse(valid) << endl;
	cout << "rmse on hidden set = " << model.rmse(hidden) << endl;
	cout << "rmse on probe set = " << model.rmse(probe) << endl;
}

