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
#define NF 50 // number of latent factors
using namespace std;

/*
#iteration = 34
rmse on valid set = 0.812972
rmse on probe set = 0.908595
output output_svdpp(f=50).dta generated.
*/


default_random_engine e(time(0));

float operator*(vector<float>& v1, vector<float>&v2){
	if(v1.size() != v2.size()){
		cout << "vector dimension error." << endl;
		exit(0);
	}
	float res = 0;
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

void load_all_data(vector<Feature> &v1, vector<Feature> &v2, vector<Feature> &v3, vector<Feature> &v4, vector<Feature> &v5, vector<vector<Feature>> &Ru, vector<vector<Feature>> &Nu){
    ifstream all_dta("all.dta");
    ifstream all_idx("all.idx");
    int index;
    int uu, mm, tt, rr;
    while(all_idx >> index){
        all_dta >> uu >> mm >> tt >> rr;
		Nu[uu-1].push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 1){ 
			//v1.push_back(Feature(uu-1, mm-1, tt-1, rr));
			Ru[uu-1].push_back(Feature(uu-1, mm-1, tt-1, rr));
		}
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
	float mu = 3.6095162;
	vector<float> bu = vector<float>(NUSER, 0); 
	vector<float> bm = vector<float>(NMOVIE, 0);
	vector<vector<float>> p= vector<vector<float>>(NUSER, vector<float>(NF,0)); 
	vector<vector<float>> q= vector<vector<float>>(NMOVIE, vector<float>(NF,0));
	vector<vector<float>> y = vector<vector<float>>(NMOVIE, vector<float>(NF,0));
	vector<vector<float>> sum_y = vector<vector<float>>(NUSER, vector<float>(NF,0));
	vector<vector<Feature>> &Ru;
	vector<vector<Feature>> &Nu;

	float eta_b = 0.007, lambda_b = 0.005;
	float eta_pq = 0.007, lambda_pq = 0.015; 
	float eta_y = 0.007, lambda_y = 0.015;

public:
    SVD(vector<vector<Feature>> &Ruu, vector<vector<Feature>> &Nuu) : Ru(Ruu), Nu(Nuu){
        uniform_real_distribution<float> u(0,1);
        for(vector<float> &v : p) for(float &x : v) x = u(e) / NF;
        for(vector<float> &v : q) for(float &x : v) x = u(e) / NF;
    }

	float r_hat(int u, int m){
		float pyq_term = 0;
		for(int l = 0; l < NF; l++){
			pyq_term += (p[u][l] + sum_y[u][l]) * q[m][l];		
		}
		return mu + bu[u] + bm[m] + pyq_term;
	}
	
	void sgd(){
		for(int u = 0; u < NUSER; u++){
			for(Feature data : Nu[u]){
				for(int l = 0; l < NF; l++){
					sum_y[u][l] += y[data.m][l] / sqrt(Nu[u].size());
				}
			}
			vector<float> sum_eq_term = vector<float>(NF, 0);
			for(Feature data : Ru[u]){
				int m = data.m;
				int r = data.r;
				float eum = r - r_hat(u, m);
				bu[u] = bu[u] + eta_b * (eum - lambda_b * bu[u]);
				bm[m] = bm[m] + eta_b * (eum - lambda_b * bm[m]);
				for(int l = 0; l < NF; l++){
					p[u][l] = p[u][l] + eta_pq * (eum * q[m][l] - lambda_pq * p[u][l]);
					q[m][l] = q[m][l] + eta_pq * (eum * (p[u][l] + sum_y[u][l]) - lambda_pq * q[m][l]);
					//y[n][l] = y[n][l] + eta_y * (eum * q[m][l] / sqrt(Nu[u].size()) - lambda_y * y[n][l]);
					sum_eq_term[l] += eum * q[m][l] / sqrt(Nu[u].size());
				}
			}
			for(Feature data : Nu[u]){
				int n = data.m;
				for(int l = 0; l < NF; l++){
					y[n][l] = y[n][l] + eta_y * (sum_eq_term[l] - lambda_y * Ru[u].size() * y[n][l]);
				}
			}
		}
	}

    void show_parameters(){
        uniform_int_distribution<int> u_random(0, NUSER-1), m_random(0, NMOVIE - 1), l_random(0, NF-1);
        for(int i = 0; i < 10; i++){
            cout << "bu=" << bu[u_random(e)] << " ,bm=" << bm[m_random(e)] <<" ,p=" << p[u_random(e)][l_random(e)] << " ,q=" << q[m_random(e)][l_random(e)];
			cout << " ,y=" << y[m_random(e)][l_random(e)] << " ,sum_y=" << sum_y[u_random(e)][l_random(e)] << endl;
        }
    }

	void learning_rate_decay(){
		eta_b *= 0.95;
		eta_pq *= 0.95;
		eta_y *= 0.95;
	}

	float rmse(vector<Feature> &set){
    	int N = set.size();
    	float res = 0;
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
	vector<vector<Feature>> Ru(NUSER), Nu(NUSER);
    //train.reserve(94362233);
    valid.reserve(1965045);
    hidden.reserve(1964391);
    probe.reserve(1374739);
    qual.reserve(2749898);
    load_all_data(train, valid, hidden, probe, qual, Ru, Nu);

	float testscore;
	float mintestscore = 0.95;
	
	SVD model(Ru, Nu);

	for(int i = 0; i < 200; i++){
		cout << "#iteration = " << i << endl;
	
		model.sgd();
	
		model.show_parameters();
		model.learning_rate_decay();
		cout << "rmse on valid set = " << model.rmse(valid) << endl;
		testscore = model.rmse(probe);
		cout << "rmse on probe set = " << testscore << endl;
		if(testscore  < mintestscore){
			mintestscore = testscore;
			model.predict(qual, "output_svdpp(f=50).dta");
		}
	}
}

