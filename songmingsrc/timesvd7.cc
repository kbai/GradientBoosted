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
#define NTIME 2243
#define NF 200 // number of latent factors
using namespace std;


/***************************
This program trains a time dependent SVD, which is a combination of Eq.(10), (13) and (15) in [Y. Koren, 2009]. 
Notice that the y_j and p_ukt terms are omitted.
****************************/

/*
Train on all data from 1 to 4
#iteration = 15
rmse on quiz set = 0.89024
output output_timesvd7(15).dta generated.
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

void load_all_data(vector<Feature> &v1, vector<Feature> &v2, vector<Feature> &v3, vector<Feature> &v4, vector<Feature> &v5, vector<float> &Tu){
    ifstream all_dta("all.dta");
    ifstream all_idx("all.idx");
    int index;
    int uu, mm, tt, rr;
	int last_uu = 1, n_uu = 0, sum_t;
    while(all_idx >> index){
        all_dta >> uu >> mm >> tt >> rr;
        if(index != 5) v1.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 2) v2.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 3) v3.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 4) v4.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 5) v5.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(uu == last_uu){
            n_uu++;
            sum_t += tt-1;
        }
        else{
            Tu[last_uu-1] = (float)sum_t / n_uu;
            sum_t = tt-1;
            n_uu = 1;
            last_uu = uu;
        }
    }
	Tu[last_uu-1] = (float)sum_t / n_uu;
	cout << v1.size() << endl;
    cout << v2.size() << endl;
    cout << v3.size() << endl;
    cout << v4.size() << endl;
    cout << v5.size() << endl;

	all_dta.close();
	all_idx.close();
}

class TimeSVD{
private:
	float mu = 3.6095162;
	vector<float> bu = vector<float>(NUSER, 0);
	vector<float> bm = vector<float>(NMOVIE, 0);
	vector<float> au = vector<float>(NUSER, 0);
	vector<float> cu = vector<float>(NUSER, 1);
	vector<vector<float>> but = vector<vector<float>>(NUSER, vector<float>(NTIME, 0));
	vector<vector<float>> bmbin = vector<vector<float>>(NMOVIE, vector<float>(30, 0)); 
	vector<vector<float>> cut = vector<vector<float>>(NUSER, vector<float>(NTIME, 0));
	vector<float> &Tu;
	
	vector<vector<float>> pu = vector<vector<float>>(NUSER, vector<float>(NF,0));
	vector<vector<float>> alphau = vector<vector<float>>(NUSER, vector<float>(NF,0));
	vector<vector<float>> qm = vector<vector<float>>(NMOVIE, vector<float>(NF,0));	
	
	//vector<vector<float>> p= vector<vector<float>>(NUSER, vector<float>(NF,0)), q= vector<vector<float>>(NMOVIE, vector<float>(NF,0));
	//float eta_pq = 0.002, lambda_pq = 0.015, eta_b = 0.002, lambda_b = 0.008;
	
	float eta_bu = 0.003, lambda_bu = 0.03;
	float eta_but = 0.0025, lambda_but = 0.005;
	float eta_au = 0.00001, lambda_au = 50;
	float eta_bm = 0.002, lambda_bm = 0.03;
	float eta_bmbin = 0.00005, lambda_bmbin = 0.1;
	float eta_cu = 0.008, lambda_cu = 0.01;
	float eta_cut = 0.002, lambda_cut = 0.005;

	float eta_pu = 0.008, lambda_pu = 0.015;
	float eta_alphau = 1E-5,  lambda_alphau = 50;
	float eta_qm = 0.008, lambda_qm = 0.015;

public:
    TimeSVD(vector<float> &Tuu) : Tu(Tuu) {
        uniform_real_distribution<float> u(0,1);
        for(vector<float> &v : pu) for(float &x : v) x = u(e) / NF;
        for(vector<float> &v : qm) for(float &x : v) x = u(e) / NF;
    }


	float dev(int u, int t){
		return t - Tu[u] > 0 ? pow(t-Tu[u], 0.4) : - pow(Tu[u]-t, 0.4);
	}
	
	int bin(int t){
		return t / 75;
	}

	float r_hat(int u, int m, int t){
		float pq_term = 0;
		for(int l = 0; l < NF; l++){
			pq_term += (pu[u][l] + alphau[u][l] * dev(u,t)) * qm[m][l];
		}
		return mu + bu[u] + au[u]*dev(u,t) + but[u][t] + (bm[m] + bmbin[m][bin(t)]) * (cu[u] + cut[u][t]) + pq_term;
	}
	
	void sgd(Feature &data){	
		int u = data.u;
		int m = data.m;
		int t = data.t;
		int r = data.r;
		float eum = r - r_hat(u, m, t);
		bu[u] = bu[u] + eta_bu * (eum - lambda_bu * bu[u]);
		au[u] = au[u] + eta_au * (eum * dev(u,t) - lambda_au * au[u]);
		but[u][t] = but[u][t] + eta_but * (eum - lambda_but * but[u][t]);
		bm[m] = bm[m] + eta_bm * (eum * (cu[u] + cut[u][t]) - lambda_bm * bm[m]);
		bmbin[m][bin(t)] = bmbin[m][bin(t)] + eta_bmbin * (eum * (cu[u] + cut[u][t]) -lambda_bmbin * bmbin[m][bin(t)]);
		cu[u] = cu[u] + eta_cu * (eum * (bm[m] + bmbin[m][bin(t)]) - lambda_cu * (cu[u] - 1));
		cut[u][t] = cut[u][t] + eta_cut * (eum * (bm[m] + bmbin[m][bin(t)]) - lambda_cut * cut[u][t]);
		for(int l = 0; l < NF; l++){
			pu[u][l] = pu[u][l] + eta_pu * (eum * qm[m][l] - lambda_pu * pu[u][l]);
			qm[m][l] = qm[m][l] + eta_qm * (eum * (pu[u][l] + alphau[u][l] * dev(u,t)) - lambda_qm * qm[m][l]);
			alphau[u][l] = alphau[u][l] + eta_alphau * (eum * qm[m][l] * dev(u,t) - lambda_alphau * alphau[u][l]);
		}
	}

    void show_parameters(){
        uniform_int_distribution<int> u_random(0, NUSER-1), m_random(0, NMOVIE-1), t_random(0, NTIME-1), l_random(0, NF-1);
        for(int i = 0; i < 10; i++){
            cout << "bu=" << bu[u_random(e)] << " ,au=" << au[u_random(e)] <<" ,but=" << but[u_random(e)][t_random(e)] << " ,bm=" << bm[m_random(e)];
			cout << " ,bmbin=" << bmbin[m_random(e)][bin(t_random(e))] << " ,cu=" << cu[u_random(e)] << " ,cut=" << cut[u_random(e)][t_random(e)];
			cout << " ,pu=" << pu[u_random(e)][l_random(e)] << " ,qm=" << qm[m_random(e)][l_random(e)] << " ,alphau=" << alphau[u_random(e)][l_random(e)] << endl;
        }
    }

	void learning_rate_decay(){
		eta_pu *= 0.90;
		eta_qm *= 0.90;
	}

	float rmse(vector<Feature> &set){
    	int N = set.size();
    	float res = 0;
    	for(Feature data : set){
        	res += pow((r_hat(data.u, data.m, data.t) - data.r), 2.0) / N;
		}
		return sqrt(res);
    }
	
	void predict(vector<Feature> &set, string s){
		int N = set.size();
		ofstream file(s);
		for(Feature data : set){
			file << r_hat(data.u, data.m, data.t) << endl;
		}
		file.close();
		cout << "output " << s << " generated." << endl;
	}

};


int main(){

	vector<Feature> train, valid, hidden, probe, qual;
	vector<float> Tu(NUSER,0);
    train.reserve(99666408);
    valid.reserve(1965045);
    hidden.reserve(1964391);
    probe.reserve(1374739);
    qual.reserve(2749898);
    load_all_data(train, valid, hidden, probe, qual, Tu);

	float testscore;
	float mintestscore = 0.90;
	
	TimeSVD model(Tu);

	for(int i = 0; i < 200; i++){
		cout << "#iteration = " << i << endl;
		for(int j = 0; j < 99666408; j++){ 
			model.sgd(train[j]);
		}
		model.show_parameters();
		//if (i > 10) model.learning_rate_decay();
		cout << "rmse on valid set = " << model.rmse(valid) << endl;
		testscore = model.rmse(probe);
		cout << "rmse on probe set = " << testscore << endl;
		if(testscore  < mintestscore){
			mintestscore = testscore;
			model.predict(qual, "output_timesvd7("+ to_string(i) +").dta");
		}
	}
}
