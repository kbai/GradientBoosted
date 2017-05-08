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
#define N1 94362233
#define N2 1965045
#define N3 1964391
#define N4 1374739
#define N5 2749898
using namespace std;


/***************************
This program trains BellKor's neighborhood model, which corresponds to Eq.(9) of [Y. Koren, 2008].
****************************/

/*
#iteration = 14
rmse on quiz set = 0.89802
output output_neighborhood_predict.dta generated.
*/



default_random_engine e(time(0));

struct Feature{
    int u;
    int m;
    int t;
    int r;
    Feature(int uu, int mm, int tt, int rr) : u(uu), m(mm), t(tt), r(rr) {}
};

void load_all_data(vector<Feature> &v1, vector<Feature> &v2, vector<Feature> &v3, vector<Feature> &v4, vector<Feature> &v5, vector<float> &Tu, vector<unordered_map<int, int>> &Fut, vector<vector<Feature>> &Ru, vector<int> &Rm, vector<vector<Feature>> &NuminusRu){
    ifstream all_dta("all.dta");
    ifstream all_idx("all.idx");
    int index;
    int uu, mm, tt, rr;
	int last_uu = 1, n_uu = 0, sum_t;
	while(all_idx >> index){
        all_dta >> uu >> mm >> tt >> rr;
        if(index == 1 || index == 2 || index == 3 || index == 4) {
			v1.push_back(Feature(uu-1, mm-1, tt-1, rr)); 
			Ru[uu-1].push_back(Feature(uu-1, mm-1, tt-1, rr)); 
			Rm[mm-1]++;
		}
		else NuminusRu[uu-1].push_back(Feature(uu-1, mm-1, tt-1, rr));
        //if(index == 2) v2.push_back(Feature(uu-1, mm-1, tt-1, rr));
        //if(index == 3) v3.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 4) v4.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 5) v5.push_back(Feature(uu-1, mm-1, tt-1, rr));
		Fut[uu-1][tt-1]++;
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

class Neighborhood{
private:
	float mu = 3.6095162;
	vector<float> bu = vector<float>(NUSER, 0);
	vector<float> bm = vector<float>(NMOVIE, 0);
	vector<vector<float>> wmn = vector<vector<float>>(NMOVIE, vector<float>(NMOVIE, 0));
	vector<vector<float>> cmn = vector<vector<float>>(NMOVIE, vector<float>(NMOVIE, 0));
	vector<float> bu_tilde = vector<float>(NUSER, 0);
	vector<float> bm_tilde = vector<float>(NMOVIE, 0);
	vector<vector<Feature>> &Ru;
	vector<int> &Rm;
	vector<vector<Feature>> &NuminusRu;

	float lambda1 = 25, lambda2 = 10;
	float eta_bu = 0.005, lambda_bu = 0.002;
	float eta_bm = 0.005, lambda_bm = 0.002;
	float eta_wmn = 0.005, lambda_wmn = 0.002;
	float eta_cmn = 0.005, lambda_cmn = 0.002;

public:
	Neighborhood(vector<vector<Feature>> &Ruu, vector<int> &Rmm, vector<vector<Feature>> &NuminusRuu) : Ru(Ruu), Rm(Rmm), NuminusRu(NuminusRuu) {}

	void initialize_btilde(vector<Feature> &set){
		for(Feature data : set){
			bm_tilde[data.m] += (float)(data.r - mu) / (lambda1 + Rm[data.m]);
		}
		for(Feature data : set){
			bu_tilde[data.u] += (float)(data.r - mu - bm_tilde[data.m]) / (lambda2 + Ru[data.u].size());
		}
	}
	
	float bum_tilde(int u, int m){
		return mu + bu_tilde[u] + bm_tilde[m];
	}
	
	float r_hat(int u, int m, int t){
		float wc_term = 0;
		for(Feature x : Ru[u]){
			int n = x.m;
			wc_term += (x.r - bum_tilde(u,n)) * wmn[m][n] / sqrt(Ru[u].size()) + cmn[m][n] / sqrt((Ru[u].size() + NuminusRu[u].size()));
		}
		for(Feature x : NuminusRu[u]){
			int n = x.m;
			wc_term += cmn[m][n] / sqrt((Ru[u].size() + NuminusRu[u].size()));
		}
		return mu + bu[u] + bm[m] + wc_term;
	}

	void sgd(Feature &data){	
		int u = data.u;
		int m = data.m;
		int t = data.t;
		int r = data.r;
		float eum = r - r_hat(u, m, t);
		bu[u] = bu[u] + eta_bu * (eum - lambda_bu * bu[u]);
		bm[m] = bm[m] + eta_bm * (eum - lambda_bm * bm[m]);
		for(Feature x : Ru[u]){
			int n = x.m;
			if(n != m){
				wmn[m][n] = wmn[m][n] + eta_wmn * (eum * (x.r - bum_tilde(u, n)) / sqrt(Ru[u].size()) - lambda_wmn * wmn[m][n]);	
				cmn[m][n] = cmn[m][n] + eta_cmn * (eum / sqrt(Ru[u].size() + NuminusRu[u].size()) - lambda_cmn * cmn[m][n]);
			}
		}
		for(Feature x : NuminusRu[u]){
			int n = x.m;
			cmn[m][n] = cmn[m][n] + eta_cmn * (eum / sqrt(Ru[u].size() + NuminusRu[u].size()) - lambda_cmn * cmn[m][n]);
		}
	}

	void learning_rate_decay(){
		float decayrate = 0.90;
		eta_bu *= decayrate;
		eta_bm *= decayrate;
		eta_cmn *= decayrate;
		eta_wmn *= decayrate;
	}

	float rmse(vector<Feature> &set){
    	int n = set.size();
    	float res = 0;
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
        for(int i = 0; i < 10; i++){
			int u = u_random(e), m = m_random(e);
            cout << "u=" << u << " ,m=" << m << " ,bu_tilde=" << bu_tilde[u] << " ,bm_tilde=" << bm_tilde[m] << " ,bu=" << bu[u] << " ,bm=" << bm[m];
			cout << " ,wmn=" << wmn[m][Ru[u][0].m] << " ,cmn=" << cmn[m][NuminusRu[u][0].m] << endl;
       	}
    }
};



int main(){

	vector<Feature> train, valid, hidden, probe, qual;
	vector<float> Tu(NUSER,0);
	vector<unordered_map<int, int>> Fut(NUSER);
	vector<vector<Feature>> Ru(NUSER);
	vector<int> Rm(NMOVIE);
	vector<vector<Feature>> NuminusRu(NUSER);
    train.reserve(N1+N2+N3+N4);
    //valid.reserve(N2);
    //hidden.reserve(N3);
    probe.reserve(N4);
    qual.reserve(N5);
    load_all_data(train, valid, hidden, probe, qual, Tu, Fut, Ru, Rm, NuminusRu);

    float testscore;
    float mintestscore = 0.94;

	Neighborhood model(Ru, Rm, NuminusRu);
	model.initialize_btilde(train);

    for(int i = 0; i <= 14; i++){
        cout << "#iteration = " << i << endl;
        for(Feature &data : train){
            model.sgd(data);
        }
        model.show_parameters();
        //model.learning_rate_decay();
        //cout << "rmse on valid set = " << model.rmse(valid) << endl;
        testscore = model.rmse(probe);
        cout << "rmse on probe set = " << testscore << endl;
        if(testscore  < mintestscore){
            mintestscore = testscore;
            model.predict(qual, "output_neighborhood_predict.dta");
        }
    }

}

