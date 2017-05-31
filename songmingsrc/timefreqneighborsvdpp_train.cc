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
#define NF 250 // number of latent factors
using namespace std;


/***************************
This program trains a time and frequency dependent SVD model, which is a combination of Eq.(11), (15) and (16) in [Y. Koren, 2009].
Notice that the y_ij and p_ukt terms are not included.
****************************/

/*
Trained on set 1 to 4.
#iteration = 15
rmse on quiz set = 0.88863
output output_timefrequencysvd6(15).dta generated.
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

void load_all_data(vector<Feature> &v1, vector<Feature> &v2, vector<Feature> &v3, vector<Feature> &v4, vector<Feature> &v5, vector<float> &Tu, vector<vector<int>> &Fut, vector<vector<Feature>> &Ru, vector<int> &Rm, vector<vector<Feature>> &Nu){
    ifstream all_dta("all.dta");
    ifstream all_idx("all.idx");
    int index;
    int uu, mm, tt, rr;
	int last_uu = 1, n_uu = 0, sum_t;
	while(all_idx >> index){
        all_dta >> uu >> mm >> tt >> rr;
		Nu[uu-1].push_back(Feature(uu-1, mm-1, tt-1, rr));
		//v1.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 1 || index == 2 || index == 3){ 
			//v1.push_back(Feature(uu-1, mm-1, tt-1, rr));
			Ru[uu-1].push_back(Feature(uu-1, mm-1, tt-1, rr));
			Rm[mm-1]++;
		}
        //if(index == 2) v2.push_back(Feature(uu-1, mm-1, tt-1, rr));
        //if(index == 3) v3.push_back(Feature(uu-1, mm-1, tt-1, rr));
        if(index == 4) v4.push_back(Feature(uu-1, mm-1, tt-1, rr));
        //if(index == 5) v5.push_back(Feature(uu-1, mm-1, tt-1, rr));
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

class TimeFreqNeighborSVDpp{
private:
	float mu = 3.6095162;
	vector<float> bu = vector<float>(NUSER, 0);
	vector<float> bm = vector<float>(NMOVIE, 0);
	vector<float> au = vector<float>(NUSER, 0);
	vector<float> cu = vector<float>(NUSER, 1);
	vector<vector<float>> but = vector<vector<float>>(NUSER, vector<float>(NTIME, 0));
	vector<vector<float>> bmbin = vector<vector<float>>(NMOVIE, vector<float>(30, 0)); 
	vector<vector<float>> cut = vector<vector<float>>(NUSER, vector<float>(NTIME, 0));
	vector<vector<float>> bmf = vector<vector<float>>(NMOVIE, vector<float>(5, 0));
	vector<float> &Tu;
	vector<vector<int>> &Fut;
	vector<vector<Feature>> &Ru;
	vector<int> &Rm;
	vector<vector<Feature>> &Nu;

    vector<vector<float>> wmn = vector<vector<float>>(NMOVIE, vector<float>(NMOVIE, 0));
    vector<vector<float>> cmn = vector<vector<float>>(NMOVIE, vector<float>(NMOVIE, 0));
    vector<float> bu_tilde = vector<float>(NUSER, 0);
    vector<float> bm_tilde = vector<float>(NMOVIE, 0);

	vector<vector<float>> pu = vector<vector<float>>(NUSER, vector<float>(NF,0));
	vector<vector<float>> alphau = vector<vector<float>>(NUSER, vector<float>(NF, 0));
	vector<vector<float>> qm = vector<vector<float>>(NMOVIE, vector<float>(NF, 0));
	vector<vector<vector<float>>> qmf = vector<vector<vector<float>>>(NMOVIE, vector<vector<float>>(5, vector<float>(NF, 0)));
	vector<vector<float>> y = vector<vector<float>>(NMOVIE, vector<float>(NF,0));
	vector<vector<float>> Y = vector<vector<float>>(NUSER, vector<float>(NF,0));
	
	float eta_bu = 0.00267, lambda_bu = 0.0255;
	float eta_but = 0.00257, lambda_but = 0.00231;
	float eta_au = 3.11E-6, lambda_au = 3.95;
	float eta_bm = 0.000488, lambda_bm = 0.0255;
	float eta_bmbin = 0.000115, lambda_bmbin = 0.0929;
	float eta_cu = 0.00564, lambda_cu = 0.0476;
	float eta_cut = 0.00103, lambda_cut = 0.019;
	float eta_bmf = 0.00236, lambda_bmf = 1.1E-8;

	float lambda1 = 25, lambda2 = 10;
	float eta_wmn = 0.001, lambda_wmn = 0.015;
	float eta_cmn = 0.001, lambda_cmn = 0.015;
	//float eta_wmn = 0.005, lambda_wmn = 0.002;
	//float eta_cmn = 0.005, lambda_cmn = 0.002;

	float eta_pu = 0.008, lambda_pu = 0.015;
	float eta_alphau = 1E-5, lambda_alphau = 50;
	float eta_qm = 0.008, lambda_qm = 0.015;
	float eta_qmf = 2E-5, lambda_qmf = 0.02;	
	float eta_y = 0.008, lambda_y = 0.015;
	
public:
    TimeFreqNeighborSVDpp(vector<float> &Tuu, vector<vector<int>> &Futt, vector<vector<Feature>> &Ruu, vector<int> &Rmm, vector<vector<Feature>> &Nuu) : Tu(Tuu), Fut(Futt), Ru(Ruu), Rm(Rmm), Nu(Nuu) {
        uniform_real_distribution<float> u(-0.005,0.005);
        for(vector<float> &v : pu) for(float &x : v) x = u(e);
        for(vector<float> &v : qm) for(float &x : v) x = u(e);
    }


	float dev(int u, int t){
		return t - Tu[u] > 0 ? pow(t-Tu[u], 0.4) : - pow(Tu[u]-t, 0.4);
	}
	
	int bin(int t){
		return t / 75;
	}

	int freq(int u, int t){
	/*
		float a = 6.76;
		if(Fut[u][t] == 0) return 0;
		return log(Fut[u][t]) / log(a);
	*/
		if(Fut[u][t] <= 6) return 0;
		if(Fut[u][t] <= 45) return 1;
		if(Fut[u][t] <= 308) return 2;
		if(Fut[u][t] <= 2088) return 3;
		else return 4;
	}

    void initialize_btilde(){
        for(int u = 0; u < NUSER; u++){
            for(Feature &data : Ru[u]){
                bm_tilde[data.m] += (float)(data.r - mu) / (lambda1 + Rm[data.m]);
            }
        }
        for(int u = 0; u < NUSER; u++){
            for(Feature &data : Ru[u]){
                bu_tilde[data.u] += (float)(data.r - mu - bm_tilde[data.m]) / (lambda2 + Ru[data.u].size());
            }
        }
    }

    float bum_tilde(int u, int m){
        return mu + bu_tilde[u] + bm_tilde[m];
    }

	float r_hat(int u, int m, int t){
		float pqy_term = 0;
		float wc_term = 0;
        float sqrtRu = sqrt(Ru[u].size()), sqrtNu = sqrt(Nu[u].size());
        for(Feature x : Ru[u]){
            int n = x.m;
            wc_term += (x.r - bum_tilde(u,n)) * wmn[m][n] / sqrtRu;
        }
        for(Feature x : Nu[u]){
            int n = x.m;
            wc_term += cmn[m][n] / sqrtNu;
        }

		for(int l = 0; l < NF; l++){
			pqy_term +=	(pu[u][l] + alphau[u][l]*dev(u,t) + Y[u][l]) * (qm[m][l] + qmf[m][freq(u,t)][l]);
		}	

		return mu + bu[u] + au[u]*dev(u,t) + but[u][t] + (bm[m] + bmbin[m][bin(t)]) * (cu[u] + cut[u][t]) + bmf[m][freq(u,t)] + wc_term + pqy_term;
	}
	
	void sgd(){	
		for(int u = 0; u < NUSER; u++){
			float sqrtRu = sqrt(Ru[u].size()), sqrtNu = sqrt(Nu[u].size());
			for(Feature &data : Nu[u]){
				for(int l = 0; l < NF; l++){
					Y[u][l] += y[data.m][l] / sqrtNu;
				}
			}
			vector<float> cached_Y = Y[u];
			for(Feature &data : Ru[u]){
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
				bmf[m][freq(u,t)] = bmf[m][freq(u,t)] + eta_bmf * (eum - lambda_bmf * bmf[m][freq(u,t)]);
				for(int l = 0; l < NF; l++){
					float cached_pu = pu[u][l], cached_qm = qm[m][l], cached_alphau = alphau[u][l], cached_qmf = qmf[m][freq(u,t)][l];
					pu[u][l] = cached_pu + eta_pu * (eum * (cached_qm + cached_qmf) - lambda_pu * cached_pu);
					qm[m][l] = cached_qm + eta_qm * (eum * (cached_pu + cached_alphau*dev(u,t) + Y[u][l]) - lambda_qm * cached_qm);
					alphau[u][l] = cached_alphau + eta_alphau * (eum * (cached_qm + cached_qmf) * dev(u,t) - lambda_alphau * cached_alphau);
					qmf[m][freq(u,t)][l] = cached_qmf + eta_qmf * (eum * (cached_pu + cached_alphau*dev(u,t) + Y[u][l]) - lambda_qmf * cached_qmf);
					Y[u][l] = Y[u][l] + eta_y * (eum * (cached_qm + cached_qmf) - lambda_y * Y[u][l]);
				}

                for(Feature &x : Ru[u]){
                    int n = x.m;
                    if(n != m){
                        wmn[m][n] = wmn[m][n] + eta_wmn * (eum * (x.r - bum_tilde(u, n)) / sqrtRu - lambda_wmn * wmn[m][n]);
                    }
                }
                for(Feature &x : Nu[u]){
                    int n = x.m;
                    if(n != m){
                        cmn[m][n] = cmn[m][n] + eta_cmn * (eum / sqrtNu - lambda_cmn * cmn[m][n]);
                    }
                }				
			}
			for(Feature &data : Nu[u]){
				int n = data.m;
				for(int l = 0; l < NF; l++){
					y[n][l] = y[n][l] + (Y[u][l]-cached_Y[l])/sqrtNu + eta_y*lambda_y* Ru[u].size()*(cached_Y[l]/sqrtNu - y[n][l]);
				}
			}				
		}
	}

    void show_parameters(){
        uniform_int_distribution<int> u_random(0, NUSER-1), m_random(0, NMOVIE-1), t_random(0, NTIME-1), f_random(0, 4), l_random(0, NF-1);
        for(int i = 0; i < 10; i++){
            cout << "bu=" << bu[u_random(e)] << " ,au=" << au[u_random(e)] <<" ,but=" << but[u_random(e)][t_random(e)] << " ,bm=" << bm[m_random(e)];
			cout << " ,bmbin=" << bmbin[m_random(e)][bin(t_random(e))] << " ,cu=" << cu[u_random(e)] << " ,cut=" << cut[u_random(e)][t_random(e)];
			cout << " ,bmf=" << bmf[m_random(e)][f_random(e)] << endl;
			cout << "pu=" << pu[u_random(e)][l_random(e)] << " ,alphau=" << alphau[u_random(e)][l_random(e)];
			cout << " ,qm=" << qm[m_random(e)][l_random(e)] << " ,qmf=" << qmf[m_random(e)][f_random(e)][l_random(e)];
			cout << " ,y=" << y[m_random(e)][l_random(e)] << " ,Y=" << Y[u_random(e)][l_random(e)];
			cout << " ,bu_tilde=" << bu[u_random(e)] << " ,bm_tilde=" << bm[m_random(e)] << " ,wmn=" << wmn[m_random(e)][m_random(e)] << " ,cmn=" << cmn[m_random(e)][m_random(e)] << endl;
        }
    }

	void learning_rate_decay(){
        eta_bu *= 0.90;
        eta_but *= 0.90;
        eta_au *= 0.90;
        eta_bm *= 0.90;
        eta_bmbin *= 0.90;
        eta_cu *= 0.90;
        eta_cut *= 0.90;
        eta_bmf *= 0.90;

        eta_pu *= 0.90;
        eta_qm *= 0.90;
        eta_y *= 0.90;
        eta_alphau *= 0.90;
        eta_qmf *= 0.90;
		
		eta_wmn *= 0.90;
		eta_cmn *= 0.90;

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
	vector<vector<int>> Fut(NUSER, vector<int>(NTIME));
	vector<vector<Feature>> Ru(NUSER), Nu(NUSER);
	vector<int> Rm(NMOVIE);
    //train.reserve(102416306);
    //valid.reserve(1965045);
    //hidden.reserve(1964391);
    probe.reserve(1374739);
    //qual.reserve(2749898);
    load_all_data(train, valid, hidden, probe, qual, Tu, Fut, Ru, Rm, Nu);


	float testscore;
	float mintestscore = 0.893;
	
	TimeFreqNeighborSVDpp model(Tu, Fut, Ru, Rm, Nu);
	model.initialize_btilde();

	for(int i = 0; i < 100; i++){
		cout << "#iteration = " << i << endl;
		model.sgd();
		model.show_parameters();
		model.learning_rate_decay();
		//cout << "rmse on valid set = " << model.rmse(valid) << endl;
		testscore = model.rmse(probe);
		cout << "rmse on probe set = " << testscore << endl;
		if(testscore  < mintestscore){
			mintestscore = testscore;
			model.predict(probe, "output_timefreqneighborsvdpp_train.dta");
		}
	}
}

