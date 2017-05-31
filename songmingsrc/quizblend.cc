#include<iostream>
#include<random>
#include<ctime>
#include<vector>
#include<string>
#include<fstream>
#define N 2749898
using namespace std;

// LUP decomposition: A -> PA=LU
void LUP_decompose(vector<vector<double>> &a, vector<vector<double>> &l, vector<vector<double>> &u, vector<int> &pi){
	int n = a.size();
    for(int k = 0; k < n - 1; k++){
    //exchange rows in case of zero pivot
        if(a[k][k] == 0){
            int kk = k + 1;
            while (a[kk][k] == 0) kk++;
            if(kk == n) std::cout << "singular matrix" << std::endl;
            else{
                swap(pi[k], pi[kk]);
                for(int i = 0; i < n; i++) swap(a[k][i], a[kk][i]);
            }
        }
    //standard LU decomposition
        for(int i = k + 1; i < n; i++){
            a[i][k] /= a[k][k];
            for(int j = k + 1; j < n; j++){
                a[i][j] -= a[i][k] * a[k][j];
            }
        }
    }
    //set the value of matrices L and U
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i < j){
                u[i][j] = a[i][j];
                l[i][j] = 0;
            }
            if(i == j){
                l[i][j] = 1;
                u[i][j] = a[i][j];
            }
            if(i > j){
                u[i][j] = 0;
                l[i][j] = a[i][j];
            }
        }
    }
}

//solve Ax=b -> PAx=Pb -> LUx=Pb -> Ly=Pb && Ux=y
void linear_system_solver(vector<vector<double>> &a, vector<double> &x, vector<double> &b){

	int n = a.size();
    vector<vector<double>> l(n, vector<double>(n, 0)), u(n, vector<double>(n, 0)); 
	vector<double> y(n, 0);
    vector<int> pi(n, 0);
    for(int i = 0; i < n; i++) pi[i] = i;
    LUP_decompose(a,l,u,pi);
    for(int i = 0; i < n; i++){
        y[i] = b[pi[i]];
        for(int j = 0; j < i; j++){
            y[i] -= l[i][j]*y[j];
        }
    }
    for(int i = n-1; i >= 0; i--){
        x[i] = y[i];
        for(int j = n-1; j > i; j--){
            x[i] -= u[i][j]*x[j];
        }
        x[i] /= u[i][i];
    }
}

void initialize_XT(vector<string> &predictors, vector<double> &rmse, vector<vector<double>> &XT){
	const int P = predictors.size();
	for(int i = 1; i <= P; i++){
		ifstream datafile(predictors[i-1]);
		int u = 0;
		double r;
		while(datafile >> XT[i][u++]);
		cout << rmse[i-1] << " " << predictors[i-1] << " loaded.";
		cout << " " << u << endl;
		datafile.close();
	}
	cout << "XT initialized." << endl;
}

void initialize_XTX(vector<vector<double>> &XT, vector<vector<double>> &XTX){
	double lambda = 0.0014;
	for(int i = 0; i < XT.size(); i++){
		for(int j = 0; j < XT.size(); j++){
			double sum = 0;
			for(int u = 0; u < N; u++) sum += XT[i][u] * XT[j][u];
			XTX[i][j] = sum;
		}
	}
	for(int i = 0; i < XT.size(); i++){
		XTX[i][i] += lambda * N;
	}
	cout << "XTX initialized." << endl;
}

void initialize_XTy(vector<double> &rmse, vector<vector<double>> &XT, vector<double> &XTy){
	XTy[0] = N * 3.674;
	for(int j = 1; j < XTy.size(); j++){
		double sum = 0;
		for(int u = 0; u < N; u++){
			sum += XT[j][u] * XT[j][u];
		}
		XTy[j] = 0.5 * (N * (pow(3.674, 2) + 1.2699) + sum - N * pow(rmse[j-1], 2));
	}
	cout << "XTy initialized." << endl;
}

void predict(vector<vector<double>> &XT, vector<double> beta, string s){
	ofstream file(s);
	for(int u = 0; u < N; u++){
		double y = 0;
		for(int i = 0; i < XT.size(); i++){
			y += beta[i] * XT[i][u];
		}
		if(y > 5) y = 5;
		if(y < 1) y = 1;
		file << y << endl;
	}
	file.close();
	cout << "output " << s << " generated." << endl;
}
		

int main(){
	vector<string> predictors = {"predictor86483.dta", "predictor86498.dta", "predictor86547.dta", "predictor86552.dta", "predictor86566.dta", "predictor86601.dta", "predictor86607.dta", "predictor86655.dta", "predictor86664.dta", "predictor86676.dta", "predictor86699.dta", "predictor86709.dta", "predictor86721.dta", "predictor86744.dta", "predictor86762.dta", "predictor86851.dta", "predictor86895.dta", "predictor87133.dta", "predictor87416.dta", "predictor87815.dta", "predictor87916.dta", "predictor88357.dta", "predictor88638.dta", "predictor88644.dta", "predictor88863.dta", "predictor88927.dta", "predictor88929.dta", "predictor89014.dta", "predictor89024.dta", "predictor89037.dta", "predictor89041.dta", "predictor89051.dta", "predictor89055.dta", "predictor89155.dta", "predictor89162.dta", "predictor89223.dta", "predictor89232.dta", "predictor89516.dta", "predictor89524.dta", "predictor89802.dta", "predictor90214.dta", "predictor90219.dta", "predictor90635.dta", "predictor92643.dta", "predictor92908.dta", "predictor95337.dta", "predictor97970.dta"};
	vector<double> rmse = {0.86483, 0.86498, 0.86547, 0.86552, 0.86566, 0.86601, 0.86607, 0.86655, 0.86664, 0.86676, 0.86699, 0.86709, 0.86721, 0.86744, 0.86762, 0.86851, 0.86895, 0.87133, 0.87416, 0.87815, 0.87916, 0.88357, 0.88638, 0.88644, 0.88863, 0.88927, 0.88929, 0.89014, 0.89024, 0.89037, 0.89041, 0.89051, 0.89055, 0.89155, 0.89162, 0.89223, 0.89232, 0.89516, 0.89524, 0.89802, 0.90214, 0.90219, 0.90635, 0.92643, 0.92908, 0.95337, 0.97970};
	if(predictors.size() != rmse.size()) exit(0);
	const int P = predictors.size();
	vector<vector<double>> XT = vector<vector<double>>(P+1, vector<double>(N,1));
	vector<vector<double>> XTX = vector<vector<double>>(P+1, vector<double>(P+1));
	vector<double> XTy = vector<double>(P+1);
	vector<double> beta = vector<double>(P+1, 0);
	initialize_XT(predictors, rmse, XT);
	initialize_XTX(XT, XTX);
	initialize_XTy(rmse, XT, XTy);
	linear_system_solver(XTX, beta, XTy);
	cout << "beta solved." << endl;
	for(double x : beta) cout << x << " ";
	cout << endl;
	predict(XT, beta, "output_quizblend.dta");	


}


