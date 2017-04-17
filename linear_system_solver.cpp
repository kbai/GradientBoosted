#include<iostream>
#include<random>
#include<ctime>
#include<vector>
using namespace std;

/* 
A Linear System Solver based on LUP decomposition

Solve: Ax = b
input: A, n*n  dimensional vector
	   b, n dimensional vector
return: x, n dimensional vector


Example: A = {{3,  2, -1},
             {2,  -2, 4},
	     {-1,0.5,-1}}

         b = {1, -2, 0}

result:  x = {1, -2, -2}

*/


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
vector<double> linear_system_solver(vector<vector<double>> a, vector<double> b){
	try{
		if(a.size() != a[0].size()) throw 1;
		if(a.size() != b.size()) throw 2;
	}
	catch(int e){
		if(e == 1){
			cout << "matrix A should be an n*n matrix" << endl;
			return {};
		}
		if(e == 2){
			cout << "matrix A and vector b should have the same dimension" << endl;
			return {};
		}
	}
	int n = a.size();
	vector<double> x(n, 0);
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
	return x;
}

int main(){
	vector<vector<double>> a({{3,2,-1},{2,-2,4},{-1,0.5,-1}});
	vector<double> b({1,-2,0});
	vector<double> x = linear_system_solver(a, b);
	for(double d : x){
		cout << d << " ";
	}
	cout << endl;
}
