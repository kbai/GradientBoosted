#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <thread>
#include <fstream>
#define NUMTH 50
using namespace std;
default_random_engine e(0);
float sigmoid(float x)
{
	return 1./(1.+exp(-x));
}
class nn
{
	int n_neuron;
	int n_pred;
	vector<vector<float>> w1;
	vector<float> w2;
	vector<float> b1;
	float b2 = 0;
	float lr = 0.01;
	float insamplerror = 5.0;
	public:
	nn(int _neuron, int _pred):
		w1(_neuron,vector<float>(_pred,0)),
		w2(_neuron, 0),
		b1(_pred,0)
	{
		n_neuron = _neuron;
		n_pred = _pred;
		uniform_real_distribution<float> u(-0.005,0.005);
		for(auto&a : w1) for(auto &b : a) b =u(e);
		for(auto&a : w2) a = u(e);	
	}
	void output()
	{
		cout << "w1:" <<  w1[e()%n_neuron][e()%n_pred] << endl;
		cout << "w2:" << w2[e()%n_pred] << endl;
		cout << "b1:" << b1[e()%n_pred] << endl;
		cout << "b2:" << b2 << endl;
		cout << "insample error" << insamplerror<< endl;
	}
	float ratedecay()
	{
		lr = lr * 0.95;
	}

	float predict(vector<float> & data)
	{
		vector<float> f1(n_neuron,0.0);
		vector<float> g1(n_neuron,0.0);
		float f2 = 0.0;
		float g2 = 0.0;
		for(int i = 0; i < n_neuron; i++)
		{
			for(int j = 0; j < n_pred; j++)
			{
				f1[i] += w1[i][j]*data[j]; 
			}
			g1[i] = sigmoid(f1[i]+b1[i]);
		}
		for(int i = 0; i < n_neuron; i++)
		{
			f2 += g1[i] * w2[i];
		}
		f2 += b2;
		g2 = sigmoid(f2);
		return 4.0 * g2 + 1.0;

	}
	float sgd(vector<float>& data,float label)
	{
		vector<float> f1(n_neuron,0.0);
		vector<float> g1(n_neuron,0.0);
		float f2 = 0.0;
		float g2 = 0.0;
		for(int i = 0; i < n_neuron; i++)
		{
			for(int j = 0; j < n_pred; j++)
			{
				f1[i] += w1[i][j]*data[j]; 
			}
			g1[i] = sigmoid(f1[i]+b1[i]);
		}
		for(int i = 0; i < n_neuron; i++)
		{
			f2 += g1[i] * w2[i];
		}
		f2 += b2;
		g2 = sigmoid(f2);
		float error = label-(g2*4.0+1.0);
		insamplerror = sqrt(0.999*insamplerror*insamplerror + 0.001 *(error*error));



		for(int i = 0 ; i < n_neuron; i++)
		{
			w2[i] += 4.0 * lr* error* g2 * (1.0 - g2) * g1[i];
		}
		b2 += 4.0 * lr * error * g2 * (1.0 - g2) ;
		for(int i = 0 ; i < n_neuron; i++)
		{
			for(int j = 0; j < n_pred; j++)
			{
				w1[i][j] += 4.0 * lr* error* g2 * (1.0 - g2) * w2[i]* g1[i] *(1-g1[i])*data[j];
			}
			b1[i] +=  4.0 * lr * error * g2 * (1.0-g2) * w2[i] * g1[i] * (1-g1[i]);
		}
	}
	double compute_rmse(vector<vector<float>>& dataset, vector<float> &label)
	{
		double error = 0.0;
		int n = dataset.size();
		for(int i = 0; i < dataset.size(); i++)
		{
			error += pow((predict(dataset[i]) - label[i]),2.0)/n; 
		}
		return sqrt(error);
	}

};

void load_data_set(ifstream& a,int n, vector<vector<float>>&x, vector<float>&y)
{
	vector<float> tmp(n,0);;
	float tmp2;
	while(!(a.eof()))
	{
		a >> tmp2;   //assume the first one is label
		for(int i = 0; i < n; i++)  a >> tmp[i];
		x.push_back(tmp);
		y.push_back(tmp2);	
	}
	cout<<"set size" << x.size() << endl;
}


int main()
{
	int n_neuron = 50;
	int n_data;
	int n_pred = 68;
//	vector<vector<float>> w1(n_predictors,vector<float>(n_neurons,0));
//	vector<float> w2(n_neuron,0);
	vector<vector<float>> trainset;
	vector<float> trainlabel;
	vector<vector<float>> testset;
	vector<float> testlabel;

	ifstream a("./trainset.txt");
	ifstream b("./testset.txt");
	load_data_set(a, n_pred, trainset, trainlabel);
	load_data_set(b, n_pred, testset,  testlabel);
	nn blender(n_neuron, n_pred);
	int epoch =100;
	for(int i=0; i < epoch; i++)
	{
		for(int j =0; j < trainset.size(); j++)
		{
			if(j%100000 == 0)
			{
				cout << j << endl;
				blender.output();
			}
			blender.sgd(trainset[j], trainlabel[j]);
		}
		cout << "rmse:" << blender.compute_rmse(testset,testlabel);
		blender.ratedecay();

	}

}
