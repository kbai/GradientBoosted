#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <thread>
#include <fstream>
#include <ctime>
#include <mutex>
#define NMOVIE 17770
#define NUSER 458293
#define NTIME 2243
#define NUMTH 1
using namespace std;
mutex mtx;

default_random_engine e(time(0));
float sigmoid(float x)
{
	return 1./(1.+exp(-x));
}
class nn
{
	int n_neuron;
	int n_pred;
	double insample;
	vector<vector<float>> w1;
	vector<float> w2;
	vector<float> b1;
	vector<vector<float>> bu;
	vector<vector<float>> bm;
	vector<vector<float>> bt;
	vector<vector<vector<float>>> bif;
	float b2 = 0;
	float lr = 0.0005;
	public:
	nn(int _neuron, int _pred):
		w1(_neuron,vector<float>(_pred,0)),
		w2(_neuron, 0),
		b1(_pred,0),
		bu(NUSER,vector<float>(_pred,0)),
		bm(NMOVIE,vector<float>(_pred,0)),
		bt(NTIME,vector<float>(_pred,0)),
		bif(NMOVIE, vector<vector<float>>(30, vector<float>(_pred,0)))
	{
		n_neuron = _neuron;
		n_pred = _pred;
		uniform_real_distribution<float> u(-0.5,0.5);
		for(auto&a : w1) for(auto &b : a) b =u(e);
		for(auto&a : w2) a = u(e);	
	}
	void output()
	{
		cout << "w1:" <<  w1[e()%n_neuron][e()%n_pred] << endl;
		cout << "w2:" << w2[e()%n_pred] << endl;
		cout << "b1:" << b1[e()%n_pred] << endl;
		cout << "b2:" << b2 << endl;
	}
	float ratedecay()
	{
		lr -= 1e-6;
	}
	void set0()
	{
		insample = 0.0;
	}
	double get_insample()
	{
		return insample;
	}
	

	float predict(vector<float> & data,int iu, int im, int it, int ife)
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
			g1[i] = sigmoid(f1[i]+b1[i] + bm[im][i] + bu[iu][i] + bt[it][i] + bif[im][ife][i]);
		}
		for(int i = 0; i < n_neuron; i++)
		{
			f2 += g1[i] * w2[i];
		}
		f2 += b2;
		g2 = sigmoid(f2);
		return 4.0 * g2 + 1.0;

	}
	float sgd(vector<float>& data,int iu, int im,int it,int ife, float label)
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
			g1[i] = sigmoid(f1[i]+b1[i]+ bm[im][i] + bu[iu][i] + bt[it][i] + bif[im][ife][i]);
		}
		for(int i = 0; i < n_neuron; i++)
		{
			f2 += g1[i] * w2[i];
		}
		f2 += b2;
		g2 = sigmoid(f2);
		float error = label-(g2*4.0+1.0);
		insample += error*error;
//		insamplerror = sqrt(0.999*insamplerror*insamplerror + 0.001 *(error*error));


		//mtx.lock();
		float con = 4.0 * lr* error* g2 * (1.0 - g2); 


		for(int i = 0 ; i < n_neuron; i++)
		{
			w2[i] += con * g1[i];
		}
		b2 += con;

		for(int i = 0 ; i < n_neuron; i++)
		{
			float con2 = con *w2[i]* g1[i] *(1-g1[i]);
			for(int j = 0; j < n_pred; j++)
			{
				w1[i][j] +=con2 * data[j];
			}
			b1[i] +=  con2;
			bm[im][i] += 0.1 * (con2 - bm[im][i] * lr * 0.01);
			bu[iu][i] += 0.5 * (con2 - bu[iu][i] * lr * 0.01);
			bif[im][ife][i] += 0.5 * (con2 - bif[im][ife][i] * lr * 0.01);
			bt[it][i] += 0.005 * (con2 - bt[it][i] * lr * 0.001);
	
		}
	//	mtx.unlock();
	}
	double compute_rmse(vector<vector<float>>& dataset,vector<int>&iu, vector<int>&im,vector<int>&it, vector<int>&ife, vector<float> &label)
	{
		double error = 0.0;
		int n = dataset.size();
		for(int i = 0; i < dataset.size(); i++)
		{
			error += pow((predict(dataset[i],iu[i],im[i], it[i],ife[i]) - label[i]),2.0)/n; 
		}
		return sqrt(error);
	}

};

void load_data_set(ifstream& a,int n, vector<vector<float>>&x,vector<int>&iu, vector<int>&im,vector<int>&it,vector<int>&ife, vector<float>&y)
{
	vector<float> tmp(n,0);;
	float tmp2;
	float tu,tm,tt,fe;
	int p = 0;
	while(!(a.eof()))
	{
		if(p%10000 ==0) cout << p << endl;
		p++;
		a >> tmp2;   //assume the first one is label
		for(int i = 0; i < n; i++)  a >> tmp[i];
		x.push_back(tmp);
		a >> tu >> tm >> tt>> fe;
		if((int) tu % 10000 == 0) cout << tu << endl;
		iu.push_back((int)tu-1);
		im.push_back((int)tm-1);
		it.push_back((int)tt-1);
		ife.push_back((int)(log(fe)/log(6.86)));
		y.push_back(tmp2);	
	}
	cout<<"set size" << x.size() << endl;
}

void load_quiz_set(ifstream& c,int n, vector<vector<float>>&x,vector<int>&iu, vector<int>&im, vector<int>&it,vector<int>& ife)
{
	vector<float> tmp(n,0);;
	int a = 0;
	float tu,tm, tt, fe;
	while(!(c.eof()))
	{
		for(int i = 0; i < n; i++)  c >> tmp[i];
		c>> tu >> tm >> tt >> fe;
		iu.push_back((int)tu-1);
		im.push_back((int)tm-1);
		it.push_back((int)tt-1);
		ife.push_back((int)(log(fe)/log(6.86)));
		a++;
		if(a%10000 == 0) cout << a<< endl;
		x.push_back(tmp);
	}
	cout<<"quiz set size" << x.size() << endl;
}

void shuffle(vector<int>& order)
{
	int ofset;
	int tmp;
	int ind;
	for(int i = 0; i < order.size(); i++)
	{
		ofset = order.size() - i;
		ind = i+e()%ofset;
		tmp = order[i];
		order[i] = order[ind];
		order[ind] = tmp;
	}
}


int main(int argc, char** argv)
{
	e.seed(atoi(argv[2]));
	int n_neuron = 60;
	int n_data;
	int n_pred = 101;
//	vector<vector<float>> w1(n_predictors,vector<float>(n_neurons,0));
//	vector<float> w2(n_neuron,0);
	vector<vector<float>> trainset;
	vector<float> trainlabel;
	vector<int> trainiu;
	vector<int> trainim;
	vector<int> trainit;
	vector<int> trainife;
	vector<vector<float>> testset;
	vector<float> testlabel;
	vector<int> testiu;
	vector<int> testim;
	vector<int> testit;
	vector<int> testife;
	vector<vector<float>> quizset;
	vector<float> quizlabel;
	vector<int> quiziu;
	vector<int> quizim;
	vector<int> quizit;
	vector<int> quizife;
	double olderror;
	double newerror;
	olderror = 1e8;
	newerror = 0.0;

	ifstream a("./trainset.txt");
	ifstream b("./testset.txt");
	ifstream c("./quizset.txt");
	load_data_set(a, n_pred, trainset, trainiu, trainim,trainit,trainife, trainlabel);
	load_data_set(b, n_pred, testset,  testiu, testim, testit,testife, testlabel);
	load_quiz_set(c, n_pred, quizset, quiziu, quizim , quizit, quizife); 
	vector<int> order(trainset.size(),0);
	int p = 0;
	for(int &u: order)
	{
		u = p;
		p++;
	}


	std::thread threads[NUMTH];
	nn blender(n_neuron, n_pred);
	int epoch =1000;
	for(int i=0; i < epoch; i++)
	{
//		shuffle(order);
		blender.set0();
		for(int j =0; j < trainset.size(); j+= NUMTH)
		{
			if(j%100000 == 0)
			{
				cout << j << endl;
				blender.output();
			}
			
//			blender.sgd(trainset[j], trainlabel[j]);
			for(int k = j ; k < j+NUMTH; k++)
			{
				if(k < trainset.size())
					blender.sgd(trainset[k], trainiu[k], trainim[k],trainit[k],trainife[k], trainlabel[k]);
	//		threads[k-j] = std::thread(&nn::sgd, &blender,std::ref(trainset[order[k]]),trainlabel[order[k]] );
			}
			//for(auto &a : threads) if(a.joinable()) a.join();
		}


		newerror =  blender.compute_rmse(testset,testiu, testim,testit,testife, testlabel);
		cout << "rmse:" << newerror<< endl;
		cout <<"in sample error:"<< sqrt(blender.get_insample()/trainset.size()) << endl;


		blender.ratedecay();
		if(i%20 == 1)
		{
			ofstream d(std::string(argv[1])+std::to_string(i)+".txt");
			for(int p = 0; p < quizset.size(); p++)
				d<< blender.predict(quizset[p], quiziu[p], quizim[p],quizit[p],quizife[p]) << endl;
			d.close();
		}
		olderror = newerror;



	}
	ofstream d(std::string(argv[1]) + std::to_string(100)+".txt");
	for(int p = 0; p < quizset.size(); p++)
				d<< blender.predict(quizset[p], quiziu[p], quizim[p], quizit[p], quizife[p]) << endl;

	cout <<"in sample error:"<< sqrt(blender.get_insample()/trainset.size()) << endl;

}
