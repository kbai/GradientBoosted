#include "commonheader.h"
#include <iostream>
#include <math.h>
#include <assert.h>
#include "utils.h"
#include "ctime"
#include <chrono>
//#include <mutex>
#define NUMTH 50
#include <thread>
#include <mutex>
using namespace std;
#ifdef MPIU
#endif
#include "RBM5.h"
#define NLEN
#include <cstring>
#include <time.h>
#include <sys/time.h>
using namespace utils;

std::mutex mtx;



double get_wall_time(){
	    struct timeval time;
		gettimeofday(&time,NULL);
		return (double)time.tv_sec + (double)time.tv_usec * .000001;
}



double compute_loss_train(
			const vector<vector<int>>&train,
			const vector<vector<int>>&train_index, 
			const RBM& testrbm)
{
	vector<double> reconstructed_X(5*NMOVIE,0);
	int counter = 0;
	double error = 0;


	for(int i = 0 ; i < train_index.size(); i++)
	{
		testrbm.reconstruct(train[i],train_index[i],reconstructed_X);
		counter += train_index[i].size();
		for(int j = 0; j < train_index[i].size(); j++)
		{
			double score = 0.0;
			double tscore = 0.0;
			for(int ns = 0; ns < 5; ns++)
			{
				score += ( ns + 1 )* (reconstructed_X[5*train_index[i][j]+ns]);
				tscore += train[i][5*j+ns] * (ns+1);
			}
			error += pow(score-tscore,2.0);
		}

	}
	error /= counter;
	error = sqrt(error);
	return error;
}

double compute_loss(	const vector<int>&trainid,
			const vector<int>&userid,
			const vector<vector<int>>&train,
			const vector<vector<int>>&train_index, 
			const vector<vector<int>>&test,
			const vector<vector<int>>& test_index,  
			const RBM& testrbm,
			const string& filename)
{
	vector<double> reconstructed_X(5*NMOVIE,0);
	int counter = 0;
	double error = 0;
	vector<int> idmap;
	int j = 0;
	ofstream output(filename);
	for(auto i: userid)
	{
		while(trainid[j] < i) j++;
		if(trainid[j] > i) 
		{
			cout << trainid[j-1]<<endl;
			cout << trainid[j]<< endl;
			cout << i << endl; 
			exit(1);
		}//some id in test set cannot be found in the training set. return an error;
		idmap.push_back(j);
	}

	for(int i = 0 ; i < test_index.size(); i++)
	{
		testrbm.reconstruct(train[idmap[i]],train_index[idmap[i]],test_index[i],reconstructed_X);
		counter += test_index[i].size()/10;
		for(int j = 0; j < test_index[i].size()/10; j++)
		{
			double score = 0.0;
			for(int ns = 0; ns < 5; ns++)
			{
				score += ( ns + 1 )* (reconstructed_X[5*test_index[i][10*j]+ns] - test[i][5*j+ns]);
				output << reconstructed_X[5*test_index[i][10*j]+ns] <<"\t" ;
			}
			output << endl;
			error += pow(score,2.0);
		}
	}
	error /= counter;
	error = sqrt(error);
	output.close();
	cout << "error:" << filename<<"\t" <<error << endl;
	return error;
}

RBM::RBM(int size, int n_v, int n_h, double **w, double *hb, double *vb) {
  N = size;
  n_visible = n_v;
  n_hidden = n_h;

    W = new double*[n_hidden];

    for(int i=0; i<n_hidden; i++)
	{
		W[i] = (double*) malloc(5*sizeof(double)*n_visible);
	}
    double a = 1.0 / n_visible;
	
	bif = new double[5*NMOVIE*30];
	but = new double[5*NLARGE];
	bu = new double[5*NUSER];
    hbias = new double[n_hidden];
    vbias = new double[5*n_visible];

	for(int i=0; i<NUSER; i++) bu[i] = 0.0;
	for(int i=0; i< NLARGE; i++) but[i] = 0.0;
	for(int i=0; i< NMOVIE*30;i++) bif[i] = 0.0;
    for(int i=0; i<n_hidden; i++) for(int j=0; j<5*n_visible; j++) W[i][j] = uniform(-a, a);
    for(int i=0; i<n_hidden; i++) hbias[i] = 0;
    for(int i=0; i< 5 * n_visible; i++) vbias[i] = 0;




    vW = new double*[n_hidden];


    for(int i=0; i<n_hidden; i++)
	{
		vW[i] = (double*) malloc(5*sizeof(double)*n_visible);
	}
	vbif = new double[5*NMOVIE*30];
	vbut = new double[5*NLARGE];
	vbu = new double[5*NUSER];
    vhbias = new double[n_hidden];
    vvbias = new double[5*n_visible];

	for(int i=0; i<NUSER; i++) vbu[i] = 0.0;
	for(int i=0; i< NLARGE; i++) vbut[i] = 0.0;
	for(int i=0; i< NMOVIE*30;i++) vbif[i] = 0.0;
    for(int i=0; i<n_hidden; i++) for(int j=0; j<5*n_visible; j++) vW[i][j] = 0.0;
    for(int i=0; i<n_hidden; i++) vhbias[i] = 0;
    for(int i=0; i< 5 * n_visible; i++) vvbias[i] = 0;



}
RBM::RBM(const RBM& a) {
  N = a.N;
  n_visible = a.n_visible;
  n_hidden = a.n_hidden;

    W = new double*[n_hidden];

    for(int i=0; i<n_hidden; i++)
	{
		W[i] = (double*) malloc(5*sizeof(double)*n_visible);
		memcpy(W[i],a.W[i],sizeof(double)*5*n_visible);
	}
	
	bif = new double[5*NMOVIE*30];
	but = new double[5*NLARGE];
	bu = new double[5*NUSER];
    hbias = new double[n_hidden];
    vbias = new double[5*n_visible];
	memcpy(bif		,a.bif		,5*NMOVIE*30*sizeof(double));
	memcpy(but		,a.but		,5*NLARGE*sizeof(double));
	memcpy(bu		,a.bu		,5*NUSER*sizeof(double));
	memcpy(hbias	,a.hbias	,n_hidden*sizeof(double));
	memcpy(vbias	,a.vbias	,5*n_visible*sizeof(double));

    vW = new double*[n_hidden];

    for(int i=0; i<n_hidden; i++)
	{
		vW[i] = (double*) malloc(5*sizeof(double)*n_visible);
		memcpy(vW[i],a.vW[i],sizeof(double)*5*n_visible);
	}
	
	vbif = new double[5*NMOVIE*30];
	vbut = new double[5*NLARGE];
	vbu = new double[5*NUSER];
    vhbias = new double[n_hidden];
    vvbias = new double[5*n_visible];
	memcpy(vbif		,a.vbif		,5*NMOVIE*30*sizeof(double));
	memcpy(vbut		,a.vbut		,5*NLARGE*sizeof(double));
	memcpy(vbu		,a.vbu		,5*NUSER*sizeof(double));
	memcpy(vhbias	,a.vhbias	,n_hidden*sizeof(double));
	memcpy(vvbias	,a.vvbias	,5*n_visible*sizeof(double));

}

RBM::~RBM() {
  for(int i=0; i<n_hidden; i++) delete[] W[i];
  delete[] W;
  delete[] bif;
  delete[] but;
  delete[] bu;
  delete[] hbias;
  delete[] vbias;
  delete[] vW;
  delete[] vbif;
  delete[] vbut;
  delete[] vbu;
  delete[] vhbias;
  delete[] vvbias;

}

void RBM::operator=(const RBM& a) {
  N = a.N;
  n_visible = a.n_visible;
  n_hidden = a.n_hidden;
    for(int i=0; i<n_hidden; i++)
	{
		memcpy(W[i],a.W[i],sizeof(double)*5*n_visible);
	}
	
	memcpy(bif		,a.bif		,5*NMOVIE*30*sizeof(double));
	memcpy(but		,a.but		,5*NLARGE*sizeof(double));
	memcpy(bu		,a.bu		,5*NUSER*sizeof(double));
	memcpy(hbias	,a.hbias	,n_hidden*sizeof(double));
	memcpy(vbias	,a.vbias	,5*n_visible*sizeof(double));
    for(int i=0; i<n_hidden; i++)
	{
		memcpy(vW[i],a.vW[i],sizeof(double)*5*n_visible);
	}
	
	memcpy(vbif		,a.vbif		,5*NMOVIE*30*sizeof(double));
	memcpy(vbut		,a.vbut		,5*NLARGE*sizeof(double));
	memcpy(vbu		,a.vbu		,5*NUSER*sizeof(double));
	memcpy(vhbias	,a.vhbias	,n_hidden*sizeof(double));
	memcpy(vvbias	,a.vvbias	,5*n_visible*sizeof(double));



}


void RBM::contrastive_divergence(const std::vector<int>& input ,const std::vector<int>& imovie, const double lr,const double momentum, const int k) {	
  double *ph_mean = new double[n_hidden];
  int *ph_sample = new int[n_hidden];
  double *nv_means = new double[5*n_visible];
  int *nv_samples = new int[5*n_visible];
  double *nh_means = new double[n_hidden];
  int *nh_samples = new int[n_hidden];

  /* CD-k */
  sample_h_given_v(&(input[0]), imovie, ph_mean, ph_sample);
	  //cout << "cd0.5" << endl;

  for(int step=0; step<k; step++) {
    if(step == 0) {
      gibbs_hvh(ph_sample,imovie, nv_means, nv_samples, nh_means, nh_samples);
    } else {
      gibbs_hvh(nh_samples,imovie,  nv_means, nv_samples, nh_means, nh_samples);
    }
  }

  //cout << "cd1" << endl;
	mtx.lock();

  for(int i=0; i<n_hidden; i++) {
    for(int jr=0; 10 * jr < imovie.size(); jr++) {
		int j = imovie[ 10 * jr];
		for(int ns = 0; ns < 5; ns ++)
		{
//      dW[i][j] += lr * (ph_mean[i] * input[jr] - nh_means[i] * nv_means[jr]) / N;
#ifdef MM
			vW[i][5*j+ns] = lr * (ph_mean[i] * input[5*jr + ns] - nh_means[i] * nv_means[5*jr+ns] - 0.001*W[i][5*j+ns]) + momentum* vW[i][5*j+ns];
			W[i][5*j+ns] += vW[i][5*j+ns];
#else
			W[i][5*j + ns] += lr * (ph_mean[i] * input[5*jr + ns] - nh_means[i] * nv_means[5*jr+ns] - 0.001*W[i][5*j+ns]);
#endif
//cout <<"\t"<<imovie.size() <<"\t" <<ph_mean[i] <<"\t" <<input[5*jr + ns] <<" "<< nh_means[i]<<"\t"<< nv_means[5*jr+ns] << endl;

		}
    }
#ifdef MM
    vhbias[i] = lr * (ph_mean[i] - nh_means[i] - 0.001 * hbias[i]) + momentum * vhbias[i];
    hbias[i] += vhbias[i]; 
#else
	hbias[i] +=  lr * (ph_mean[i] - nh_means[i] - 0.001 * hbias[i]); 
#endif	
  }
  //cout << "cd2" << endl;

  for(int ir=0; 10 * ir<imovie.size(); ir++) {
	int im  = imovie[10 * ir + 0];
	int ife = imovie[10 * ir + 3];
	int iut = imovie[10 * ir + 4];
	int it  = imovie[10 * ir + 2];
	int iu  = imovie[10 * ir + 1];
	
	  for(int ns = 0; ns < 5; ns++)
	  {
		double bias =  input[5*ir + ns] - nv_means[5 * ir + ns];
#ifdef MM
		vvbias[5*im + ns]       =  lr * (bias - 0.001*vbias[5*im + ns]) + momentum*vvbias[5*im+ns];
		vbias[5*im + ns]		+= vvbias[5*im + ns];
		vbu[5*iu + ns]			=  0.25*lr*(bias - 0.001* bu[5*iu+ns]) + momentum * vbu[5*iu + ns];
		bu[5*iu + ns]           += vbu[5*iu+ns];
		vbut[5*iut+ns]			=  0.80 *lr*(bias - 0.001* but[5*iut+ns]) + momentum * vbut[5*iut+ns];
		but[5*iut+ns]			+=	vbut[5*iut+ns];
		vbif[5*((30*im)+ife)+ns]=  0.02*lr*(bias - 0.001* bif[5*((30*im)+ife)+ns]) + momentum * vbif[5*((30*im)+ife)+ns]; 
		bif[5*((30*im)+ife)+ns] += vbif[5*((30*im)+ife)+ns];
#else
		vbias[5*im + ns]        +=  lr * (bias - 0.001 *vbias[5*im + ns]);
		bu[5*iu + ns]			+=  0.25*lr*(bias - 0.001* bu[5*iu+ns]);
		but[5*iut+ns]			+=  0.80*lr*(bias - 0.001* but[5*iut+ns]);
		bif[5*((30*im)+ife)+ns] +=  0.02*lr*(bias - 0.001* bif[5*((30*im)+ife)+ns]);
#endif
	  }
  }
  mtx.unlock();
 //cout << "cd3"<< endl; 
  delete[] ph_mean;
  delete[] ph_sample;
  delete[] nv_means;
  delete[] nv_samples;
  delete[] nh_means;
  delete[] nh_samples;
}
void RBM:: mpigath()
{
#ifdef MPIU

	for(int  i = 0; i < n_hidden; i++)
	{
		MPI_Allreduce(dW[i],sW[i],n_visible,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	}
  MPI_Allreduce(dh,sh,n_hidden,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(dv,sv,n_visible,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  for(int i=0; i<n_hidden; i++) {
    for(int j=0; j< n_visible; j++) {
      W[i][j] += sW[i][j] ;
	  W[i][j] *= 0.99;
	  dW[i][j] = 0.0;
	  sW[i][j] = 0.0;
    }
	hbias[i] += sh[i];
	hbias[i] *= 0.99;
	sh[i] = 0.0;dh[i] = 0.0;
  }

  for(int i=0; i < n_visible; i++) {
    vbias[i] += sv[i];
	vbias[i] *= 0.99;
	sv[i] = 0.0;dv[i] =0.0;
  }
#endif
}

void RBM::sample_h_given_v(const int *v0_sample, const vector<int>& imovie,double *mean, int *sample) const
{
  for(int i=0; i<n_hidden; i++) {
    mean[i] = propup(v0_sample,imovie, W[i], hbias[i]);
    sample[i] = binomial(1, mean[i]);
  }
}

void RBM::sample_v_given_h(const int *h0_sample,const vector<int>& imovie,double *mean, int *sample) const
{
  for(int ir=0; 10 * ir<imovie.size(); ir++) {

	int im  = imovie[10 * ir + 0];
	int ife = imovie[10 * ir + 3];
	int iut = imovie[10 * ir + 4];
	int it  = imovie[10 * ir + 2];
	int iu  = imovie[10 * ir + 1];
	
	for(int ns = 0; ns < 5; ns++)
	{
		propdown(h0_sample, im, vbias + 5 * im, bif + 5*( 30 * im + ife), but + 5* iut, bu+5*iu , mean + 5*ir);
		binomial5( mean+5*ir, sample+5*ir);
	}
  }
}

double RBM::propup(const int *v,const vector<int> &imovie, double *w, double b) const 
{

  double pre_sigmoid_activation = 0.0;

  for(int ns = 0; ns < 5; ns++)
  {
	for(int j=0; 10 * j < imovie.size(); j++) {
		pre_sigmoid_activation += w[5*imovie[ 10 * j ] + ns] * v[5*j + ns];
	}
  }
  pre_sigmoid_activation += b;

  return sigmoid(pre_sigmoid_activation);
}

void RBM::propdown(const int *h, const int i, const double* b,const double *b1, const double *b2,const  double*b3,  double *means) const
{
  double pre_sigmoid_activation[5] = {0.0};
//  cout <<"setp" << pre_sigmoid_activation[4] << endl;
  for(int ns = 0; ns < 5; ns++)
  {
	for(int j=0; j<n_hidden; j++) {
		pre_sigmoid_activation[ns] += W[j][5*i+ns] * h[j];
	}
	pre_sigmoid_activation[ns] += b[ns]+ b1[ns] + b2[ns] + b3[ns] ;
  }
	 //cout << pre_sigmoid_activation[0] <<" " << pre_sigmoid_activation[4] << endl;
	 sigmoid5(pre_sigmoid_activation,means);
	 return;
}

void RBM::gibbs_hvh(const int *h0_sample,const vector<int>& imovie, double *nv_means, int *nv_samples, \
                    double *nh_means, int *nh_samples) const
{
  sample_v_given_h(h0_sample, imovie, nv_means, nv_samples);
  sample_h_given_v(nv_samples, imovie,  nh_means, nh_samples);
}

void RBM::reconstruct(const vector<int>& vv,const vector<int>& imovie, vector<double>& reconstructed_v) const
{
  const int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation[5] = {0.0};

  for(int i=0; i<n_hidden; i++) {
    h[i] = propup(v,imovie,  W[i], hbias[i]);
  }

  for(int i=0; i<n_visible; i++) {
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] = 0.0;
	}

    for(int j=0; j<n_hidden; j++) {
		for(int ns = 0; ns < 5; ns++)
		{
			pre_sigmoid_activation[ns] += W[j][5*i+ns] * h[j];
		//	cout << W[j][5*i+ns] <<" \t"<<h[j] << endl;
		}
    }
	for(int ns=0;ns<5;ns++)  pre_sigmoid_activation[ns] += vbias[5*i + ns];
	//cout << pre_sigmoid_activation[0] << endl;
	sigmoid5(pre_sigmoid_activation,&(reconstructed_v[5*i]));

  }
    
  delete[] h;
}

void RBM::reconstruct(const vector<int>& vv,const vector<int>& imovie, const vector<int>& imovie2, vector<double>& reconstructed_v) const
{
  const int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation[5];

  for(int i=0; i<n_hidden; i++) 
  {
    h[i] = propup(v,imovie,  W[i], hbias[i]);
  }

  for(int i=0; 10 * i<imovie2.size(); i++) 
  {
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] = 0.0;
	}
	int im  = imovie2[10 * i + 0];
	int ife = imovie2[10 * i + 3];
	int iut = imovie2[10 * i + 4];
	int it  = imovie2[10 * i + 2];
	int iu  = imovie2[10 * i + 1];
	
    for(int j=0; j<n_hidden; j++) 
	{
		for(int ns = 0; ns < 5; ns++)
		{
			pre_sigmoid_activation[ns] += W[j][5*im + ns] * h[j];
		}
    }
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] += vbias[5*im + ns] + bu[ 5*iu + ns ] + but[ 5*iut + ns ] + bif[5 * (im * 30 + ife) + ns];
	}

    sigmoid5(pre_sigmoid_activation, &(reconstructed_v[5*im]));
  }

  delete[] h;
}


void load_all_data(
				vector<vector<int>>& train,
				vector<vector<int>>& train_index,  
				vector<int>& uid_train,
				vector<vector<int>>& test,	
				vector<vector<int>>& test_index,
				vector<int>& uid_test,
				vector<vector<int>>& quiz,
				vector<vector<int>>& quiz_index,
				vector<int>& uid_quiz)
{
	ifstream traindat("../../data/1map.dta");
	ifstream testdat("../../data/4map.dta");
	ifstream testdata2("../../data/5map.dta");
	vector<int> a;
	vector<int> b;
	vector<int> tmp(5,0);
	int iu = -1;

#define LOADDATA(A,B,C,D)\
    while(A >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] >> tmp[5])\
    {\
		if(iu == -1) iu = tmp[0];\
		if(iu < tmp[0])\
		{\
			B.push_back(a);\
			C.push_back(b);\
			D.push_back(iu-1);\
			a.clear();\
			b.clear();\
			iu = tmp[0];\
		}\
		for(int f = 1; f< 6; f++) a.push_back(tmp[5] == f); \
		b.push_back(tmp[1]-1);\
		b.push_back(tmp[0]-1);\
		b.push_back(tmp[2]-1);\
		b.push_back((int)(log(tmp[3])/log(2.718)));\
		b.push_back(tmp[4]);\
		for(int tmp=0; tmp < 5; tmp++) b.push_back(-1e8);\
    }\
	B.push_back(a);\
	C.push_back(b);\
	D.push_back(iu-1);\
	a.clear();\
	b.clear();\
	iu = -1;

	LOADDATA(traindat, train, train_index, uid_train)
	cout << "Training set loaded!" << train.size() <<"data" << endl;
	LOADDATA(testdat, test, test_index,uid_test)
	cout << "Test set loaded!" << test.size() <<"data" << endl;
	LOADDATA(testdata2, quiz, quiz_index, uid_quiz)
	cout << "Quiz set loaded!" << quiz.size() <<"data" << endl;
}

void test_rbm(string& f1, string& f2) {
	int test_N = 2;
	int n_visible = NMOVIE;
	int n_hidden = 200;
	int training_epochs = 10000;
	int k = 7;
	vector<int> testuid, trainuid, quizuid;
	vector<vector<int>> trainset,testset,train_index, test_index;
	vector<vector<int>> quizset, quiz_index;
	load_all_data(trainset, train_index,  trainuid, testset, test_index, testuid, quizset, quiz_index, quizuid);
	int train_N = trainset.size();
	double learning_rate = 0.01;

	RBM rbm(train_N,  n_visible, n_hidden, NULL, NULL, NULL);
	RBM rbm2(rbm);
	cout << "1111" << endl;
	int iu = 0;
//	int len = trainset.size();
	int len = train_N;

//    	 cout <<"error" << compute_loss(trainuid,testuid,trainset,train_index,testset,test_index,rbm) << endl;

	 std::thread threads[NUMTH];
	 std::thread writedata[2];

	auto begin = chrono::steady_clock::now();
	
	for(int iepo = 0 ; iepo < 100; iepo++)
	{

	//	if(iepo%2 == 0 && k < 9) k++;
		
		for(int y = 0 ; y < len; y+=NUMTH)
		{

		if(y%1000 == 0)
		{
			auto finish = chrono::steady_clock::now();
			double elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - begin).count();
			cout << y <<":" << elapsed <<"seconds"<< endl;
		}

			for(int it = 0; it < NUMTH; it++)
			{
				int iy = (y+it)%len;
				threads[it] = std::thread(&RBM::contrastive_divergence,&rbm,std::ref(trainset[iy]),std::ref(train_index[iy]), learning_rate,0.9 ,k);
			}
			for(auto& th:threads) th.join();
		//	cout << "elapsed time:"<< 1.0*(end-begin)/CLOCKS_PER_SEC/NUMTH<< endl;

		}
		learning_rate *= 0.91;
		//cout <<"4545" <<endl;
		 if(writedata[0].joinable()) writedata[0].join();
		 if(writedata[1].joinable()) writedata[1].join();
		 rbm2 = rbm;
		 writedata[0] = std::thread(compute_loss,std::ref(trainuid),std::ref(testuid),std::ref(trainset),std::ref(train_index),std::ref(testset),std::ref(test_index),std::ref(rbm2),std::ref(f1));
		 writedata[1] = std::thread(compute_loss,std::ref(trainuid),std::ref(quizuid),std::ref(trainset),std::ref(train_index),std::ref(quizset),std::ref(quiz_index),std::ref(rbm2),std::ref(f2));

	}
		 if(writedata[0].joinable()) writedata[0].join();
		 if(writedata[1].joinable()) writedata[1].join();
	
	double maxv = -999.0,minv =999.0;
}


int main(int argc, char** argv) {
	string f1 = std::string(argv[1]);
	string f2 = std::string(argv[2]);
	cout << f1 << f2 <<endl;

	  test_rbm(f1,f2);
  return 0;
}


