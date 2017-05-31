#include "commonheader.h"
#include <iostream>
#include <math.h>
#include "utils.h"
using namespace std;

#ifdef MPIU
#endif
#include "RBM.h"
using namespace utils;

double compute_loss_train(
			vector<vector<int>>&train,
			vector<vector<int>>&train_index, 
			RBM& testrbm)
{
	vector<double> reconstructed_X(NMOVIE,0);
	int counter = 0;
	double error = 0;

	for(int i = 0 ; i < train.size(); i++)
	{
		testrbm.reconstruct(train[i],train_index[i],reconstructed_X);
		counter += train[i].size();
		for(int j = 0; j < train[i].size(); j++)
		error += pow((train[i][j] - reconstructed_X[train_index[i][j]]),2.0);
	}
	error /= counter;
	error = sqrt(error);
	return error;
}


double compute_loss(vector<int>&userid,
			vector<vector<int>>&train,
			vector<vector<int>>&train_index, 
			vector<vector<int>>&test,
			vector<vector<int>>& test_index,  
			RBM& testrbm)
{
	vector<double> reconstructed_X(NMOVIE,0);
	int counter = 0;
	double error = 0;

	for(int i = 0 ; i < test.size(); i++)
	{
		testrbm.reconstruct(train[userid[i]],train_index[userid[i]],test_index[i],reconstructed_X);
		counter += test[i].size();
		for(int j = 0; j < test[i].size(); j++)
		error += pow((test[i][j] - reconstructed_X[test_index[i][j]]),2.0);
	}
	error /= counter;
	error = sqrt(error);
	return error;
}

RBM::RBM(int size, int n_v, int n_h, double **w, double *hb, double *vb) {
  N = size;
  n_visible = n_v;
  n_hidden = n_h;

#ifdef MPIU
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#endif
  if(w == NULL) {
    W = new double*[n_hidden];
	dW = new double*[n_hidden];
	sW = new double*[n_hidden];

    for(int i=0; i<n_hidden; i++)
	{
		try{
		W[i] = (double*) malloc(sizeof(double)*n_visible);
		dW[i] = (double*) malloc(sizeof(double)*n_visible);
		sW[i] = (double*) malloc(sizeof(double)*n_visible);
		}
		catch(std::bad_alloc&)
		{
			cout << "bad allocation error!" << endl;
		}
	}
    double a = 1.0 / n_visible;

    for(int i=0; i<n_hidden; i++) {
      for(int j=0; j<n_visible; j++) {
        W[i][j] = uniform(-a, a);
		dW[i][j] = 0.0;
		sW[i][j] = 0.0;
      }
    }
  } else {
    W = w;

  }
#ifdef MPIU
    for(int u = 0; u < n_hidden; u++) MPI_Bcast(W[u],n_v,MPI_DOUBLE, 0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	
	if(hb == NULL) {
    hbias = new double[n_hidden];
	dh = new double[n_hidden];
	sh = new double[n_hidden];

    for(int i=0; i<n_hidden; i++){
	   	hbias[i] = 0;
		dh[i]=0;
		sh[i]=0;
	}
  } else {
    hbias = hb;
  }

  if(vb == NULL) {
    vbias = new double[n_visible];
	dv = new double[n_visible];
	sv = new double[n_visible];
    for(int i=0; i<n_visible; i++) {vbias[i] = 0;
		dv[i] = 0;
		sv[i] = 0;
	}
  } else {
    vbias = vb;
  }
}

RBM::~RBM() {
  for(int i=0; i<n_hidden; i++) delete[] W[i];
  delete[] W;
  delete[] hbias;
  delete[] vbias;
}


void RBM::contrastive_divergence(std::vector<int>& input ,std::vector<int>& imovie, double lr, int k) {
  double *ph_mean = new double[n_hidden];
  int *ph_sample = new int[n_hidden];
  double *nv_means = new double[n_visible];
  int *nv_samples = new int[n_visible];
  double *nh_means = new double[n_hidden];
  int *nh_samples = new int[n_hidden];

  /* CD-k */
  sample_h_given_v(&(input[0]), imovie, ph_mean, ph_sample);

  for(int step=0; step<k; step++) {
    if(step == 0) {
      gibbs_hvh(ph_sample,imovie, nv_means, nv_samples, nh_means, nh_samples);
    } else {
      gibbs_hvh(nh_samples,imovie,  nv_means, nv_samples, nh_means, nh_samples);
    }
  }

  for(int i=0; i<n_hidden; i++) {
    for(int jr=0; jr<imovie.size(); jr++) {
		int j = imovie[jr];
//      dW[i][j] += lr * (ph_mean[i] * input[jr] - nh_means[i] * nv_means[jr]) / N;
      W[i][j] += lr * (ph_mean[i] * input[jr] - nh_means[i] * nv_means[jr]) / N;

//cout <<"\t"<<imovie.size() <<"\t" <<ph_mean[i] <<"\t" <<input[jr] <<" "<< nh_means[i]<<"\t"<< nv_means[jr] << endl;
    }
    hbias[i] += lr * (ph_mean[i] - nh_means[i]) / N;
	
  }

  for(int ir=0; ir<imovie.size(); ir++) {
	  int i = imovie[ir];
    vbias[i] += lr * (input[ir] - nv_means[ir]) / N;
  }

  
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

void RBM::sample_h_given_v(int *v0_sample, vector<int>& imovie,double *mean, int *sample) {
  for(int i=0; i<n_hidden; i++) {
    mean[i] = propup(v0_sample,imovie, W[i], hbias[i]);
    sample[i] = binomial(1, mean[i]);
  }
}

void RBM::sample_v_given_h(int *h0_sample, vector<int>& imovie,double *mean, int *sample) {
	int i ;
  for(int ir=0; ir<imovie.size(); ir++) {
	i = imovie[ir];
    mean[ir] = propdown(h0_sample, i, vbias[i]);
    sample[ir] = binomial(1, mean[ir]);
  }
}

double RBM::propup(int *v,vector<int> &imovie, double *w, double b) {
  double pre_sigmoid_activation = 0.0;
  for(int j=0; j<imovie.size(); j++) {
    pre_sigmoid_activation += w[imovie[j]] * v[j];
  }
  pre_sigmoid_activation += b;
  return sigmoid(pre_sigmoid_activation);
}

double RBM::propdown(int *h, int i, double b) {
  double pre_sigmoid_activation = 0.0;
  for(int j=0; j<n_hidden; j++) {
    pre_sigmoid_activation += W[j][i] * h[j];
  }
  pre_sigmoid_activation += b;
  return sigmoid(pre_sigmoid_activation);
}

void RBM::gibbs_hvh(int *h0_sample,vector<int>& imovie, double *nv_means, int *nv_samples, \
                    double *nh_means, int *nh_samples) {
  sample_v_given_h(h0_sample, imovie, nv_means, nv_samples);
  sample_h_given_v(nv_samples, imovie,  nh_means, nh_samples);
}

void RBM::reconstruct(vector<int>& vv,vector<int>& imovie, vector<double>& reconstructed_v) {
  int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation;

  for(int i=0; i<n_hidden; i++) {
    h[i] = propup(v,imovie,  W[i], hbias[i]);
  }

  for(int i=0; i<n_visible; i++) {
    pre_sigmoid_activation = 0.0;
    for(int j=0; j<n_hidden; j++) {
      pre_sigmoid_activation += W[j][i] * h[j];
    }
    pre_sigmoid_activation += vbias[i];

    reconstructed_v[i] = sigmoid(pre_sigmoid_activation);
  }

  delete[] h;
}




void RBM::reconstruct(vector<int>& vv,vector<int>& imovie, vector<int>& imovie2, vector<double>& reconstructed_v) {
  int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation;

  for(int i=0; i<n_hidden; i++) {
    h[i] = propup(v,imovie,  W[i], hbias[i]);
  }

  for(int i=0; i<imovie2.size(); i++) {
	  int im = imovie2[i];
    pre_sigmoid_activation = 0.0;
    for(int j=0; j<n_hidden; j++) {
      pre_sigmoid_activation += W[j][im] * h[j];
    }
    pre_sigmoid_activation += vbias[im];
    reconstructed_v[i] = sigmoid(pre_sigmoid_activation);
  }

  delete[] h;
}


void load_all_data(
				vector<vector<int>>& train,
				vector<vector<int>>& train_index,  
				vector<vector<int>>& test,	
				vector<vector<int>>& test_index,
				vector<int>& uid_test)
{
	ifstream traindat("../../data/1.dta");
	ifstream testdat("../../data/2.dta");
	vector<int> a;
	vector<int> b;
	vector<int> tmp(5,0);
	int iu = 1;

    while(traindat >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] && iu < 100)
    {
		if(iu < tmp[0])
		{
			train.push_back(a);
			train_index.push_back(b);
			a.clear();
			b.clear();
			iu = tmp[0];
		}
		a.push_back(tmp[4]>4.5);
		b.push_back(tmp[1]-1);
    }
	train.push_back(a);
	train_index.push_back(b);
	a.clear();
	b.clear();
	iu = -1;

	cout << "Training set loaded!" << train.size() <<"data" << endl;

    while(testdat >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4])
    {
		if(iu == -1) iu = tmp[0];
		if(iu < tmp[0])
		{
			uid_test.push_back(iu-1);
			test.push_back(a);
			test_index.push_back(b);
			a.clear();
			b.clear();
			iu = tmp[0];
		}
		a.push_back(tmp[4]-1);
		b.push_back(tmp[1]-1);
    }
	test.push_back(a);
	test_index.push_back(b);
	a.clear();
	b.clear();
	cout << "Test set loaded!" << test.size() <<"data" << endl;
}

void test_rbm() {
#ifdef MPIU
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#endif
	int train_N = 100;
	int test_N = 2;
	int n_visible = NMOVIE;
	int n_hidden = 100;
	double learning_rate = 1.0;
	int training_epochs = 10000;
	int k = 7;
	vector<int> testuid;
	vector<vector<int>> trainset,testset,train_index, test_index;
	load_all_data(trainset, train_index, testset, test_index, testuid);
	RBM rbm(train_N,  n_visible, n_hidden, NULL, NULL, NULL);
#ifdef MPIU
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#endif
	int iu = 0;
//	int len = trainset.size();
	int len = 100;
	for(int y = 0 ; y < 6000; y++)
	{
//		cout << y << endl;
		rbm.contrastive_divergence(trainset[y%len], train_index[y%len], learning_rate, k);

#ifdef MPIU
		MPI_Barrier(MPI_COMM_WORLD); 
		rbm.mpigath();
#endif
		if(y%100 == 0) cout <<compute_loss_train(trainset,train_index,rbm) << endl;

	}



	int counter = 0;
	vector<double> reconstructed_X(NMOVIE,0.0);
	rbm.reconstruct(trainset[0], train_index[0] , reconstructed_X);
	for(int i = 0; i < trainset[0].size(); i++)
	{
	   cout << trainset[0][i] << "\t" << reconstructed_X[train_index[0][i]]<< endl;
	}



	double maxv = -999.0,minv =999.0;
}


int main() {
#ifdef MPIU
	MPI_Init(NULL,NULL);
	char pron[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(pron,&name_len);
	cout << pron << endl;
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	cout << myrank << endl;
#endif
	  test_rbm();
#ifdef MPIU
	MPI_Barrier(MPI_COMM_WORLD);
#endif
  return 0;
}


