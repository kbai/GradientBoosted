#include "commonheader.h"
#include <iostream>
#include <math.h>
#include <assert.h>
#include "utils.h"
//#include <mutex>
#define NUMTH 100
#include <thread>
#include <mutex>
using namespace std;

#ifdef MPIU
#endif
#include "RBM5.h"
#define NLEN
using namespace utils;

std::mutex mtx;

double compute_loss_train(
			vector<vector<int>>&train,
			vector<vector<int>>&train_index, 
			RBM& testrbm)
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

double compute_loss(	vector<int>&trainid,
			vector<int>&userid,
			vector<vector<int>>&train,
			vector<vector<int>>&train_index, 
			vector<vector<int>>&test,
			vector<vector<int>>& test_index,  
			RBM& testrbm,
			string& filename)
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
		counter += test_index[i].size();
		for(int j = 0; j < test_index[i].size(); j++)
		{
			double score = 0.0;
			for(int ns = 0; ns < 5; ns++)
			{
				score += ( ns + 1 )* (reconstructed_X[5*test_index[i][j]+ns] - test[i][5*j+ns]);
				output << reconstructed_X[5*test_index[i][j]+ns] <<"\t" ;
			}
			output << endl;
			error += pow(score,2.0);
		}
	}
	error /= counter;
	error = sqrt(error);
	output.close();
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
	D = new double*[n_hidden];
	dW = new double*[n_hidden];
	sW = new double*[n_hidden];

    for(int i=0; i<n_hidden; i++)
	{
		try{
		W[i] = (double*) malloc(5*sizeof(double)*n_visible);
		D[i] = (double*) malloc(sizeof(double)*n_visible);
#ifdef MPIU
		dW[i] = (double*) malloc(5*sizeof(double)*n_visible);
		sW[i] = (double*) malloc(5*sizeof(double)*n_visible);
#endif
		}
		catch(std::bad_alloc&)
		{
			cout << "bad allocation error!" << endl;
		}
	}
    double a = 1.0 / n_visible;

    for(int i=0; i<n_hidden; i++) {
      for(int j=0; j<5*n_visible; j++) {
        W[i][j] = uniform(-a, a);
#ifdef MPIU
		dW[i][j] = 0.0;
		sW[i][j] = 0.0;
#endif
      }
	  for(int j = 0; j < n_visible; j++) D[i][j] = 0.0;

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
#ifdef MPIU
	dh = new double[n_hidden];
	sh = new double[n_hidden];
#endif
    for(int i=0; i<n_hidden; i++){
	   	hbias[i] = 0;
#ifdef MPIU
		dh[i]=0;
		sh[i]=0;
#endif
	}
  } else {
    hbias = hb;
  }

  if(vb == NULL) {
    vbias = new double[5*n_visible];
#ifdef MPIU
	dv = new double[5*n_visible];
	sv = new double[5*n_visible];
#endif
    for(int i=0; i< 5 * n_visible; i++) {vbias[i] = 0;
#ifdef MPIU
		dv[i] = 0;
		sv[i] = 0;
#endif
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


void RBM::contrastive_divergence(std::vector<int>& input ,std::vector<int>& imovie, double lr, double lbd,int k) {
  double *ph_mean = new double[n_hidden];
  int *ph_sample = new int[n_hidden];
  double *nv_means = new double[5*n_visible];
  int *nv_samples = new int[5*n_visible];
  double *nh_means = new double[n_hidden];
  int *nh_samples = new int[n_hidden];

  /* CD-k */
  sample_h_given_v(&(input[0]), imovie, ph_mean, ph_sample);
//  cout << "cd0.5" << endl;

  for(int step=0; step<k; step++) {
    if(step == 0) {
      gibbs_hvh(ph_sample,imovie, nv_means, nv_samples, nh_means, nh_samples);
    } else {
      gibbs_hvh(nh_samples,imovie,  nv_means, nv_samples, nh_means, nh_samples);
    }
  }

//  cout << "cd1" << endl;
	mtx.lock();

  for(int i=0; i<n_hidden; i++) {
    for(int jr=0; jr<imovie.size(); jr++) {
		int j = imovie[jr];
		for(int ns = 0; ns < 5; ns ++)
		{
			W[i][5*j + ns] += lr * (ph_mean[i] * input[5*jr + ns] - nh_means[i] * nv_means[5*jr+ns] - lbd * W[i][5*j + ns]) ;
		}
		D[i][j] += lr * (ph_mean[i] - nh_means[i] - lbd * D[i][j]);

    }
    hbias[i] += lr * (ph_mean[i] - nh_means[i]  - lbd * hbias[i]);
	
  }
//  cout << "cd2" << endl;

  for(int ir=0; ir<imovie.size(); ir++) {
	  int i = imovie[ir];
	  for(int ns = 0; ns < 5; ns++)
	  {
		vbias[5*i + ns] += lr * (input[5*ir + ns] - nv_means[5 * ir + ns]  - lbd * vbias[5*i + ns]);
	  }
  }
  mtx.unlock();

// cout << "cd3"<< endl; 
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
    mean[i] = propup(v0_sample,imovie, W[i], D[i], hbias[i]);
    sample[i] = binomial(1, mean[i]);
  }
}

void RBM::sample_v_given_h(int *h0_sample, vector<int>& imovie,double *mean, int *sample) {
	int i ;
  for(int ir=0; ir<imovie.size(); ir++) {
	i = imovie[ir];
		propdown(h0_sample, i, vbias + 5*i, mean + 5*ir);
		binomial5( mean+5*ir, sample+5*ir);
  }
}

double RBM::propup(int *v,vector<int> &imovie, double *w, double *d, double b) {

  double pre_sigmoid_activation = 0.0;

  for(int ns = 0; ns < 5; ns++)
  {
	for(int j=0; j<imovie.size(); j++) {
		pre_sigmoid_activation += w[5*imovie[j] + ns] * v[5*j + ns];
	}
  }

  for(int j=0; j<imovie.size(); j++) 
  {
		pre_sigmoid_activation += d[imovie[j]] ;
  }


  pre_sigmoid_activation += b;

  return sigmoid(pre_sigmoid_activation);
}

void RBM::propdown(int *h, int i, double* b, double *means) {
  double pre_sigmoid_activation[5] = {0.0};
//  cout <<"setp" << pre_sigmoid_activation[4] << endl;
  for(int ns = 0; ns < 5; ns++)
  {
	for(int j=0; j<n_hidden; j++) {
		pre_sigmoid_activation[ns] += W[j][5*i+ns] * h[j];
	}
	pre_sigmoid_activation[ns] += b[ns];
  }
	 //cout << pre_sigmoid_activation[0] <<" " << pre_sigmoid_activation[4] << endl;
	 sigmoid5(pre_sigmoid_activation,means);
	 return;
}

void RBM::gibbs_hvh(int *h0_sample,vector<int>& imovie, double *nv_means, int *nv_samples, \
                    double *nh_means, int *nh_samples) {
  sample_v_given_h(h0_sample, imovie, nv_means, nv_samples);
  sample_h_given_v(nv_samples, imovie,  nh_means, nh_samples);
}

void RBM::reconstruct(vector<int>& vv,vector<int>& imovie, vector<double>& reconstructed_v) {
  int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation[5] = {0.0};

  for(int i=0; i<n_hidden; i++) {
    h[i] = propup(v,imovie,  W[i],D[i], hbias[i]);
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
/*
void RBM::reconstruct_and_predict(
		vector<int>& vv,
		vector<int>& imovie, 
		vector<int>& imovie2,
		vector<int>& imovie3,
	   	vector<double>& reconstructed_v) {

  int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation[5];

  for(int i=0; i<n_hidden; i++) {
    h[i] = propup(v,imovie,  W[i], hbias[i]);
  }

  for(int i=0; i<imovie2.size(); i++) {
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] = 0.0;
	}
	  int im = imovie2[i];
    for(int j=0; j<n_hidden; j++) {
		for(int ns = 0; ns < 5; ns++)
		{
			pre_sigmoid_activation[ns] += W[j][5*im + ns] * h[j];
		}
    }
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] += vbias[5*im + ns];
	}

    sigmoid5(pre_sigmoid_activation, &(reconstructed_v[5*im]));
  }

  for(int i=0; i<imovie3.size(); i++) {
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] = 0.0;
	}
	  int im = imovie2[i];
    for(int j=0; j<n_hidden; j++) {
		for(int ns = 0; ns < 5; ns++)
		{
			pre_sigmoid_activation[ns] += W[j][5*im + ns] * h[j];
		}
    }
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] += vbias[5*im + ns];
	}

    sigmoid5(pre_sigmoid_activation, &(reconstructed_v[5*im]));
  }


  delete[] h;
}

*/


void RBM::reconstruct(vector<int>& vv,vector<int>& imovie, vector<int>& imovie2, vector<double>& reconstructed_v) {
  int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation[5];

  for(int i=0; i<n_hidden; i++) {
    h[i] = propup(v,imovie,  W[i],D[i], hbias[i]);
  }

  for(int i=0; i<imovie2.size(); i++) {
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] = 0.0;
	}
	  int im = imovie2[i];
    for(int j=0; j<n_hidden; j++) {
		for(int ns = 0; ns < 5; ns++)
		{
			pre_sigmoid_activation[ns] += W[j][5*im + ns] * h[j];
		}
    }
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] += vbias[5*im + ns];
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
	ifstream traindat("../../data/1234.dta");
	ifstream testdat("../../data/4.dta");
	ifstream testdata2("../../data/5.dta");
	vector<int> a;
	vector<int> b;
	vector<int> tmp(5,0);
	int iu = 1;

    while(traindat >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4])
    {
		if(iu < tmp[0])
		{
			train.push_back(a);
			train_index.push_back(b);
			uid_train.push_back(iu-1);
			a.clear();
			b.clear();
			iu = tmp[0];
		}
		for(int f = 1; f< 6; f++) a.push_back(tmp[4] == f); //push into the array a 5d vector
		b.push_back(tmp[1]-1);
    }
	train.push_back(a);
	train_index.push_back(b);
	uid_train.push_back(iu-1);
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
			//cout << iu << endl;
		}
		for(int f = 1; f< 6; f++) a.push_back(tmp[4] == f); //push into the array a 5d vector
		b.push_back(tmp[1]-1);
    }
	test.push_back(a);
	test_index.push_back(b);
	uid_test.push_back(iu -1);
	a.clear();
	b.clear();
	cout << "Test set loaded!" << test.size() <<"data" << endl;

	iu = -1;

    while(testdata2 >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4])
    {
		if(iu == -1) iu = tmp[0];
		if(iu < tmp[0])
		{
			uid_quiz.push_back(iu-1);
			quiz.push_back(a);
			quiz_index.push_back(b);
			a.clear();
			b.clear();
			iu = tmp[0];
			//cout << iu << endl;
		}
		for(int f = 1; f< 6; f++) a.push_back(tmp[4] == f); //push into the array a 5d vector
		b.push_back(tmp[1]-1);
    }
	quiz.push_back(a);
	quiz_index.push_back(b);
	uid_quiz.push_back(iu -1);
	a.clear();
	b.clear();
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
	double lambda = 0.001;

	RBM rbm(train_N,  n_visible, n_hidden, NULL, NULL, NULL);
	cout << "1111" << endl;
	int iu = 0;
//	int len = trainset.size();
	int len = train_N;

//    	 cout <<"error" << compute_loss(trainuid,testuid,trainset,train_index,testset,test_index,rbm) << endl;

	 std::thread threads[NUMTH];

	
	for(int iepo = 0 ; iepo < 100; iepo++)
	{
		
		for(int y = 0 ; y < len; y+=NUMTH)
		{
		if(y%1000 == 0) cout << y << endl;

			for(int it = 0; it < NUMTH; it++)
			{
				int iy = (y+it)%len;
				threads[it] = std::thread(&RBM::contrastive_divergence,&rbm,std::ref(trainset[iy]),std::ref(train_index[iy]), learning_rate,lambda, k);
//		rbm.contrastive_divergence(trainset[y], train_index[y], learning_rate, k);
			}
			for(auto& th:threads) th.join();

		}
    	 cout<<"error:" <<compute_loss(trainuid,testuid,trainset,train_index,testset,test_index,rbm,f1) << endl;
		 cout<<"error,quiz:" <<compute_loss(trainuid,quizuid,trainset,train_index,quizset,quiz_index,rbm,f2) << endl;

//	cout <<compute_loss_train(trainset,train_index,rbm) << endl;

		 learning_rate *= 0.91;
	}


/*
	int counter = 0;
	vector<double> reconstructed_X(5*NMOVIE,0.0);
	rbm.reconstruct(trainset[0], train_index[0] , reconstructed_X);
	for(int i = 0; i < train_index[0].size(); i++)
	{
	   cout << trainset[0][5*i] << "\t" << reconstructed_X[5*train_index[0][i]]<< endl;
	}

*/

	double maxv = -999.0,minv =999.0;
}


int main(int argc, char** argv) {
	string f1 = std::string(argv[1]);
	string f2 = std::string(argv[2]);
	cout << f1 << f2 <<endl;

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
	  test_rbm(f1,f2);
#ifdef MPIU
	MPI_Barrier(MPI_COMM_WORLD);
#endif
  return 0;
}


