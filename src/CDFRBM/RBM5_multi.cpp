#include "commonheader.h"
#include <random>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <cstring>
#include "utils.h"
#include <stdlib.h>
//#include <mutex>
#define NUMTH 50
#include <thread>
#include <mutex>
using namespace std;

#ifdef MPIU
#endif
#include "RBM5.h"
using namespace utils;
default_random_engine e(0);


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
//		cout << "cl1" << endl;
		testrbm.reconstruct(train[i],train_index[i],reconstructed_X);
//		cout << "cl2" << endl;
		counter += train_index[i].size();
		for(int j = 0; j < train_index[i].size(); j++)
		{
			double score = 0.0;
			double tscore = 0.0;
			for(int ns = 0; ns < 5; ns++)
			{
				score += ( ns + 1 )* (reconstructed_X[5*train_index[i][j]+ns]);
				tscore += train[i][5*j+ns] * (ns+1);
//				cout << 5*train_index[i][j] + ns << "  "<< 5*NMOVIE << endl;
			}
//			cout << score <<" : " << tscore << endl;
			error += pow(score-tscore,2.0);
		}
//		cout << "cl3" << endl;

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
//		assert((idmap[i] == i));
//		assert((idmap[i] <= train.size()));
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
	cout << "error" << error << endl;
	output.close();
	return error;
}

RBM::RBM(int size,int nlat, int n_v, int n_h, double **w, double *hb, double *vb) {
  N = size;
  n_lat = nlat;
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
	A = new double*[n_hidden];
	B = new double*[n_lat];
	dW = new double*[n_hidden];
	sW = new double*[n_hidden];
	for(int i=0; i < n_lat; i++)
		B[i] = (double*) malloc(5*sizeof(double)*n_visible);

    for(int i=0; i<n_hidden; i++)
	{
		try{

		A[i] = (double*) malloc(n_lat*sizeof(double));
		}
		catch(std::bad_alloc&)
		{
			cout << "bad allocation error!" << endl;
		}
	}
    double a = 1.0 / n_visible;

    for(int i=0; i<n_hidden; i++) 
	{




	  for(int j = 0; j < n_lat; j++) 
		  A[i][j] = uniform(-a,a);
    }

  for(int i = 0; i < n_lat; i++)
	  for(int j =0; j < 5*n_visible; j++) 
		  B[i][j] = uniform(-a,a);
  
	
    hbias = new double[n_hidden];
    for(int i=0; i<n_hidden; i++){
	   	hbias[i] = 0;
	}

    vbias = new double[5*n_visible];
    for(int i=0; i< 5 * n_visible; i++) {
		vbias[i] = 0;
	}
}
}
void RBM::output() 
{
  cout << "W:"<<A[e()%n_hidden][e()%n_lat] <<"\t hbias:" << hbias[e()%n_hidden] <<"\t vbias:" << 
	  vbias[e()%(5*n_visible)] << endl;


}


RBM::~RBM() {
  for(int i=0; i<n_hidden; i++) delete[] A[i];
  delete[] A;
  delete[] hbias;
  delete[] vbias;
}
void RBM::contrastive_divergence(std::vector<int>& input ,std::vector<int>& imovie, double lr, int k) {
  double *ph_mean = new double[n_hidden];
  int *ph_sample = new int[n_hidden];
  double *nv_means = new double[5*n_visible];
  int *nv_samples = new int[5*n_visible];
  double *nh_means = new double[n_hidden];
  int *nh_samples = new int[n_hidden];
  double *pinter = new double[n_lat];
  double *ninter = new double[n_lat];
  double *pinterb = new double[n_lat];
  double *ninterb = new double[n_lat];


  /* CD-k */
  sample_h_given_v(&(input[0]), imovie, ph_mean, ph_sample);

  for(int step=0; step<k; step++) {
    if(step == 0) {
      gibbs_hvh(ph_sample,imovie, nv_means, nv_samples, nh_means, nh_samples);
    } else {
      gibbs_hvh(nh_samples,imovie,  nv_means, nv_samples, nh_means, nh_samples);
    }
  }

  //cout << "cccccc"<< endl;
  for(int i=0; i<n_lat; i++) 
  {
	  pinter[i] = 0;
	  ninter[i] = 0;
    for(int j=0; j< n_hidden; j++) 
	{
		pinter[i] += ph_mean[j] * A[j][i];
		ninter[i] += nh_means[j] * A[j][i];

    }
	
  }

  for(int i=0; i<n_lat; i++) 
  {
	  pinterb[i] = 0;
	  ninterb[i] = 0;
    for(int jr=0; jr<imovie.size(); jr++) 
	{
		int j = imovie[jr];
		for(int ns = 0; ns <5 ; ns++)
		{
			pinterb[i] += input[5*jr+ns] * B[i][5*j+ns];
			ninterb[i] += nv_means[5*jr+ns] * B[i][5*j+ns];
		}

    }
	
  }
 
  mtx.lock();

  for(int i=0; i<n_lat; i++) 
  {
    for(int jr=0; jr<imovie.size(); jr++) 
	{
		int j = imovie[jr];
		for(int ns = 0; ns < 5; ns ++)
		{
			B[i][5*j + ns] +=  lr * (pinter[i] * input[5*jr + ns] - ninter[i] * nv_means[5*jr+ns] - 1E-3*B[i][5*j + ns]);
		}
    }
    for(int j=0; j<n_hidden; j++) 
	{
		A[j][i] += 0.1 * lr * (pinterb[i] * ph_mean[j] - ninterb[i] * nh_means[j] - 1E-4*A[j][i]);
    }
  }
     for(int jr=0; jr<imovie.size(); jr++) 
	{
		int j = imovie[jr];
		for(int ns = 0; ns < 5; ns ++)
		{
			vbias[5*j + ns] += 4.0 * lr * (input[5*jr + ns] - nv_means[5*jr+ns] - 1E-3*vbias[5*j + ns]);
		}
    }
    for(int j=0; j<n_hidden; j++) 
	{
		hbias[j] += 0.25 * lr * ( ph_mean[j] - nh_means[j] - 1E-5*hbias[j]);
    }
 
  
  mtx.unlock();

// cout << "cd3"<< endl; 
  delete[] ph_mean;
  delete[] ph_sample;
  delete[] nv_means;
  delete[] nv_samples;
  delete[] nh_means;
  delete[] nh_samples;
  delete[] pinterb;
  delete[] ninterb;
  delete[] pinter;
  delete[] ninter;
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

void RBM::sample_h_given_v(int *v0_sample, vector<int>& imovie,double *mean, int *sample) 
{
  double inter[n_lat];
  for(int i = 0; i < n_lat; i++) inter[i] = 0.0;
  for(int i = 0 ; i < n_lat; i++) 
	for(int ns = 0; ns < 5; ns++)
		for(int j=0; j<imovie.size(); j++) 
		{
			inter[i] += B[i][5*imovie[j] + ns] * v0_sample[5*j + ns];
		}

  for(int i=0; i<n_hidden; i++) 
  {
    mean[i] = propup(inter, A[i], hbias[i]);
    sample[i] = binomial(1, mean[i]);
  }
}

void RBM::sample_v_given_h(int *h0_sample, vector<int>& imovie,double *mean, int *sample) {
	int i ;
  double inter[n_lat];

  for(int i = 0; i < n_lat; i++) inter[i] = 0.0;

  for(int i = 0 ; i < n_lat; i++) 
		for(int j=0; j<n_hidden; j++) 
			inter[i] += A[j][i] * h0_sample[j];

  //cout << "kkkkkkkk"<< endl;

  for(int ir=0; ir<imovie.size(); ir++) 
  {
	  //cout <<imovie.size() << endl;
		i = imovie[ir];
		//cout << i << endl;
	    propdown(inter,i,vbias+5*i, mean+5*ir);
		//cout << "11111"<< endl;
		binomial5( mean+5*ir, sample+5*ir);
  }

}

double RBM::propup(double *inter, double *A, double b) {

  double pre_sigmoid_activation = 0.0;

	for(int j=0; j < n_lat ; j++) pre_sigmoid_activation += A[j] * inter[j];
  pre_sigmoid_activation += b;

  return sigmoid(pre_sigmoid_activation);
}

void RBM::propdown(double *h, int i, double* b, double *means) {
  double pre_sigmoid_activation[5] = {0.0};

//  cout <<"setp" << pre_sigmoid_activation[4] << endl;
  for(int ns = 0; ns < 5; ns++)
  {
	for(int j=0; j<n_lat; j++) {
		pre_sigmoid_activation[ns] += B[j][5*i+ns] * h[j];
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
  ///cout <<"ddddaaaaa"<< endl;
  sample_h_given_v(nv_samples, imovie,  nh_means, nh_samples);
}

void RBM::reconstruct(vector<int>& vv,vector<int>& imovie, vector<double>& reconstructed_v) {
   
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
void RBM::operator=(const RBM& a) {
	cout << "copy constructor called!" << endl;
  N = a.N;
  n_visible = a.n_visible;
  n_hidden = a.n_hidden;
  n_lat = a.n_lat;
    for(int i=0; i<n_hidden; i++)
	{
		memcpy(A[i],a.A[i],sizeof(double)*n_lat);
	}
	for(int i=0; i<n_lat; i++)
	{
		memcpy(B[i],a.B[i],sizeof(double)*5*n_visible);
	}
	
	memcpy(hbias	,a.hbias	,n_hidden*sizeof(double));
	memcpy(vbias	,a.vbias	,5*n_visible*sizeof(double));

}

void RBM::reconstruct(vector<int>& vv,vector<int>& imovie, vector<int>& imovie2, vector<double>& reconstructed_v) {
  int * v = &(vv[0]);
  double *h = new double[n_hidden];
  double pre_sigmoid_activation[5];
  double inter[n_lat];

  for(int i = 0; i < n_lat; i++) inter[i] = 0.0;
  for(int i = 0 ; i < n_lat; i++) 
	for(int ns = 0; ns < 5; ns++)
		for(int j=0; j<imovie.size(); j++) 
			inter[i] += B[i][5*imovie[j] + ns] * v[5*j + ns];

  //compute B*v
 
  for(int i=0; i<n_hidden; i++) {
    h[i] = propup(inter,  A[i], hbias[i]);
  }
  //compute sigmoid(A*B*v)

  for(int i = 0; i < n_lat; i++) inter[i] = 0.0;

  for(int i = 0 ; i < n_lat; i++) 
		for(int j=0; j<n_hidden; j++) 
			inter[i] += A[j][i] * h[j];

  //compute A^T*h

  for(int i=0; i<imovie2.size(); i++) 
  {
	for(int ns = 0; ns < 5; ns++)
	{
		pre_sigmoid_activation[ns] = 0.0;
	}
	  int im = imovie2[i];
    for(int j=0; j<n_lat; j++) 
	{
		for(int ns = 0; ns < 5; ns++)
		{
			pre_sigmoid_activation[ns] += B[j][5*im + ns] * inter[j];
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
	ifstream traindat("../../data/1.dta");
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
			//cout << iu << endl;
			//assert(b.size() < NMOVIE);
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
	int n_lat = 20;
	int n_hidden = 100;
	int training_epochs = 10000;
	int k = 1;
	vector<int> testuid, trainuid, quizuid;
	vector<vector<int>> trainset,testset,train_index, test_index;
	vector<vector<int>> quizset, quiz_index;
	load_all_data(trainset, train_index,  trainuid, testset, test_index, testuid, quizset, quiz_index, quizuid);
	int train_N = trainset.size();
	double learning_rate = 0.01;

	RBM rbm(train_N, n_lat,  n_visible, n_hidden, NULL, NULL, NULL);
	RBM rbm2(train_N, n_lat,  n_visible, n_hidden, NULL, NULL, NULL);

	cout << "1111" << endl;
	int iu = 0;
//	int len = trainset.size();
	int len = train_N;

//    	 cout <<"error" << compute_loss(trainuid,testuid,trainset,train_index,testset,test_index,rbm) << endl;

	 std::thread threads[NUMTH];
	 std::thread writedata[2];


	
	for(int iepo = 0 ; iepo < 30; iepo++)
	{
		if(k < 7)	k = k + 1;
			for(int y = 0 ; y < len; y+=NUMTH)
			{
			//cout << y<< endl;
			if(y%10000 == 0) 
			{
				cout << y << endl;
				rbm.output();
			}
	
			for(int it = 0; it < NUMTH; it++)
			{
				
				int iy = y+it;
				//cout <<"iy:"<< iy << endl;
				if(iy < len) threads[it] = std::thread(&RBM::contrastive_divergence,&rbm,std::ref(trainset[iy]),std::ref(train_index[iy]), learning_rate,k);
			}
			for(auto& th:threads) if(th.joinable()) th.join();
		}
		if(learning_rate > 1e-3) learning_rate *= 0.91;
		 if(writedata[0].joinable()) writedata[0].join();
		 if(writedata[1].joinable()) writedata[1].join();
		 rbm2 = rbm;
		 writedata[0] = std::thread(compute_loss,std::ref(trainuid),std::ref(testuid),std::ref(trainset),std::ref(train_index),std::ref(testset),std::ref(test_index),std::ref(rbm2),std::ref(f1));
		 writedata[1] = std::thread(compute_loss,std::ref(trainuid),std::ref(quizuid),std::ref(trainset),std::ref(train_index),std::ref(quizset),std::ref(quiz_index),std::ref(rbm2),std::ref(f2));
	}
		 if(writedata[0].joinable()) writedata[0].join();
		 if(writedata[1].joinable()) writedata[1].join();
	

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


