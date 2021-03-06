struct bkmodel_gpu: bkmodel
{
	float* d_bm;
	float* d_bu;
	float* d_pm;
	float* d_pu;
	float* d_pu1;
	float* d_ptu;
	float* d_bt;
	float* d_bta;
	float* d_btu;
	float* d_btm;
	float* d_bf;
	float* d_bfm;
// gpu data storage:
	int* d_iu;
	int* d_im;
	int* d_it;
	int* d_if;
	float* d_tb;
	int* d_rate;
	int* d_iu1;
	int* d_im1;
	int* d_it1;
	int* d_if1;
	float* d_tb1;
	int* d_rate1;
	int* d_iu2;
	int* d_im2;
	int* d_it2;
	int* d_if2;
	float* d_tb2;
	int* d_rate2;
    int sz0;
	int sz1;
	int sz2;
	int beg=0;

//	bool alternating;
	bkmodel_gpu();
	void test(float lr);
	void retrieve_gpu();
	void loaddata(feature &a, feature &a1, feature &a2);
	double compute_error();
	virtual float g(int a,int b,int c,int d,float e){return bkmodel::g(a,b,c,d,e);};
	virtual  void update_param_sgd(feature &a){return;};
	virtual ~bkmodel_gpu(){};
};
