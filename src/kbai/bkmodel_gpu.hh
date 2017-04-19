struct bkmodel_gpu: bkmodel
{
	float* d_bm;
	float* d_bu;
	float* d_bt;
	float* d_btu;
	float* d_btm;
//	bool alternating;
	bkmodel_gpu();
	void test();
	virtual float g(int,int,int){return 0;};
	virtual  void update_param_sgd(feature &a){return;};
	virtual ~bkmodel_gpu(){};
};
