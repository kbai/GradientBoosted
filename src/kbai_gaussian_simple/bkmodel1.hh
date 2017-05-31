struct bkmodel1: bkmodel
{
	vector<vector<vector<float>>> bst;
//	bool alternating;
	bkmodel1();
	virtual float g(int,int,int,int,float);
	virtual  void update_param_sgd(feature &a);
	virtual ~bkmodel1(){};
};
