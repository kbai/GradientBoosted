#define NLAT 100
struct bkmodel
{
	bool alternating;
	vector<int> uuser;
	vector<int> umovie;
	vector<float> bm;
	vector<float> bu;
	vector<vector<float>> pm; // 50 latent factors
	vector<vector<float>> pu;
	vector<vector<float>> pu1;
	vector<float> bt;
	vector<float> btu;
	vector<vector<float>> btm;
	vector<float> bf;
	vector<vector<float>> bfm;
	float mean;
	float _lr;
	float _alpha;
	float _lambda;

	bkmodel();
	void half_lr(); // setting learning rate
	virtual float g(int iu, int im, int it,int ife,float tt);
	virtual void update_param_sgd(feature &a);
	virtual ~bkmodel(){}; // have to have a virtual destructor
};


