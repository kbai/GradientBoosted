#define NLAT 200
#define MEAN 3.6095162 
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
	vector<vector<float>> ptu;
	vector<float> au;
	vector<float> htu;
	vector<float> bt;
	vector<float> bta;// time drifting term containing 2243 days
	vector<float> btu;
	vector<vector<float>> btm;
	vector<float> bf;
	vector<float> cut;
	vector<float> but;
	vector<vector<float>> bfm;
	vector<vector<float>> pmf;
	float mean;
	float _lr;
	float _alpha;
	float _lambda;

	bkmodel();
	void half_lr(); // setting learning rate
	virtual float g(int iu,int itu, int im, int it,int ife,float tt);
	virtual void update_param_sgd(feature &a);//no longer maintained.
	virtual ~bkmodel(){}; // have to have a virtual destructor
};


