struct bkmodel
{
	vector<int> uuser;
	vector<int> umovie;
	vector<vector<float>> bm; // 50 latent factors
	vector<float> bu;
	vector<vector<float>> bs; // 50 latent factors
	vector<float> btu;
	vector<vector<float>> bt;
	float mean;
	float _lr;
	float _alpha;
	float _lambda;

	bkmodel();
	void half_lr(); // setting learning rate
	virtual float g(int iu, int im, int it);
	void update_param_sgd(feature &a);
	virtual ~bkmodel(){}; // have to have a virtual destructor
};


