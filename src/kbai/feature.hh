struct bkmodel;
struct feature
{
	int _iu;
	int _im;
	int _it;
	int _if;
	int _rate;
	int _ita;
	float _tb;
	float* residue;

	vector<int> viu;
	vector<int> vim;
	vector<int> vit;
	vector<int> vita;
	vector<int> vif;
	vector<int> vrate;
	vector<float> vtb;


	bool _loadall;
	vector<vector<int>> _testset;
	vector<float> tu;
	ifstream & _ifile;
	feature(ifstream& a, bool b=false);
	void retrieve_feature();
	void load_alldata();
	void load_gpu();
	double compute_RMSE(bkmodel &bk);
	double compute_QUAL(bkmodel &bk, string ofilename);
    void compute_residue(bkmodel &abk);
};


