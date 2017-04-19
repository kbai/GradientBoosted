struct bkmodel;
struct feature
{
	int _iu;
	int _im;
	int _it;
	int _rate;
	bool _loadall;
	vector<vector<int>> _testset;
	ifstream & _ifile;
	feature(ifstream& a);
	feature(ifstream& a, bool b);
	void retrieve_feature();
	void load_alldata();
	double compute_RMSE(bkmodel &bk);
	double compute_QUAL(bkmodel &bk, string ofilename);

};


