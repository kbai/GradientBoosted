struct bkmodel;
struct feature
{
	int iu;
	int im;
	int it;
	int rate;
	bool loadall;
	vector<vector<int>> testset;
	ifstream & _ifile;
	feature(ifstream& a);
	feature(ifstream& a, bool b);
	void retrieve_feature();
	void load_alldata();
	double compute_RMSE(bkmodel &bk);
	double compute_QUAL(bkmodel &bk, string ofilename);

};


