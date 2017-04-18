#include "commonheader.h"
using namespace std;
#include "feature.hh"
#include "bkmodel.hh"

feature::feature(ifstream &a):_ifile(a),loadall(false)
{}

feature::feature(ifstream &a, bool b):_ifile(a),loadall(b)
{
	if(loadall) load_alldata();
}


void feature::retrieve_feature()
{
	int a;
	if(!(_ifile >> iu  >> im  >> it >> a>> rate))
	{
		cout << "here is the end" << endl;
		_ifile.clear();
		_ifile.seekg(0,ios::beg);// if read fails go to beginning
		_ifile >> iu  >> im  >> it >> a>> rate;
	}
}

void feature::load_alldata()
{
	loadall = true;
	vector<int> data(5,0);
	_ifile.clear();
	_ifile.seekg(0,ios::beg); 
	while(_ifile>>data[0]>>data[1]>>data[2]>>data[3]>>data[4])
	{
		testset.push_back(data);
	}
	cout << testset.size() << endl;
}


double feature::compute_RMSE(bkmodel &bk)
{
	if(!loadall)
	{
		cout << "cannot compute RMSE" << endl;
		return -1.0;
	}
	int im,iu,rate,it;
	int n = testset.size();
	double rmse = 0.0;
	for(auto & p: testset)
	{
		iu = p[0]-1;
		im = p[1]-1;
		it = p[2];
		rate = p[4];
		rmse +=pow( (bk.g(iu,im,it) - rate),2.0)/n;
	}
	rmse = sqrt(rmse);
	return rmse;
}


double feature::compute_QUAL(bkmodel &bk, string ofilename)
{
	if(!loadall)
	{
		cout << "cannot compute RMSE" << endl;
		return -1.0;
	}
	int im,it,iu,rate;
	int n = testset.size();
	double rmse = 0.0;
	ofstream qual(ofilename);
	for(auto & p: testset)
	{
		iu = p[0]-1;
		im = p[1]-1;
		it = p[2];
		rate = p[4];
		qual << bk.g(iu,im,it) << endl;
	}
	rmse = sqrt(rmse);
	qual.close();
	return rmse;
}






