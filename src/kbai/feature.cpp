#include "commonheader.h"
using namespace std;
#include "feature.hh"
#include "bkmodel.hh"

feature::feature(ifstream &a):_ifile(a),_loadall(false)
{}

feature::feature(ifstream &a, bool b):_ifile(a),_loadall(b)
{
	if(_loadall) load_alldata();
	srand(0);
}


void feature::retrieve_feature()
{
	int a;
	if(_loadall)
	{
		int i = rand()%_testset.size();
		_iu = _testset[i][0];
		_im = _testset[i][1];
		_it = _testset[i][2];
		a = _testset[i][3];
		_rate = _testset[i][4];
	}
	else // streaming mode
	{
		if(!(_ifile >> _iu  >> _im  >> _it >> a>> _rate))
		{
			cout << "here is the end" << endl;
			_ifile.clear();
			_ifile.seekg(0,ios::beg);// if read fails go to beginning
			_ifile >> _iu  >> _im  >> _it >> a>> _rate;
		}
	}
}

void feature::load_alldata()
{
	_loadall = true;
	vector<int> data(5,0);
	_ifile.clear();
	_ifile.seekg(0,ios::beg); 
	while(_ifile>>data[0]>>data[1]>>data[2]>>data[3]>>data[4])
	{
		_testset.push_back(data);
	}
	cout << _testset.size() << endl;
}


double feature::compute_RMSE(bkmodel &bk)
{
	if(!_loadall)
	{
		cout << "cannot compute RMSE" << endl;
		return -1.0;
	}
	int _im,_iu,_rate,_it;
	int n = _testset.size();
	double rmse = 0.0;
	for(auto & p: _testset)
	{
		_iu = p[0]-1;
		_im = p[1]-1;
		_it = p[2];
		_rate = p[4];
		rmse +=pow( (bk.g(_iu,_im,_it) - _rate),2.0)/n;
	}
	rmse = sqrt(rmse);
	return rmse;
}


double feature::compute_QUAL(bkmodel &bk, string ofilename)
{
	if(!_loadall)
	{
		cout << "cannot compute RMSE" << endl;
		return -1.0;
	}
	int _im,_it,_iu,_rate;
	int n = _testset.size();
	double rmse = 0.0;
	ofstream qual(ofilename);
	for(auto & p: _testset)
	{
		_iu = p[0]-1;
		_im = p[1]-1;
		_it = p[2];
		_rate = p[4];
		qual << bk.g(_iu,_im,_it) << endl;
	}
	rmse = sqrt(rmse);
	qual.close();
	return rmse;
}






