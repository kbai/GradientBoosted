#include <assert.h>
#include "commonheader.h"
using namespace std;
#include "feature.hh"
#include "bkmodel.hh"



feature::feature(ifstream &a,bool b):_ifile(a),_loadall(b),tu(NUSER,0)
{
	ifstream tufile(std::string(DATAPATH)+"./averagetime.dta");
	for(auto &i : tu) tufile>>i;
	tufile.close();
	if(_loadall) load_alldata();
	srand(0);
}
int f(int a)
{
	int o =  (int)(log(a)/log(2.718));
	o = (o>=7)?7:o;
	return o; 
}

void feature::retrieve_feature()
{
	int a;
	if(_loadall)
	{
		int i = rand()%_testset.size();
		_iu = _testset[i][0];
		_im = _testset[i][1];
		_it = (int)(_testset[i][2]/75);
		assert(_it<30);
		_if = f(_testset[i][3]);
		_rate = _testset[i][4];
		_tb = _testset[i][2]-tu[_iu-1];
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
		_tb = _it - tu[_iu-1];
		_it = _it/75;
		_if = f(a);
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
void feature::load_gpu()
{
	for(auto & p: _testset)
	{
		_iu = p[0]-1;
		_im = p[1]-1;
		_it = (int)(p[2]/75);
		_if = f(p[3]);
		_tb = p[2] - tu[_iu];
		_rate = p[4];
		viu.push_back(_iu);
		vim.push_back(_im);
		vit.push_back(_it);
		vif.push_back(_if);
		vrate.push_back(_rate);
		vtb.push_back(_tb);
	}
	return;
}

double feature::compute_RMSE(bkmodel &bk)
{
	if(!_loadall)
	{
		cout << "cannot compute RMSE" << endl;
		return -1.0;
	}
	int _im,_iu,_rate,_it,_tb;
	int n = _testset.size();
	double rmse = 0.0;
	for(auto & p: _testset)
	{
		_iu = p[0]-1;
		_im = p[1]-1;
		_it = (int)(p[2]/75);
		_if = f(p[3]);
		_tb = p[2] - tu[_iu];
		_rate = p[4];
		rmse +=pow( (bk.g(_iu,_im,_it,_if,_tb) - _rate),2.0)/n;

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
	int _im,_it,_iu,_if,_rate;
	int n = _testset.size();
	double rmse = 0.0;
	ofstream qual(ofilename);
	for(auto & p: _testset)
	{
		_iu = p[0]-1;
		_im = p[1]-1;
		_it = (int)(p[2]/75);
		_if = f(p[3]);
		_tb = p[2] - tu[_iu];
		_rate = p[4];
		qual << bk.g(_iu,_im,_it,_if,_tb) << endl;
	}
	rmse = sqrt(rmse);
	qual.close();
	return rmse;
}






