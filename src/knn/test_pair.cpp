#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
using namespace std;
int main()
{
	vector< pair<int,float> > a(5,std::make_pair(0,0.0));
	for(int i = 0; i < 5;i++)
	{
		a[i].first = i;
		a[i].second = 0.6-0.1*i;
	}
	std::sort(a.begin(), a.end(), [](pair<int,float>&a, pair<int,float>&b){return a.second < b.second;});
	for(int i = 0; i < 5; i++) cout << a[i].first << "\t" << a[i].second <<"\t"<< endl;
}
