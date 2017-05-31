#include "commonheader.h"
using namespace std;
int main() {
  ifstream infile(std::string(DATAPATH)+"1.dta");
  vector<int> tmp(5,0);
  int iu = 1;
  int nproc=50;

  ofstream out[nproc];

  for(int i =0; i < nproc; i++)
  {
	  out[i].open(std::string(DATAPATH)+"MPI"+std::to_string(i) +".dta");
  }

  while((infile >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4]))
  {
		int	  p = tmp[0]%nproc;
	  out[p] << tmp[0] <<"\t"<<tmp[1] <<"\t" << tmp[2] <<"\t" << tmp[3] << "\t" << tmp[4] << endl;
  }

for(int i =0; i < nproc; i++)
  {
	  out[i].close(); 
  }

  return 0;
}
