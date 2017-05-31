#include <iostream>
#include <fstream>
// decompose file into 5 parts
using namespace std;
typedef int Number;
int main()
{
	Number u[4],counter[5];
	int a = 0;
	int idx  = 0;
	streampos begin,end;
	ifstream myfile("./all.dta");
	ifstream id("./all.idx");
	ofstream out1("./1.dta");//,ios::binary);
	ofstream out2("./2.dta");//,ios::binary);
	ofstream out3("./3.dta");//,ios::binary);
	ofstream out4("./4.dta");//,ios::binary);
	ofstream out5("./5.dta");//,ios::binary);
	ofstream *p ;
	
	cout << sizeof(Number) << endl;
	
	
	
	while(!myfile.eof())
	{
		cout << a++ << endl;
		myfile >>u[0]>>u[1] >> u[2] >> u[3];
		id >> idx;
		switch(idx)
		{
			case 1:  p = &out1; counter[0]++; break;
			case 2:  p = &out2; counter[1]++; break;
			case 3:  p = &out3; counter[2]++; break;
			case 4:  p = &out4; counter[3]++; break;
			case 5:  p = &out5; counter[4]++; break;
		}	
		//p->write(reinterpret_cast<const char *>(u),4*sizeof(Number));
		(*p) << u[0] <<"\t" << u[1] << "\t" << u[2] << "\t" << u[3] << endl;
	}
	cout << a << endl;
	cout << counter[0]<< "\t" << counter[1]<<"\t" << counter[2]<<"\t" << counter[3]<<"\t" << counter[4] << endl;
	myfile.close();	
}
