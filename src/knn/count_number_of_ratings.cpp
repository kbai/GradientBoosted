#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
using namespace std;

main ( int argc, char **argv ) {
     int a,b,c,d,e;
  
     // double sum = 0.0;
     // double counter = 0.0;
     // double movies[17771];
     int users[458294] = {0};
     std::ifstream infile("../../data/1.dta");
//     std::ifstream outputfile("../../qual.dta");
     

     while(infile >> a >> b >> c >> d >> e)
     {
     if (a % 1000 == 0) {
          cout << a << endl;
     }

          users[a] += 1;

     }


     std::ofstream myfile ("../../data/user_to_num_movies.dta");
     for(int i = 1; i < 458294; i++)
     {    
          myfile << i << "\t" << users[i] << endl;
     }
     

     return 0; // Indicates that everything went well.
}