#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
using namespace std;


bool sortcol(const vector<int>& v1, const vector<int>& v2)
{
     return v1[1] > v2[1];
}

int main ( int argc, char **argv ) {
     int a,b;
     //int old_b = 1;
     //int old_a = 1;
     
     // double sum = 0.0;
     // double counter = 0.0;
     // double movies[17771];
     //int users[458294] = {0.};
     vector <vector<int> > user_count;
     std::ifstream infile("../../data/user_to_num_movies.dta");
//     std::ifstream outputfile("../../qual.dta");

     vector <int> temp;
     while(infile >> a >> b)
     {
          //printf("%i\n", a);
          temp.push_back(a);
          temp.push_back(b);
          user_count.push_back(temp);
          temp.clear();
          
     }

     sort(user_count.begin(), user_count.end(), sortcol);

     std::ofstream myfile ("../../data/sorted_users.dta");
     for(int i = 0; i < 458293; i++)
     {    
          myfile << user_count[i][0] << ' ' << user_count[i][1] << "\n";
     }
     

     return 0; // Indicates that everything went well.
}