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

vector<int> get_top_q_users (int q) {
     int a,b;
     int total = 0;
     vector<int> top_q_users;
     std::ifstream infile("../../data/sorted_users.dta");

     while(infile >> a >> b)
     {
          total += 1;
         // / printf("%i\n", a);
          top_q_users.push_back(a);
          if (total == q)
          {
               break;
          }
          
     }

     

     return top_q_users; // Indicates that everything went well.
}

int main(int argc, char const *argv[])
{
     vector <int> a = get_top_q_users(10);
     for (int i = 0; i < 10; i++)
     cout << a[i] << "\n";

     return 0;
}