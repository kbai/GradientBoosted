#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <math.h>

using namespace std;



vector<vector < int > > movies_users_rating (17770);
vector<vector < int > > movies_users_id (17770);
vector<vector < double > > movies_users_time (17770);
// vector<double> global_movie (17770);
// vector<double> global_user (458294);
vector <double> min_time (458294, 2500);





int main(int argc, char **argv) {
     int a,b,c,d,e;
     int old_a = 1;
     double Xui;
     //std::ifstream infile("../../data/k_n_n.dta");
     std::ifstream infile("../../data/1.dta");
     
     printf("hi\n");
     // movies start at 1
     while(infile >> a >> b >> c >> d >> e)
     {
        movies_users_id[b].push_back(a);
        movies_users_rating[b].push_back(e);
        movies_users_time[b].push_back(c);
        if(min_time[a] > c)
        {
            min_time[a] = c;
        }
             
        
    }

    
     
    for(int i = 1; i < 17770; i++)
     {
        for(int j = 0; j < movies_users_id[i].size(); j++)
        {
            Xui = movies_users_time[i][j] - min_time[movies_users_id[i][j]];
            movies_users_time[i][j] = sqrt(Xui);
        }
     } 


    std::ofstream myfile2 ("../../data/user_times.dta");

    for(int i = 1; i < 17770; i++)
     {
        for(int j = 0; j < movies_users_id[i].size(); j++)
        {

            myfile2 << movies_users_id[i][j] << ' ' << movies_users_time[i][j] << endl;


        }
     } 




     return 0;
}