#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <string.h>
#include <sstream>

using namespace std;



vector<vector < int > > movies_users_rating (17770);
vector<vector < int > > movies_users_id (17770);
vector<vector < int > > movies_users_time (17770);
vector<double> global_movie (17770);
vector<double> global_user (458294);
vector<double> num_movie (17770);
vector<double> num_user (458294);







int main(int argc, char **argv) {
     int a,b,c,d,e;
     int old_a = 1;
     double sum = 0;
     double summs = 0;
     double counter = 0;
     double counter2 = 0;

     double global_mean;
     //std::ifstream infile("../../data/k_n_n.dta");
     std::ifstream infile("../../data/1.dta");
     
     printf("hi\n");
     // movies start at 1
     while(infile >> a >> b >> c >> d >> e)
     {
        movies_users_id[b].push_back(a);
        movies_users_rating[b].push_back(e);
        movies_users_time[b].push_back(c);
        summs += e;
        counter++;
        if (a == old_a)
        {      
            sum += e;
            counter2 += 1;   
            
        }
        else
          {
               global_user[old_a] = sum/counter2;
               num_user[old_a] = counter2;
               sum = e;
               counter2 = 1;
               old_a = a;
          }
    }

    printf("asdfasdf\n");
     for(int i = 1; i < 17770; i++)
     {
        sum = 0;
        counter2 = 0;

          for(int j = 0; j < movies_users_rating[i].size(); j++)
          {
            sum += movies_users_rating[i][j];
            counter2++;
          }
          global_movie[i] = sum/counter2;
          num_movie[i] = counter2;

     }
     global_mean = summs/counter;
     cout << global_mean << endl;
     cout << global_movie[1] << endl;
     cout << global_user[1] << endl;
     


     std::ofstream myfile ("../../data/user_means.dta");
     

     std::ofstream myfile2 ("../../data/movie_means.dta");



     std::ofstream myfile3 ("../../data/movie_nums.dta");


     std::ofstream myfile4 ("../../data/user_nums.dta");

     for(int i = 0; i < global_user.size(); i++)
     {   
        myfile << global_user[i] << endl;          
     }

     for(int i = 0; i < global_movie.size(); i++)
     {   
        myfile2 << global_movie[i] << endl;          
     }

     for(int i = 0; i < num_movie.size(); i++)
     {   
        myfile3 << num_movie[i] << endl;          
     }


     for(int i = 0; i < num_user.size(); i++)
     {   
        myfile4 << num_user[i] << endl;          
     }


     return 0;
}