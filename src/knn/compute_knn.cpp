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

vector<int> get_top_q_users (int q) {
     int a,b;
     int total = 0;
     vector<int> top_q_users;
     std::ifstream infile("../../data/sorted_users.dta");

     while(infile >> a >> b)
     {
          total += 1;
          printf("%i\n", a);
          top_q_users.push_back(a);
          if (total == q)
          {
               break;
          }
          
     }

     return top_q_users; // Indicates that everything went well.
}

int main ( int argc, char **argv ) {
     int a,b,c,d,e;
     int q = 2000;
     int k =100;
     printf("hi\n");
     vector <int> top_q = get_top_q_users(q);


     vector<vector < int > > movies_users_rating (17770);
     vector<vector < int > > movies_users_id (17770);
     //std::ifstream infile("../../data/k_n_n.dta");
     std::ifstream infile("../../data/1.dta");
     
     printf("hi\n");
     while(infile >> a >> b >> c >> d >> e)
     {
          printf("%i\n", a);
          for(int i = 0; i < q; i++)
          {
               if(a == top_q[i])
               {
                    movies_users_id[b].push_back(a);
                    movies_users_rating[b].push_back(e);     
               }
                    
          }
          

     }
     printf("asdf\n");
     string line;
     std::vector<int> k_users(k);
     vector<vector<int> > k_users_all(458293, vector<int>(k));
     std::ifstream otherfile("../../data/richierozay.dta");
     while(getline(otherfile, line)) {
          int current_user = 1;

          istringstream iss(line);
          string word;
          int word_index = 0;
          while(getline(iss, word, '\t')){
               k_users[word_index] = stoi(word);
               word_index++;
          }
          for(int i = 0; i < k; i++)
          {
               k_users_all[current_user][i] = k_users[i];
          }
     }

     


     std::ifstream outputfile("mu/qual.dta");

     std::ofstream myfile ("first_try.dta");
     while(outputfile >> a >> b >> c)
     {   
         int total = 0;
         int counter = 0; 
         int user = a;
         int movie = b;
         k_users = k_users_all[k];
         for (int i; i < k; i++)
         {
               for (int j; j < movies_users_id[b].size(); j++)
               {
                    if(k_users[i] == movies_users_id[b][i])
                    {
                         total += movies_users_rating[b][i];
                         counter++;
                    }
               }
         }




          myfile << total/counter << "\n";
     }





     // std::ofstream myfile ("../../data/user_movie_rating_matrix.dta");
     // for(int i = 0; i < 17771; i++)
     // {    myfile << "Movie " << i << ":    ";
     //      for(int j = 0; j < movies_users_id[i].size(); j++)
     //      {
     //           myfile << "( " << movies_users_id[i][j] << ", " << movies_users_rating[i][j] << " ), ";     
     //      }
     //      myfile << "\n";
     // }
     

     return 0; // Indicates that everything went well.
}