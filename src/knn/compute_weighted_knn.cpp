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
     int q = 1500;
     int k = 100;
     printf("hi\n");
     vector <int> top_q = get_top_q_users(q);


     vector<vector < int > > movies_users_rating (17770);
     vector<vector < int > > movies_users_id (17770);
     //std::ifstream infile("../../data/k_n_n.dta");
     std::ifstream infile("../../data/1.dta");
     
     printf("hi\n");
     while(infile >> a >> b >> c >> d >> e)
     {
          if(a%50000 ==0)
          {
          printf("%i\n", a);
               
          }
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
     vector<vector<int> > k_users_all(458295, vector<int>(k));
     std::ifstream otherfile("../../data/k_nearest_neighbors.dta");
     int current_user = 1;

     while(getline(otherfile, line)) {
          

        if (current_user % 50000 == 0) {
            cout << "calculated " << current_user << " users" << endl;
        }


        
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

          current_user++;
     }

     cout << " finished calculating users" << endl;
     
     // Read in correlation matrix and store the weights

     std::ifstream corr("../../data/tanger.dta");
     std::vector<float> k_weights(q, 0);
     vector<vector<float>> k_users_weights(458295, vector<float>(q, 0));
     current_user = 1;
     while (getline(corr, line)) {

        if (current_user % 50000 == 0) {
            cout << "read " << current_user << " correlations" << endl;
        }

        istringstream iss(line);
        string word;
        int word_index = 0;
        while(getline(iss, word, '\t')) {
            k_weights[word_index] = stof(word);
            word_index++;
        }
        for (int i = 0; i < q; i++) {
             k_users_weights[current_user][i] = (k_weights[i] + 1) / 2;
        }

        current_user++;
     }
     
     printf("hi there\n");

     for (int i = 0; i < q; i++) {
        cout << k_users_weights[1][i] << endl;
     }

     std::ifstream outputfile("../../data/qual.dta");

     std::ofstream myfile ("weighted_knn.dta");
     while(outputfile >> a >> b >> c)
     {   
         float total = 0;
         float counter = 0; 
         int user = a;
         int movie = b;
         int curr_neighbor;
         float weight;
         k_users = k_users_all[a];
         for (int i = 0; i < k; i++)
         {
               for (int j = 0; j < movies_users_id[b].size(); j++)
               {
                    if(k_users[i] == movies_users_id[b][j])
                    {
                         curr_neighbor = k_users[i];
                         weight = k_users_weights[a][curr_neighbor];
                         total +=  weight * movies_users_rating[b][j];
                         counter += abs(weight);
                    }
               }
         }

         if(counter != 0)
         {
          myfile << total/counter << endl;
         }
         else
         {
          myfile << 3 << endl;
         }

          
     }





     std::ofstream mfile ("../../data/user_movie_rating_matrix.dta");
     for(int i = 0; i < 17771; i++)
     {    mfile << "Movie " << i << ":    ";
          for(int j = 0; j < movies_users_id[i].size(); j++)
          {
               mfile << "( " << movies_users_id[i][j] << ", " << movies_users_rating[i][j] << " ), ";     
          }
          mfile << "\n";
     }
     

     return 0; // Indicates that everything went well.
}