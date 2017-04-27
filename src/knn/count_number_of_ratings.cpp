#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>

main ( int argc, char **argv ) {
     int a,b,c,d;
     int old_b = 1;
     int old_a = 1;
     
     // double sum = 0.0;
     // double counter = 0.0;
     // double movies[17771];
     int users[458294] = {0.};
     std::ifstream infile("../../data/1.dta");
//     std::ifstream outputfile("../../qual.dta");
     
     while(infile >> a >> b >> c >> d)
     {
          printf("%i\n", a);
          users[a] += 1;
          
     }
     printf("%i\n", users[1]);

     std::ofstream myfile ("../../data/user_to_num_movies.dta");
     for(int i = 1; i < 458294; i++)
     {    
          myfile << i << ' ' << users[i] << "\n";
     }
     

     return 0; // Indicates that everything went well.
}