#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#define NMOVIE 17770
#define NUSER 458294
#define DATAPATH "../../data/"
#define MUDATAPATH "../../data/mu/"
#define TRAINFILE "3.dta"


using namespace std;

/* Creates a correlation matrix between all users and the top Q users 
 * Takes vector of user ids of top Q users, filename of data file, 
 * and filename to write matrix to. 
 */
void create_user_correlation(vector<int>& ids, std::string data_file, std::string matrix_file) { 
	cout << "Intializing..." << endl;
	size_t Q = ids.size();
	vector<float> means(NUSER, 0.0);
	vector<float> stdev(NUSER, 0.0);


	vector<vector<float>> cov(NUSER, vector<float>(Q, 0.0));
	vector<vector<float>> corr(NUSER, vector<float>(Q, 0.0));


	int u[6];
	int nratings = 0;
	int curr_user = 0;
	vector<float> ratings;
	cout << "Initialization successful!" << endl;

	// Load data from training file, write into file matrix
	ifstream infile(std::string(DATAPATH) + data_file);
	ofstream outfile(std::string(DATAPATH) + matrix_file);

	cout << "Reading..." << endl;
	// Read the file and populate array of means
	while (!infile.eof()) {
		infile >> u[0] >> u[1] >> u[2] >> u[3] >> u[4];
		// Scan user by user until file finished, populate means
		// ID numbers are indexed starting at 1, subtract 1 to match 0 indexing in C++
		while (u[0] - 1 > curr_user){
			// If user has not rated any movies, default average is 3
			// If no ratings, default the standard deviation to 1
			if (nratings == 0) {
				means[curr_user] = 3;
				stdev[curr_user] = 1;
			}
			else {
				means[curr_user] /= nratings;

				// Compute and store the standard deviation for the current user
				// If only 1 rating, default the standard deviation to a uniform stdev
				if (nratings == 1) {
					stdev[curr_user] = 4 / sqrt(12);
				}
				else {
					for (int i = 0; i < nratings; i++) {
						stdev[curr_user] += pow((ratings[i] - means[curr_user]), 2);
					}
				stdev[curr_user] = sqrt(stdev[curr_user] / nratings); 
				}
				

			}
			curr_user += 1;
			ratings.clear();
			nratings = 0;
		}
	
		means[curr_user] += u[4];
		ratings.push_back(u[4]);
		nratings += 1;

	}
	infile.close();
	cout << "Read successful!" << endl;


	cout << "Computing correlations..." << endl;
	// Using stored means and standard deviations, compute correlation coefficients
	// and divide each element by the product of standard deviations
	ifstream mufile(string(MUDATAPATH) + data_file);
	int curr_movie = 0;
	vector<int> users;
	ratings.clear();
	nratings = 0;
	vector<int> M(NUSER, 0);
	while (!mufile.eof()) {
		mufile >> u[0] >> u[1] >> u[2] >> u[3] >> u[4];
		// Check all the users who have watched a movie, and then update correlations
		// with running average
		while (u[1] - 1 > curr_movie) {
			for (int i = 0; i < nratings; i++) {
				int ind_i = users[i] - 1;
				for (size_t j = 0; j < Q; j++) {
					// If a user was in the top Q, then update column in correlations
					// with the running average
					if ((size_t) ind_i == j) {
						for (int k = 0; k < nratings; k++) {
							int ind_k = users[k] - 1;
							corr[ind_k][j] = (corr[ind_k][j] * M[ind_k] + (ratings[k] - means[ind_k]) * (ratings[i] - means[ind_i])) / (M[ind_k] + 1);
							M[ind_k] += 1;
						}
					}
				}
			}

			users.clear();
			ratings.clear();
			nratings = 0;
			curr_movie += 1;

			if (curr_movie % 100 == 0) {
				cout << curr_movie << endl;
			}
		}

		// Store all users who have watched a given movie
		users.push_back(u[0]);
		ratings.push_back(u[4]);
		nratings += 1;



	}

	for (int i = 0; i < NUSER; i++) {
		for (size_t j = 0; j < Q; j++) {
			corr[i][j] /= (stdev[i] * stdev[j]);
		}
	}
	cout << "Correlations computed!" << endl;

	cout << "Writing..." << endl;
	// Write correlation matrix to file 
	
	for (int i = 0; i < NUSER; i++) {
		std::stringstream ss;
		for(size_t j = 0; j < Q; ++j) {
		 	if(j != 0)
				ss << "\t";
			ss << to_string(corr[i][j]);
		}
		std::string s = ss.str();
		outfile << s << endl;
	} 
	outfile.close();
	cout << "Write successful!" << endl;

	cout << "DEBUG: Writing means and stdev..." << endl;
	ofstream meanfile(string(DATAPATH) + "means.dta");
	ofstream stdevfile(string(DATAPATH) + "stdevs.dta");
	for (int i = 0; i < NUSER; i++) {
		meanfile << to_string(means[i]) << endl;
		stdevfile << to_string(stdev[i]) << endl;
	}
	meanfile.close();
	stdevfile.close();

}

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


int main(int argc, char ** argv) {
	vector<int> Q = get_top_q_users(10);
	create_user_correlation(Q, "1.dta", "tanger.dta");
	return 0;
}