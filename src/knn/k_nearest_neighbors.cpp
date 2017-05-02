#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#define NMOVIE 17770
#define NUSERS 458293

using namespace std;


int main(int argc, char **argv) {

    int Q = atoi(argv[1]);    // desired "top" users
    int K = atoi(argv[2]);    // desired nearest neighbors
    // reading in to find Q top users
    ifstream user_file("../../data/user_file.dta");
    int user, ratings;
    int current_user = 0;
    vector<int> top_users(Q);
    while (!user_file.eof()) {
        if (current_user >= Q) {
            break;
        }
        user_file >> user >> ratings;
        top_users[current_user] = user;
        current_user++;
    }
    user_file.close();

    ofstream outfile("../../k_nearest_neighbors.dta");
    // will store N x K matrix of nearest neighbors for all users
    vector<vector<float> > k_nearest(NUSERS, vector<float>(K));
    
    string line;
    vector<float> correlations(Q);
    // tanger.dta contains N x Q correlation matrix
    ifstream corr_file("../../data/tanger.dta");
    // go line by line through correlation matrix
    while ( getline(corr_file, line) ) {
        cout << "asd" << endl;
        int current_user = 0;
        // read in line as string and then split on tabs
        istringstream iss(line);
        string word;
        int word_index = 0;
        while ( getline(iss, word, '\t') ) {
            // store our correlations
            correlations[word_index] = stof(word);
            word_index++;
        }
        // sort to find nearest neighbors
        sort(correlations.begin(), correlations.end(), greater<int>());

        // MAP CORRELATIONS TO USER
        for (int i = 0; i < K; i++) {
            k_nearest[current_user][i] = correlations[i];
        }

        stringstream ss;
        for (size_t i = 0; i < correlations.size(); i++) {
            if (i != 0) {
                ss << '\t';
            }
            ss << correlations[i];
        }
        cout << ss.str() << endl;
        outfile << ss.str() << endl;

        current_user++;

    }

    outfile.close();


}