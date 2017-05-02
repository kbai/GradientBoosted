#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#define NMOVIE 17770
#define NUSERS 458293

using namespace std;

// helper function for sorting vector of vectors by column
bool sortcol(const vector<float>& v1, const vector<float>& v2) {
     return v1[1] > v2[1];
}


int main(int argc, char **argv) {

    cout << "Finding K-nearest neighbors" << endl;
    //int Q = atoi(argv[1]);    // desired "top" users
    int Q = 10;
    int K = 5;
    //int K = atoi(argv[2]);    // desired nearest neighbors
    // reading in to find Q top users
    ifstream user_file("../../data/sorted_users.dta");
    float user;
    int ratings;
    int current_user = 0;
    vector<float> top_users(Q);
    while (!user_file.eof()) {
        if (current_user >= Q) {
            break;
        }
        user_file >> user >> ratings;
        top_users[current_user] = user;
        cout << user << endl;
        current_user++;
    }
    user_file.close();

    ofstream outfile("../../data/k_nearest_neighbors.dta");
    // will store N x K matrix of nearest neighbors for all users
    //vector<vector<float> > k_nearest(NUSERS, vector<float>(K));
    
    string line;
    // will store Q x 2 matrix where col1 is top Q users and col2 is correlations
    vector<vector<float> > correlations(Q, vector<float>(2));

    cout << "loading correlation matrix" << endl;
    // tanger.dta contains N x Q correlation matrix
    ifstream corr_file("../../data/tanger.dta");
    // go line by line through correlation matrix
    current_user = 0;
    while ( getline(corr_file, line) ) {

        if (current_user % 50000 == 0) {
            cout << "calculated " << current_user << " users" << endl;
        }
        // read in line as string and then split on tabs
        istringstream iss(line);
        string word;
        int word_index = 0;
        while ( getline(iss, word, '\t') ) {
            // store our user and correlation
            correlations[word_index][0] = top_users[word_index];
            correlations[word_index][1] = stof(word, NULL); 
           word_index++;
        }
        // sort by second column find nearest neighbors
        sort(correlations.begin(), correlations.end(), sortcol);

        if (current_user == 0) {
            for (int i = 0; i < K; i++) {
                cout << correlations[i][0] << endl;
            }
        }


        stringstream ss;
        for (int i = 0; i < K; i++) {
            if (i != 0) {
                ss << '\t';
            }
            ss << correlations[i][0];
        }
        outfile << ss.str() << endl;

        current_user++;

    }

    outfile.close();
    cout << "we done tho" << endl;
}