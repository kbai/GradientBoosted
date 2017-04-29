#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#define NMOVIE 17770

using namespace std;


int main(int argc, char **argv) {

    Q = argv[0];    // desired "top" users
    K = argv[1];    // desired nearest neighbors

    // will store N x K matrix of nearest neighbors for all users
    vector<vector<float> > k_nearest(NMOVIE, vector<float>(K));

    string line;
    vector<float> correlations(Q);
    // tanger.dta contains N x Q correlation matrix
    ifstream infile("./tanger.dta");}
    // go line by line through correlation matrix
    while( getline(infile, line) ) {

        int current_user = 0;
        // read in line as string and then split on tabs
        istringstream iss(line);
        string word;
        int word_index = 0;
        while( getline(iss, word, '\t') ) {
            // store our correlations
            correlations[word_index] = stof(word);
            word_index++;
        }
        sort(correlations.begin(), correlations.end(), greater<>());

        for (int i = 0; i < K; i++) {
            k_nearest[current_user][i] = correlations[i]
        }

        
    }



} 