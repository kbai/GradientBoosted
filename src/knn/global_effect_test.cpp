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
#define NMOVIES 17770
#define NUSERS 458293
#define GLOBALMEAN 3.60861
using namespace std;


int main(int argc, char **argv) {
    ifstream in_user_means("../../data/user_means.dta");
    vector<float> user_means(NUSERS);
    int current_user = 0;
    float mean;
    while (!in_user_means.eof()) {
        in_user_means >> mean;
        user_means[current_user] = mean;
    }
    in_user_means.close();

    ifstream in_movie_means("../../data/movie_means.dta");
    vector<float> movie_means(NMOVIES);
    int current_movie = 0;
    while (!in_movie_means.eof()) {
        in_movie_means >> mean;
        movie_means[current_movie] = mean;
    }
    in_movie_means.close();

    ifstream in_num_movie_ratings("../../data/num_movie_ratings.dta");
    vector<int> num_movie_ratings(NMOVIES);
    current_movie = 0;
    int num_ratings;
    while (!in_num_movie_ratings.eof()) {
        in_num_movie_ratings >> num_ratings;
        num_movie_ratings[current_movie] = num_ratings;
    }
    in_num_movie_ratings.close();

    ifstream in_num_user_ratings("../../data/num_user_ratings.dta");
    vector<int> num_user_ratings(NUSERS);
    current_user = 0;
    while (!in_num_user_ratings.eof()) {
        in_num_user_ratings >> num_ratings;
        num_user_ratings[current_user] = num_ratings;
    }
    in_num_user_ratings.close();


    // calculating movie effect
    ifstream infile("../../data/mu/1.dta");
    int user, movie, time, random_shit, rating;
    vector<float> avg_movie(NMOVIES, 0.0);
    current_movie = 0;
    while (!infile.eof()) {
        infile >> user >> movie >> time >> random_shit >> rating;
        if ((current_movie + 1) == movie + 1) {
            avg_movie[current_movie] += (rating - GLOBALMEAN);
        }
        else if ((current_movie + 1) != movie) {
            current_movie++;
            avg_movie[current_movie] += (rating - GLOBALMEAN);
        }

        
    }

    // shrinking average movie effect
    float alpha = 22.0;
    for (int i = 0; i < NMOVIES; i++) {
        avg_movie[i] = avg_movie[i] / float(num_movie_ratings[i]);
        avg_movie[i] = avg_movie[i] * num_movie_ratings[i] / (num_movie_ratings[i] + alpha);
    }
    infile.close();

    cout << "movie effect done?"

    vector<float> reg_coef(NUSERS, 0.0);
    current_user = 0;
    ifstream infile2("../../data/1.dta");
    float numerator = 0;
    float denominator = 0;
    while (!infile2.eof()) {
        infile >> movie >> user >> time >> random_shit >> rating;
        if ((current_user + 1) == user) {
            numerator += (rating - GLOBALMEAN) * avg_movie[movie - 1];
            denominator += avg_movie[movie - 1] * avg_movie[movie - 1];
        }
        else if ((current_user + 1) != user) {
            reg_coef[current_user] = numerator / denominator;
            numerator = 0;
            denominator = 0;
            current_user++;
            numerator += (rating - GLOBALMEAN) * avg_movie[movie - 1];
            denominator += avg_movie[movie - 1] * avg_movie[movie - 1];
        }
    }

    alpha = 7.5;
    for (int i = 0; i < NUSERS; i++) {
        reg_coef[i] = reg_coef[i] * num_user_ratings[i] / (num_user_ratings[i] + alpha);
        cout << reg_coef[i];
    }

}