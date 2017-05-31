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



// vector<vector < int > > movies_users_rating (17770);
// vector<vector < int > > movies_users_id (17770);
// vector<vector < double > > movies_users_time (17770);
// vector<double> numer (458293, 0.0);
// vector<double> denom (458293, 0.0);
vector <double> min_time (458293, 2500);








int main(int argc, char **argv) {
    ifstream in_user_means("../../data/user_means.dta");
    vector<float> user_means(NUSERS);
    int current_user = 0;
    float mean;
    while (!in_user_means.eof()) {
        in_user_means >> mean;
        user_means[current_user] = mean;
        current_user++;
    }
    in_user_means.close();

    ifstream in_movie_means("../../data/movie_means.dta");
    vector<float> movie_means(NMOVIES);
    int current_movie = 0;
    while (!in_movie_means.eof()) {
        in_movie_means >> mean;
        movie_means[current_movie] = mean;
        current_movie++;
    }
    in_movie_means.close();

    ifstream in_num_movie_ratings("../../data/num_movie_ratings.dta");
    vector<int> num_movie_ratings(NMOVIES);
    current_movie = 0;
    int num_ratings;
    while (!in_num_movie_ratings.eof()) {
        in_num_movie_ratings >> num_ratings;
        num_movie_ratings[current_movie] = num_ratings;
        current_movie++;
    }
    in_num_movie_ratings.close();

    ifstream in_num_user_ratings("../../data/num_user_ratings.dta");
    vector<int> num_user_ratings(NUSERS);
    current_user = 0;
    while (!in_num_user_ratings.eof()) {
        in_num_user_ratings >> num_ratings;
        num_user_ratings[current_user] = num_ratings;
        current_user++;
    }
    in_num_user_ratings.close();


    // calculating movie effect
    ifstream infile1("../../data/mu/2bethere.dta");
    int user, movie, time, random_shit, rating;
    vector<float> movie_effect(NMOVIES, 0.0);
    current_movie = 1;
    bool print = false;
    while (!infile1.eof()) {
        infile1 >> user >> movie >> time >> random_shit >> rating;
        if (current_movie == movie) {
            movie_effect[current_movie - 1] += (rating - GLOBALMEAN);
        }
        else if (current_movie != movie) {
            current_movie++;
            movie_effect[current_movie - 1] += (rating - GLOBALMEAN);
            print = true;
        }
        if (current_movie % 1000 == 0 && print) {
            cout << "movie effect " << current_movie << endl;
            print = false;
        }
    }
    infile1.close();

    // shrinking average movie effect
    float alpha = 22.0;
    for (int i = 0; i < NMOVIES; i++) {
        movie_effect[i] = movie_effect[i] / float(num_movie_ratings[i] + alpha);
    }

    cout << "movie effect done" << endl;


    vector<float> reg_coef_movie(NUSERS, 0.0);
    current_user = 1;
    ifstream infile2("../../data/2bethere.dta");
    float numerator = 0.0;
    float denominator = 0.0;
    while (!infile2.eof()) {
        infile2 >> user >> movie >> time >> random_shit >> rating;
        if (current_user == user) {
            numerator += (rating - GLOBALMEAN) * movie_effect[movie - 1];
            denominator += movie_effect[movie - 1] * movie_effect[movie - 1];
        }
        else if (current_user != user) {
            reg_coef_movie[current_user - 1] = numerator / denominator;
            numerator = 0;
            denominator = 0;
            current_user++;
            numerator += (rating - GLOBALMEAN) * movie_effect[movie - 1];
            denominator += movie_effect[movie - 1] * movie_effect[movie - 1];
            print = true;
        }
        if (current_user % 30000 == 0 && print) {
            cout << "movie reg " << current_user << endl;
            print = false;
        }
        if(min_time[user] > time)
        {
            min_time[user] = time;
        }

    }
    reg_coef_movie[NUSERS - 1] = numerator / denominator;
    infile2.close();
    cout << "movie regression done" << endl;

    for (int i = 0; i < NUSERS; i++) {
        reg_coef_movie[i] = reg_coef_movie[i] * float(num_user_ratings[i]) / float(num_user_ratings[i] + alpha);
    }


    vector<float> user_effect(NUSERS, 0.0);
    ifstream infile4("../../data/2bethere.dta");
    current_user = 1;
    float estimate;
    while (!infile4.eof()) {
        infile4 >> user >> movie >> time >> random_shit >> rating;
        estimate = GLOBALMEAN + reg_coef_movie[user - 1] * movie_effect[movie - 1];
        if (current_user == user) {
            user_effect[current_user - 1] += (rating - estimate);
        }
        else if (current_user != user) {
            current_user++;
            user_effect[current_user - 1] += (rating - estimate);
            print = true;
        }
        if (current_user % 30000 == 0 && print) {
            cout << "user effect " << current_user << endl;
            print = false;
        }
    }
    infile4.close();
    alpha = 7.5;
    for (int i = 0; i < NUSERS; i++) {
        user_effect[i] = user_effect[i] / float(num_user_ratings[i] + alpha);
    }

    cout << "movie effect done" << endl;


    vector<float> reg_coef_user(NUSERS, 0.0);
    current_user = 1;
    ifstream infile3("../../data/2bethere.dta");
    numerator = 0.0;
    denominator = 0.0;
    while (!infile3.eof()) {
        infile3 >> user >> movie >> time >> random_shit >> rating;
        estimate = GLOBALMEAN + reg_coef_movie[user - 1] * movie_effect[movie - 1];
        if (current_user == user) {
            numerator += (rating - estimate) * user_effect[user - 1];
            denominator += user_effect[user - 1] * user_effect[user - 1];
        }
        else if (current_user != user) {
            reg_coef_user[current_user - 1] = numerator / denominator;
            numerator = 0;
            denominator = 0;
            current_user++;
            numerator += (rating - estimate) * user_effect[user - 1];
            denominator += user_effect[user - 1] * user_effect[user - 1];
            print = true;
        }
        if (current_user % 30000 == 0 && print) {
            cout << "user reg " << current_user << endl;
            print = false;
        }
    }
    reg_coef_user[NUSERS - 1] = numerator / denominator;
    infile3.close();


    for (int i = 0; i < NUSERS; i++) {
        reg_coef_user[i] = reg_coef_user[i] * float(num_user_ratings[i]) / float(num_user_ratings[i] + alpha);
    }
    cout << "user regression done" << endl;





    // ifstream infile("../../data/qual.dta");
    // ofstream outfile("../../data/global_results.dta");
    // int current = 0;
    // while (!infile.eof()) {
    //     if (current % 100000 == 0) {
    //         cout << "generating qual " << current << endl;
    //     }
    //     infile >> user >> movie >> time;
    //     estimate = GLOBALMEAN + reg_coef_movie[user - 1] * movie_effect[movie - 1] + 
    //                 reg_coef_user[user - 1] * user_effect[user - 1];
    //     outfile << estimate << endl;
    //     current++;
    // }
    // infile.close();
    // outfile.close();
    // cout << "done" << endl;

    alpha = 160;
    vector<float> reg_coef_time(NUSERS, 0.0);
    current_user = 1;
    ifstream time_data("../../data/2bethere.dta");
    numerator = 0.0;
    denominator = 0.0;
    float xui;
    while (!time_data.eof()) {
        time_data >> user >> movie >> time >> random_shit >> rating;
        estimate = GLOBALMEAN + reg_coef_movie[user - 1] * movie_effect[movie - 1] + 
                 reg_coef_user[user - 1] * user_effect[user - 1];
        if (current_user == user) {
            xui = sqrt(time - min_time[user]) * float(num_user_ratings[user-1])/ float(num_user_ratings[user-1] + alpha);
            numerator += (rating - estimate) * xui ;
            denominator += xui * xui;
        }
        else if (current_user != user) {
            if (denominator == 0)
            {
                denominator = 1;    
            }
            reg_coef_time[current_user - 1] = numerator / denominator;
            numerator = 0;
            denominator = 0;
            current_user++;
            xui = sqrt(time - min_time[user]) * float(num_user_ratings[user-1])/ float(num_user_ratings[user-1] + alpha);
            numerator += (rating - estimate) * xui ;
            denominator += xui * xui;
            print = true;
        }
        if (current_user % 30000 == 0 && print) {
            cout << "user time reg " << current_user << endl;
            print = false;
        }
    }
    reg_coef_time[NUSERS - 1] = numerator / denominator;
    time_data.close();


    for (int i = 0; i < NUSERS; i++) {
        reg_coef_time[i] = reg_coef_time[i] * float(num_user_ratings[i]) / float(num_user_ratings[i] + alpha);
    }


    ifstream infile("../../data/qual.dta");
    ofstream outfile("../../data/global_results_with_probe.dta");
    int current = 0;
    while (!infile.eof()) {
        if (current % 100000 == 0) {
            cout << "generating qual " << current << endl;
        }
        infile >> user >> movie >> time;
        xui = sqrt(time - min_time[user]) * float(num_user_ratings[user-1])/ float(num_user_ratings[user-1] + alpha);
        estimate = GLOBALMEAN + reg_coef_movie[user - 1] * movie_effect[movie - 1] + 
                    reg_coef_user[user - 1] * user_effect[user - 1];
        estimate += reg_coef_time[user-1] * xui;
        outfile << estimate << endl;
        current++;
    }
    infile.close();
    outfile.close();
    cout << "done" << endl;





     return 0;













}
