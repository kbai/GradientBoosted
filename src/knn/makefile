CC=icc -std=c++11 -O3 
all: knn cc knn_fast knn_fast_raw
knn:compute_knn_movie_test.cpp
	$(CC) -o knn ./compute_knn_movie_test.cpp	
cc:compute_correlation_between_movies.cpp
	$(CC) -o cc ./compute_correlation_between_movies.cpp
knn_fast:compute_knn_movie_hash.cpp
	$(CC) -o knn_fast ./compute_knn_movie_hash.cpp	
knn_fast_raw:compute_knn_movie_hash_raw.cpp
	$(CC) -o knn_fast_raw ./compute_knn_movie_hash_raw.cpp	

