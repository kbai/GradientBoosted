all: nn fastnn tdnn tfnn
nn:neuronnetwork.cc
	icc -o nn -std=c++11 -O3 ./neuronnetwork.cc
fastnn:fastnn.cc
	icc -o fastnn -pthread -std=c++11 -O3 ./fastnn.cc
tdnn:timedepnn.cc
	icc -o tdnn -pthread -std=c++11 -O3 ./timedepnn.cc
tfnn:timefreqnn.cc
	icc -o tfnn -pthread -std=c++11 -O3 ./timefreqnn.cc
