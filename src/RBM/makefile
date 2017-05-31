CC= mpicxx -std=c++11
executable: rbm_mpi
rbm.o:
	$(CC) -DMPIU  -c  RBM5_MPI.cpp
rbm_mpi: rbm.o
	$(CC) -DMPIU -o rbm_mpi RBM5_MPI.o



