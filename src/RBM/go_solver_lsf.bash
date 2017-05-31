#!/bin/bash

## job name and output file
#BSUB -J go_solver
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -e OUTPUT_FILES/%J.e

###########################################################
# USER PARAMETERS

## 4 CPUs ( 4  ), walltime 1 hour
#BSUB -n 4
#BSUB -R span[ptile=4]
#BSUB -q short_parallel

#BSUB -a openmpi

###########################################################
echo $pd    
#pd is defined as the current directory in GoSimulation 
cd $pd
#cd $currentpath
# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed

# total number of nodes is the product of the values read
numnodes=50


echo starting run in current directory $PWD
echo " "

sleep 2
#mpirun nvprof -o main%p.nvprof ./xspecfem3D
#mpirun  ./xdecompose_mesh_mpi 125 ../MESH ../DATABASES_MPI 5 5 5
mpirun ./rbm_mpi

echo "finished successfully"

