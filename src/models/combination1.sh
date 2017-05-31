#!/bin/bash
FOLDER=../kbai3
F1=$FOLDER/testtest400_1.dta
F2=$FOLDER/testtest500_3.dta
F3=../kbai_gaussian/testtestgaussian.dta
F4=../RBM/test.txt
F5=$FOLDER/test400_1.dta
F6=$FOLDER/test500_3.dta
F7=../kbai_gaussian/testgaussian.dta
F8=../RBM/quiz.txt 
python compute_rmse.py $F1 $F2 $F3 $F4 $F5 $F6 $F7 $F8 4
