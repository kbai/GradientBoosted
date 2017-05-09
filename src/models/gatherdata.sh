#!/sh  
#python gather_results.py 
F1=../../output_svdpp\(f\=50\)_train.dta 
F2=../../output_neighborhood_train.dta  
F3=../FREQRBM/testrbmfreqcdf.dta  
F4=../FREQRBM/testrbmfreq200.dta 
F5=../RBM/test.txt  
F6=../svd++3/testtestsvd++3log676.dta  
F7=../svd++3/testtestsvd++3400.dta 
F8=../svd++3/testtestsvd++3withhtu.dta  
F9=../svd++3/testtestsvd++3biglrsmallregu.dta 
F10=../svd++3/testtestsvd++3withfreq.dta 
F11=../svd++2/testtest200_2.dta  
F12=../../output_timefrequencybaseline_train.dta 
F13=../../output_timebaseline_train.dta  
cp $F1 ./
cp $F2 ./
cp $F3 ./
cp $F4 ./
cp $F5 ./
cp $F6 ./
cp $F7 ./
cp $F8 ./
cp $F9 ./
cp $F10 ./
cp $F11 ./
cp $F12 ./
cp $F13 ./
#python gather_results.py $F1 $F2 $F3 $F4 $F5 $F6 $F7 $F8 $F9 $F10 $F11 $F12 $F13 13
