#!/sh  
#python gather_results.py 
F[1]=../../output_svdpp\(f\=50\)_train\(0.89524\).dta
F[2]=../../output_neighborhood_train\(0.89802\).dta
F[3]=../FREQRBM/testrbmfreqcdf.dta  
F[4]=../FREQRBM/testrbmfreq200.dta 
F[5]=../RBM/test.txt  
F[12]=../../output_timefrequencybaseline_train\(0.92643\).dta
F[13]=../../output_timebaseline_train\(0.95337\).dta
F[14]=../../output_svdpp\(f\=200\)_train\(0.89223\).dta
F[20]=../FREQRBM/testrbmfreq50.dta 
F[22]=../knn/knnpredict.txt
F[23]=../../output_neighborhoodsvdpp_train\(0.89155\).dta
F[24]=../../output_svd\(f\=500\)_train\(0.89516\).dta
F[25]=../../output_svdpp\(f\=500\)_train\(0.89116\).dta
F[26]=../../output_svdpp\(f\=100\)_train\(0.89332\).dta
F[27]=../knn/knnpredict_svd50_30_neib.txt
F[28]=../knn/knnpredict_svd200_neib30.txt
F[30]=../CRBM/test100
F[31]=../../output_timefrequencysvdpp_train\(0.88357\).dta
F[33]=../../output_resneighborsvdpp50_train\(0.88927\).dta
F[34]=../../output_svdpp\(f\=1000\)_train\(0.89051\).dta
F[35]=../../output_svdpp\(f\=20\)_train\(0.90091\).dta
O=../../data/4.dta
echo ${F[*]}
wc -l ${F[*]} 
#python gather_results.py $F1 $F2 $F3 $F4 $F5 $F6 $F7 $F8 $F9 $F10 $F11 $F12 $F13 $F14 $F15 $F16  $F18  $F20 $F21 19
FLG=--delimiters\="\t"
paste  ${F[*]}  $O> gradientboostinput.txt
