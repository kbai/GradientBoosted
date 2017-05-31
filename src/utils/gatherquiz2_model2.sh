#!/bin/sh
F[1]=../../output_svdpp\(f\=50\)_predict\(0.89524\).dta 
F[2]=../../output_neighborhood_predict\(0.89802\).dta  
F[3]=../FREQRBM/quizrbmfreq100_withprobe.dta 
F[4]=../FREQRBM/quizrbmfreq200_withprobe.dta 
F[5]=../RBM/quizall.dta  
F[12]=../../output_timefrequencybaseline_predict\(0.92643\).dta  
F[13]=../../output_timebaseline_predict\(0.95337\).dta  
F[14]=../../output_svdpp\(f\=200\)_predict\(0.89223\).dta
F[20]=../FREQRBM/quizrbmfreq50_withprobe.dta 
F[22]=../knn/knnpredict_withprobe.txt
F[23]=../../output_neighborhoodsvdpp_predict\(0.89155\).dta
F[24]=../../output_svd\(f\=500\)_predict\(0.89516\).dta
F[25]=../../output_svdpp\(f\=500\)_predict\(0.89116\).dta
F[26]=../../output_svdpp\(f\=100\)_predict\(0.89332\).dta
F[27]=../knn/knnpredict_svd50_30_neib_withprobe.txt
F[28]=../knn/knnpredict_svd200_neib30_withprobe.txt
F[30]=../CRBM/quiz_withprobe.dta
F[31]=../../output_timefrequencysvdpp_predict\(0.88357\).dta
F[33]=../../output_resneighborsvdpp50_predict\(0.88927\).dta
F[34]=../../output_svdpp\(f\=1000\)_predict\(0.89051\).dta
F[35]=../../output_svdpp\(f\=20\)_predict\(0.90091\).dta

#F[27]=../kbai/test400_all.dta
#F[28]=../../data/weekday5.dta
O=../../data/5nolabel.dta
wc -l ${F[*]}


#python gather_quiz.py $F1 $F2 $F3 $F4 $F5 $F10 $F7 $F7 $F7 $F10 $F11 $F12 $F13 $F14 $F15 $F16  $F18  $F20 $F21 19
paste --delimiters='\t' ${F[*]} $O|tr "\t" ","|sed 's/,,/,/g'> quizinput.txt
# ../../output_svdpp\(f\=50\)_predict\(0.89524\).dta	0      772
# ../../output_neighborhood_predict\(0.89802\).dta		1      595
# ../FREQRBM/quizrbmfreq100_withprobe.dta				2-6   1507
# ../FREQRBM/quizrbmfreq200_withprobe.dta				7-11  2389
# ../RBM/quizall.dta									12-16      2606
# ../svd++3/testsvd++3log676.dta						17         609
# ../svd++3/testsvd++3400_with_probe.dta				18         618
# ../svd++3/testsvd++3withhtu.dta						19         513
# ../svd++3/testsvd++3biglrsmallregu.dta				20         482
# ../svd++3/testsvd++3withprobe_freq.dta				21         497
# ../svd++2/test200_2.dta								22         560
# ../../output_timefrequencybaseline_predict\(0.92643\).dta	23     991
# ../../output_timebaseline_predict\(0.95337\).dta  13		24     687


