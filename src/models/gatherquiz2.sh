#!/bin/sh
F1=../../output_svdpp\(f\=50\)_predict\(0.89524\).dta 
F2=../../output_neighborhood_predict\(0.89802\).dta  
F3=../FREQRBM/quizrbmfreq100_withprobe.dta 
F4=../FREQRBM/quizrbmfreq200_withprobe.dta 
F5=../RBM/quizall.dta  
F6=../svd++3/testsvd++3log676.dta 
F7=../svd++3/testsvd++3400_with_probe.dta  
F8=../svd++3/testsvd++3withhtu.dta  
F9=../svd++3/testsvd++3biglrsmallregu.dta 
F10=../svd++3/testsvd++3withprobe_freq.dta 
F11=../svd++2/test200_2.dta  
F12=../../output_timefrequencybaseline_predict\(0.92643\).dta  
F13=../../output_timebaseline_predict\(0.95337\).dta  
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

#python gather_quiz.py $F1 $F2 $F3 $F4 $F5 $F7 $F7 $F7 $F7 $F10 $F11 $F12 $F13 13
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


