#!/bin/sh

#F[25]=${F[36]}
#F[26]=${F[36]}

#F[27]=../kbai/test400_all.dta
#F[28]=../../data/weekday5.dta
O=../../data/1map.dta


#python gather_quiz.py $F1 $F2 $F3 $F4 $F5 $F10 $F7 $F7 $F7 $F10 $F11 $F12 $F13 $F14 $F15 $F16  $F18  $F20 $F21 19
paste $O| awk -v OFS="," '{print $6, $1, $2, $3, $4, $5}'> training_gbdt.txt
paste ../../data/4map.dta| awk -v OFS=',' '{print $6 ,$1 ,$2 ,$3 ,$4 ,$5}'> test_gbdt.txt


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


