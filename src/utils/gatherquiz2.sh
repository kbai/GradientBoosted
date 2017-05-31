#!/bin/sh
F[1]=../../output_svdpp\(f\=50\)_predict\(0.89524\).dta 
F[2]=../../output_neighborhood_predict\(0.89802\).dta  
F[3]=../FREQRBM/quizrbmfreq100_withprobe.dta 
F[4]=../FREQRBM/quizrbmfreq200_withprobe.dta 
F[5]=../RBM/quizall.dta  
F[6]=../svd++3/testsvd++3log676.dta 
F[7]=../svd++3/testsvd++3400_with_probe.dta  
F[8]=../svd++3/testsvd++3withhtu.dta  
F[9]=../svd++3/testsvd++3biglrsmallregu.dta 
F[10]=../svd++3/testsvd++3withprobe_freq.dta 
F[11]=../svd++2/test200_2.dta  
F[12]=../../output_timefrequencybaseline_predict\(0.92643\).dta  
F[13]=../../output_timebaseline_predict\(0.95337\).dta  
F[14]=../../output_svdpp\(f\=200\)_predict\(0.89223\).dta
F[15]=../svd++3/testsvd++3100_withouthtu_withprobe.dta
F[16]=../svd++3/testsvd++350_withouthtu_withprobe.dta
#[F1]7=../svd++3/testsvd++3200_withouthtu_withprobe.dta
F[18]=../svd++3/testsvd++350withhtu_probe.dta
F[19]=../svd++3/testsvd++3200_withhtu_probe.dta
F[20]=../FREQRBM/quizrbmfreq50_withprobe.dta 
F[21]=../kbai2/testsvd++3200if0.1.dta
F[22]=../knn/knnpredict_withprobe.txt
F[23]=../../output_neighborhoodsvdpp_predict\(0.89155\).dta
F[24]=../../output_svd\(f\=500\)_predict\(0.89516\).dta
F[25]=../../output_svdpp\(f\=500\)_predict\(0.89116\).dta
F[26]=../../output_svdpp\(f\=100\)_predict\(0.89332\).dta
#F[27]=../knn/knnpredict_svd50_30_neib_withprobe.txt
F[28]=../knn/knnpredict_svd200_neib30_withprobe.txt
F[30]=../CRBM/quiz_withprobe.dta
F[31]=../../output_timefrequencysvdpp_predict\(0.88357\).dta
F[29]=../../data/weekday5.dta
F[33]=../../output_resneighborsvdpp50_predict\(0.88927\).dta
F[34]=../../output_svdpp\(f\=1000\)_predict\(0.89051\).dta
F[35]=../../output_svdpp\(f\=20\)_predict\(0.90091\).dta
F[36]=../songmingdu200_withprobe/qual200_100.dta
F[37]=../songming500_probe/qual500_110.dta
F[38]=../../output_baselinetilde_predict\(0.9797\).dta
F[39]=../../output_svd\(f\=20\)_predict\(0.90446\).dta
F[40]=../../output_neighborhoodsvdpp100_predict\(0.89037\).dta
F[41]=../../output_resneighborsvdpp200_predict\(0.88931\).dta
F[42]=../songming1000_probe/qual1000_90.dta
#F[43]=../../data/featuresvd_probe.txt
F[44]=../../data/5addinfo.txt
#F[45]=../CRBM/quiz200_withprobe.dta
F[46]=../RBM_2/quiz50_withprobe.dta

F[47]=../RBM_2/quiz100_withprobe.dta
F[48]=../FREQRBM_float/quiz50_withprobe.dta

F[49]=../FREQRBM100/quiz100probe90.217.dta
F[50]=../../output_timefreqneighbor_predict\(0.88638\).dta
F[51]=../../output_svdpp\(f\=2000\)_predict\(0.89014\).dta
F[52]=../../output_neighborhoodsvdpp100_predict\(0.89037\).dta
F[53]=../nnbaseline/probe81quiz.txt
F[54]=../knn/knnpredictrbm.dta
F[55]=../knn/knnrbmpredict.dta
F[56]=../knn/knnpredicttfrbm.txt
F[57]=../knn/ftrbmknnpredict.dta
F[58]=../../output_timeneighborhood2_predict\(0.90635\).dta
F[59]=../../output_neighborhoodsvdpp200_predict\(0.88929\).dta
F[60]=../knn/knnpredicttfrbm100.txt
F[61]=../knn/tfrbmpredictknn100.txt
F[62]=../../output_resneighbortfrbm_predict\(0.88865\).dta
#F[63]=../songming2000_probe/qual1000_40.dta
F[64]=../../output_neighborhoodsvdpp200_predict\(0.88929\).dta
F[65]=../../output_timeneighborhood2_predict\(0.90635\).dta
F[66]=../../output_timeneighbor_predict\(0.88816\).dta
F[68]=../knn/svd50plusknn_predict.dta
F[67]=../FREQRBM50_11/quiz50.dta
F[69]=../knn/svd200knn_withprobe.dta
#



#F[48]=../../output_timezerobaseline_predict\(0.92908\).dta
#F[20]=${F[48]}
#F[3]=${F[48]}
#F[4]=${F[48]}


#


#F[25]=${F[36]}
#F[26]=${F[36]}
F[6]=${F[10]}
F[7]=${F[10]}
F[8]=${F[10]}
F[9]=${F[10]}

#F[27]=../kbai/test400_all.dta
#F[28]=../../data/weekday5.dta
O=../../data/5nolabel.dta
wc -l ${F[*]}
ls -lt ${F[*]}
ls -lt ${F[*]}| wc -l


#python gather_quiz.py $F1 $F2 $F3 $F4 $F5 $F10 $F7 $F7 $F7 $F10 $F11 $F12 $F13 $F14 $F15 $F16  $F18  $F20 $F21 19
cp ${F[*]}  ./
#paste --delimiters='\t' ${F[*]} $O|tr "\t" ","|sed 's/,,/,/g'> quizinput.txt
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


