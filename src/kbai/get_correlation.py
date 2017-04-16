import numpy as np
import os
import sys
import math
p1 = "awk '{if($2 =="+str(sys.argv[1])+" ) print }' ./1.dta > file1"
p2 = "awk '{if($2 =="+str(sys.argv[2])+" ) print }' ./1.dta > file2"

os.system(p1)
os.system(p2)
set1 = np.loadtxt('file1')
set2 = np.loadtxt('file2')
i = 0
j = 0
sumx = 0
sumy = 0
n = 0
sumxx = 0
sumyy = 0
sumxy = 0
while((i < set1.shape[0]) & (j < set2.shape[0])):
    if(set1[i,0] < set2[j,0]):
        i = i+1
        continue
    if(set1[i,0] > set2[j,0]):
        j = j+1
        continue
    if(set1[i,0] == set2[j,0]):
        print set1[i,4] , set2[j,4]
        x = set1[i,4]
        y = set2[j,4]
        j = j+1
        i = i+1
    sumx = sumx + x
    sumy = sumy + y
    n = n+1
    sumxx = sumxx + x*x
    sumyy = sumyy + y*y
    sumxy = sumxy + x*y
re = ( n*sumxy -  sumx * sumy )/math.sqrt(1e-6+sumxx*n - sumx*sumx)/math.sqrt(1e-6+sumyy*n - sumy*sumy)
print n
print re
