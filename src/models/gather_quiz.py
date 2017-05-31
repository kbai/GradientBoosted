import numpy as np
from numpy.linalg import inv,norm
import sys
#python compute_rmse.py testtest200_2.dta testtest200_1.0.dta  ../RBM//rbmtest100.dta  ./test200_2.dta ./test200_1.0.dta ../RBM/quiz100.dta 3
print sys.argv
datam = np.array([])
quiz = np.array([])
n = int(sys.argv[-1])
for a in sys.argv[1:n+1]:
    print "load test set"
    print a
    dat = np.loadtxt(a)
    print dat.shape
    print datam.shape
    datam = np.vstack((dat.T,datam)) if datam.size else dat.T
#for q in sys.argv[n+1:2*n+1]:
#    print "load quiz set"
#    print q
#    dat = np.loadtxt(q)
#    quiz = np.vstack((dat.T,quiz)) if quiz.size else dat.T


b = np.loadtxt('../../data/5.dta')
b = b[:,0:-1]
print b.shape
print datam.shape
datam = np.vstack((datam,b.T))
np.savetxt( './quizinput.txt',datam.T,delimiter='\t')
#print norm(error)/np.sqrt(n)
#
#
