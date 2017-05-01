import numpy as np
from numpy.linalg import inv,norm
import sys
print sys.argv
datam = np.array([])
quiz = np.array([])
n = int(sys.argv[-1])
for a in sys.argv[1:n+1]:
    print "load test set"
    print a
    dat = np.loadtxt(a)
    datam = np.vstack((dat,datam)) if datam.size else dat
for q in sys.argv[n+1:2*n+1]:
    print "load quiz set"
    print q
    dat = np.loadtxt(q)
    quiz = np.vstack((dat,quiz)) if quiz.size else dat

b = np.loadtxt('../../data/4.dta')
b = b[:,4]
print b.shape
print b
print datam.shape
datam = datam.T
aTa = np.dot(datam.T,datam)
print aTa.shape
aTb = np.dot(datam.T,b)
print aTb.shape
ainv = inv(aTa)
print ainv.shape
c=np.dot(ainv,aTb)
print c.shape
estimator = np.dot(c.T,datam.T)
print estimator.shape
print b.shape
error = estimator - b
print c
print norm(error)/np.sqrt(error.size)
equize = np.dot(c.T,quiz)
np.savetxt('./testblend.txt',equize,delimiter='\n')
#dat= np.loadtxt('./test.txt')
#print dat.shape
#dat2 = np.loadtxt('../../data/4.dta')
#print dat2.shape
#b = dat2[:,4]
#error = 0
#n = dat.shape[0]
#for i in range(n):
#    error = error + (dat2[i,4] - (1*dat[i,0]+2*dat[i,1]+3*dat[i,2]+4*dat[i,3]+5*dat[i,4]))**2
#print np.sqrt(error/n)
#aTa = np.dot(dat.T,dat)
#aTb = np.dot(dat.T,b)
#ainv = inv(aTa);
#c = np.dot(ainv,aTb);
#print c
#estimator = np.dot(c.T,dat.T)
#error = estimator - b
#print error.shape
#print norm(error)/np.sqrt(n)



