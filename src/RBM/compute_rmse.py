import numpy as np
from numpy.linalg import inv,norm
dat= np.loadtxt('./test500.txt')
print dat.shape
dat2 = np.loadtxt('../../data/4.dta')
print dat2.shape
b = dat2[:,4]
error = 0
n = dat.shape[0]
for i in range(n):
    error = error + (dat2[i,4] - (1*dat[i,0]+2*dat[i,1]+3*dat[i,2]+4*dat[i,3]+5*dat[i,4]))**2
print np.sqrt(error/n)
aTa = np.dot(dat.T,dat)
aTb = np.dot(dat.T,b)
ainv = inv(aTa);
c = np.dot(ainv,aTb);
print c
estimator = np.dot(c.T,dat.T)
error = estimator - b
print error.shape
print norm(error)/np.sqrt(n)



