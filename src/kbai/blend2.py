import numpy as np
import sys
print sys.argv[1]
print sys.argv[2]
print sys.argv[3]
print sys.argv[4]
a = float(sys.argv[3])
b = float(sys.argv[4])
dat1 = np.loadtxt(sys.argv[1])
dat2 = np.loadtxt(sys.argv[2])
dat3 = dat1 * a + dat2 * b
np.savetxt('./blendmodel.txt',dat3,delimiter='\n')
