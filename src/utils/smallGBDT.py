import numpy as np
from numpy.linalg import inv,norm
from numpy import genfromtxt 
from sklearn.ensemble import GradientBoostingRegressor
import numpy.random
data=np.loadtxt("./gradientboostinput.txt");
N = data.shape[1]
print N
print data.shape
data = data[0:N,:]
Y = data[0:N,-1]
A = numpy.random.random((N,))
print A>0.1
x_train = data[(A>0.1),0:-1]
y_train = data[(A>0.1),-1]
x_test = data[(A<0.1),0:-1]
y_test = data[(A<0.1),-1]
print y_test.size
print y_train.size
print x_train[0,:]
print y_train[0]
print norm(x_train[:,2] - y_train)/np.sqrt(y_train.size)
model1 = GradientBoostingRegressor(n_estimators = 150, max_depth = 20, learning_rate = 0.18, subsample = 1.0,verbose=2)
model2 = GradientBoostingRegressor(n_estimators = 150, max_depth = 20, learning_rate = 0.20, subsample = 1.0,verbose=2)
model1.fit(x_train,y_train)
#print model1.score(x_test, y_test)
y_predict = model1.predict(x_test)
print norm(y_predict-y_test)/np.sqrt(y_test.size)




