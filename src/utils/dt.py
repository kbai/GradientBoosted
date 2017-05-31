import numpy as np
from numpy import genfromtxt 
from sklearn.ensemble import GradientBoostingRegressor
data=np.loadtxt("./gradientboostinput.txt");
print data.shape
x_train = data[:,[6,8,9,10,11,12]]
y_train = data[:,-1]
print x_train[0,:]
print y_train[0]

model1 = GradientBoostingRegressor(n_estimators = 1, max_depth = 20, learning_rate = 0.18, subsample = 0.9, verbose=2)
model2 = GradientBoostingRegressor(n_estimators = 1, max_depth = 20, learning_rate = 0.20, subsample = 1.0, verbose=2)
model1.fit(x_train, y_train)
#y_predict = model1.predict(x_predict)




