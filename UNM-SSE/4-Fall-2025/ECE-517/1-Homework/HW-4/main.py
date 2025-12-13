import numpy as np                # General math
import numpy.matlib as matlib     # Matrix manipulation

import matplotlib.pyplot as plt   # Plotting stuff
from sklearn import svm           # The SVM
from sklearn.linear_model import Ridge

def data(N,sigma):
  w=np.ones(10)/np.sqrt(10)
  w1=[1., 1., 1., 1., 1., -1., -1., -1., -1., -1.]/np.sqrt(10)
  w2=[-1., -1., 0, 1., 1., -1., -1., 0, -1., -1.]/np.sqrt(8)
  x=np.zeros((4,10))
  x[1,:]=x[0,:]+sigma*w1
  x[2,:]=x[0,:]+sigma*w2
  x[3,:]=x[2,:]+sigma*w1
  X1=x+sigma*matlib.repmat(w,4,1)/2
  X2=x-sigma*matlib.repmat(w,4,1)/2
  X1=matlib.repmat(X1,2*N,1)
  X2=matlib.repmat(X2,2*N,1)
  X=np.concatenate((X1, X2), axis=0)
  Y=np.concatenate((np.ones(4*2*N), -np.ones(4*2*N)),axis=0)
  Z=np.random.permutation(16*N)
  Z=Z[:N]
  X=X[Z,:]
  X=X+0.2*sigma*np.random.randn(N,10)
  Y=Y[Z]
  return X,Y

def data2(Ntr,Ntst,sigma):
  Xtr,ytr=data(Ntr,sigma)
  Xtst,ytst=data(Ntr,sigma)
  return Xtr,ytr,Xtst,ytst

C = 0.01
sigma = 1   # Standard deviation of the clusters. Try 0.5,1,3

Xtr,ytr,Xtst,ytst=data2(100,1000,sigma)

clf = svm.SVC(kernel="linear", C=C) # Declare a SVM object
clf.fit(Xtr, ytr)                   # We call the training method named .fit to traini the SVM
ytr_= clf.predict(Xtr)              # We call the test method to predict labels ytr
ytst_= clf.predict(Xtst)            # We call the test method to predict labels ytst

Etr = np.mean(abs(ytr-ytr_)/2)      # Empirical risk (error rate or the estimation of the error probability)
Etst = np.mean(abs(ytst-ytst_)/2)   # Estimation of the actual risk

print("Train error: ", Etr)
print("Test error: ", Etst)

print(clf.dual_coef_)
