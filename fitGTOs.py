from soaplite import genBasis
import matplotlib.pyplot as plt
from soaplite import genBasis
import numpy as np
from numpy import linalg

trapSize = 1300;
nSize =40
rCut = 8.0
#myAlphas, myBetas = genBasis.getBasisFunc(rCut, nSize);
#myAlphas = np.zeros(nSize)
#myAlphas[0:10] = np.linspace(0.0001,0.005,10)
#myAlphas[0:4] = np.linspace(0.0051,0.05,4)
#myAlphas[4:8] = np.linspace(0.051,0.5,4)
#myAlphas[8:12] = np.linspace(0.51,5,4)
#myAlphas[20:30] = np.linspace(1.1,5,10)
#myAlphas[12:16] = np.linspace(1.1,5,4)
#myAlphas[12:16] = np.linspace(1.1,3,4)
myAlphas = np.zeros(nSize)
#intSize = int(myAlphas.shape[0]/10);
for i in range(0,1000000):
  #  myAlphas[:10] =0.01*np.random.rand(int(nSize/3))
 # myAlphas[10:20] =0.1*np.random.rand(int(nSize/3))
#  myAlphas[20:nSize] =2*np.random.rand(int(nSize/3))
  myAlphas = 3*np.random.rand(int(nSize))

  a = np.zeros([nSize,trapSize]);
  r = np.linspace(0,13.0,trapSize);
  f = np.cos(4*r)*np.exp(-0.5*(r-3)*(r-3));
  f[800:] = 0;
  b = np.zeros(nSize);
  S = np.zeros([nSize,nSize]);
  
  for i in range(0,nSize):
    a[i,:] = np.exp(-myAlphas[i]*r*r);
    b[i] = np.trapz(f*a[i,:], x=r);
  
  for i in range(0,nSize):
    for j in range(0,nSize):
      S[i,j] = np.trapz(a[i,:]*a[j,:],x=r)
  
  #print("S",np.exp(-myAlphas[0]*r*r))
  print("a", myAlphas)
  try:
    mulMe =linalg.inv(S)
    c = np.matmul(linalg.inv(S),b)
  except:
    print("something")
    continue
  for i in range(0,nSize):
    a[i,:] = a[i,:]*c[i]
  finalF = np.sum(a,axis=0) 
  print(sum(np.abs(finalF - f))/1500)
  if(sum(np.abs(finalF - f)/1500)< 0.01):
    plt.plot(r,finalF,r,f)
    plt.show()
print(finalF.shape)
print(c)
print(myAlphas)
