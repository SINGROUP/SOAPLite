import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import *

def myFunctions(x,xt2,nMax):
  alphas = np.loadtxt("alphasPyl1.dat")
  y = np.zeros([nMax,len(x)]);
  yt2 = np.zeros([nMax,len(x)]);
  y[0,:] = x*np.exp(-alphas[0]*x*x);
  y[1,:] = x*np.exp(-alphas[1]*x*x);
  y[2,:] = x*np.exp(-alphas[2]*x*x);
  y[3,:] = x*np.exp(-alphas[3]*x*x);
  y[4,:] = x*np.exp(-alphas[4]*x*x);
  y[5,:] = x*np.exp(-alphas[5]*x*x);
  yt2[0,:] = xt2*np.exp(-alphas[0]*xt2*xt2);
  yt2[1,:] = xt2*np.exp(-alphas[1]*xt2*xt2);
  yt2[2,:] = xt2*np.exp(-alphas[2]*xt2*xt2);
  yt2[3,:] = xt2*np.exp(-alphas[3]*xt2*xt2);
  yt2[4,:] = xt2*np.exp(-alphas[4]*xt2*xt2);
  yt2[5,:] = xt2*np.exp(-alphas[5]*xt2*xt2);
  return y,yt2

def getGns(y,yt2,wi,x,xt2,rCut,nMax):
  mat = np.zeros([nMax,nMax]);
  gn = np.zeros([nMax,len(x)]);
  for i in range(0,nMax):
    for j in range(0,nMax):
      mat[i,j] = 2*rCut*0.5*np.sum(wi*xt2*xt2*yt2[i,:]*yt2[j,:]);

  print("M:",mat)
  invMat = sqrtm(inv(mat));
  print("inv:", invMat)
  for n in range(0,nMax):
    for a in range(0,nMax):
      gn[n,:] = gn[n,:] + invMat[n,a]*y[a,:] 

#  for i in range(0,nMax): # Check orthnorm 
#    for j in range(0,nMax):
#      print(rCut*0.5*np.sum(wi*x*x*gn[i,:]*gn[j,:]))

  return gn

if __name__=='__main__':

  lMax  = 7;
  nMax  = 6;
  rCut  = 15.0;
  rSize = 100; #Do not touch

  xi = np.loadtxt("GL.x");
  wi = np.loadtxt("GL.w");

  xt2 = 0.5*2*rCut*(xi + 1);
  x = 0.5*rCut*(xi + 1);
  f,f2 = myFunctions(x,xt2,nMax);
  gns = getGns(f,f2,wi,x,xt2,rCut,nMax);
  np.savetxt("gss.dat",gns)
  np.savetxt("rw.dat",x)
  for i in range(0,nMax):
    plt.plot(x,gns[i,:])
#    plt.plot(x,f[i,:])
    print(gns)
  plt.show()
