from scipy.special import *
from scipy.linalg import *
from scipy.optimize import *
import numpy as np

#--------------------------------------------------
def myGamma(l,x):
    return gamma(l)*gammaincc(l,x)
#--------------------------------------------------
def intAllMat(l,a):
    m = np.zeros((a.shape[0],a.shape[0]))
    m[:,:] = a
    m = m + m.transpose()
#    print(m)
    return 0.5*gamma(l + 3.0/2.0)*m**(-l-3.0/2.0)
#--------------------------------------------------
def intAllSqr(l,a):
    return 0.5*gamma(l + 3.0/2.0)*(2*a)**(-l-3.0/2.0)
#--------------------------------------------------
def intPartSqr(l,a,x):
    return x**(2*l - 1)/2./(2*a)**2*(2*a*x**2)**(0.5 - l)*(gamma(l + 3.0/2.0) - myGamma(l + 3.0/2.0,2*a*x**2))
#--------------------------------------------------
def minimizeMe(alpha,l,x):
    return np.abs(intPartSqr(l,alpha,x)/intAllSqr(l,alpha) - .99)
#--------------------------------------------------
def findAlpha(l,a, alphaSpace):
    alphas = np.zeros(a.shape[0])
#    alphaSpace = np.linspace(0.001, 100, 100000)
    for i,j in enumerate(a):
      initG = alphaSpace[np.argmin(minimizeMe(alphaSpace,l,j))]
#      print(initG)
      alphas[i] = fmin(minimizeMe, x0=initG, args=(l,j), disp=False)
     #initG = alphas[i]
    return alphas 
#--------------------------------------------------
def getOrthNorm(X):
    x = sqrtm(inv(X))
#    print("inv",x)
    return x
#--------------------------------------------------
def getBasisFunc(rcut, n):
    a = np.linspace(1,rcut,n)
    alphasFull = np.array([])
    betas = np.array([])
    betasFull = np.array([])
    alphaSpace = np.array([])
    val = 0.00001
    for i in range(1100):
        val = val*1.015
#        print(val)
        alphaSpace = np.append(alphaSpace, val)
    for l in range(0,10):
       alphas = findAlpha(l,a, alphaSpace)
       print(alphas)
       betas = getOrthNorm(intAllMat(l,alphas))
       alphasFull = np.append(alphasFull, alphas)
       betasFull  = np.append(betasFull , betas)
#    np.savetxt('alphasPy.dat', alphasFull)
#    np.savetxt('betasPy.dat',  betasFull)
    return  alphasFull, betasFull
#--------------------------------------------------
if __name__ == '__main__':
    alphas, betas = getBasisFunc(5.0, 10)
#    print(alphas, betas)
    np.savetxt('alphasPy.dat', alphas)
    np.savetxt('betasPy.dat',  betas)










