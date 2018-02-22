import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import *
from scipy.special import *
import sys

l,a  = list(map(float,sys.argv[1:]))

def myGamma(l,x):
  return gamma(l)*gammaincc(l,x)
def intPartSqr(l,a,x):
  return x**(2*l - 1)/2./(2*a)**2*(2*a*x**2)**(0.5 - l)*(gamma(l + 3.0/2.0) - myGamma(l + 3.0/2.0,2*a*x**2))
def func(alpha,l,x):
  return abs(intPartSqr(l,alpha,x)/intAllSqr(l,alpha) - .99)
def intAllSqr(l,a):
  return 0.5*gamma(l + 3.0/2.0)*(2*a)**(-l - 3.0/2.0)

alpha = np.linspace(0.001,10, num = 1000)
y = np.zeros(alpha.size)
for i,j in enumerate(alpha):
  y[i] = func(np.array([j]),l,a)

plt.plot(alpha,y)
plt.savefig('test.png')
