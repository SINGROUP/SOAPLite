import soaplite
import datetime
import ase
import ase.io
import numpy as np
from soaplite import genBasis
from numpy import genfromtxt
import time

x = 2.00
for i in range(1000):
  for j in range(2, 12):
    print(i,j)
    myAlphas, myBetas = genBasis.getBasisFunc(x, j)# input:(rCut, nMax)
    np.savetxt("myAlphasR"+str(i)+"N"+str(j), myAlphas)
    np.savetxt("myBetasR"+str(i)+"N"+str(j), myBetas.flatten())
  x = x + 0.01
