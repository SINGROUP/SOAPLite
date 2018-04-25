import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import soapPy
import datetime
import ase
import numpy as np
import genBasis
from numpy import genfromtxt

# example structure
atoms = ase.io.read("reduced_mos2h_n6_cut0_s0_10000.xyz")

print("local soap based on given positions")
Hpos = genfromtxt('reduced_mos2H.dat').tolist()

for i in range(6,13):
  for j in range(5,10):
    myAlphas, myBetas = genBasis.getBasisFunc(6.0, i)
    a = datetime.datetime.now()
    x = soapPy.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCutHard=6.0, NradBas=i, Lmax=j) #rCutSoft = rCutHard - 3.0
    b = datetime.datetime.now()
    c = b - a
    np.savetxt('testMoS2Hs.txt_R' + str(i)+ '_L' + str(j) ,x)
    print("The runtime for soaps was ",c.microseconds, "micro seconds.")


