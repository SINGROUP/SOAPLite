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
#atoms = ase.Atoms('CH4', positions=[[-1.5, 0.0, 0.0],[ 1.5, 0.0, 0.0],[0.0, 1.5, 0.0],[ -1.5, 0.0, 0.0]])

# example structure
#atoms = ase.io.read("mos2_51.xyz")
#atoms = ase.io.read("au40cu40.xyz")
atoms = ase.io.read("mosHFake.xyz")
#atoms = ase.Atoms('H2O', positions = [[0.0, 0.1, 0.2 ], [3.2,4.5,5.5],[5.0, 6.0, 7.0]])

#print("local soap based on given positions")
#Hpos = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0],[10.0,11.0,12.0]]
#Hpos = genfromtxt('H.dat').tolist()
#Hpos = genfromtxt('mos2H.xyz').tolist()

myAlphas, myBetas = genBasis.getBasisFunc(6.0, 10)
#a = datetime.datetime.now()
#x = soapPy.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCutHard=5.0, NradBas=5, Lmax=j) #rCutSoft = rCutHard - 3.0
#b = datetime.datetime.now()
#c = b - a
#np.savetxt('testMoS2Hs.txt_N' + str(i) + '_L' + str(j) ,x)
#print("The runtime for soaps was ",c.microseconds, "micro seconds.")

print("soap for each atom in structure")
y = soapPy.get_soap_structure(atoms, myAlphas, myBetas, rCut=6.0, NradBas=10, Lmax=9) 

#print("XXX")
# SOAP solution: x
#print("soap size: ", np.shape(x))
np.savetxt('structure_test.txt',y)
