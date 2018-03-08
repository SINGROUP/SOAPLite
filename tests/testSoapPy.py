import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import soapPy
import ase
import numpy as np
from numpy import genfromtxt

# example structure
atoms = ase.io.read("au40cu40.xyz")
#atoms = ase.Atoms('H2O', positions = [[0.0, 0.1, 0.2 ], [3.2,4.5,5.5],[5.0, 6.0, 7.0]])

print("local soap based on given positions")
Hpos = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0],[10.0,11.0,12.0]]
#Hpos = genfromtxt('H.dat').tolist()

x = soapPy.get_soap_locals(atoms, Hpos, rCutHard=8.0, NradBas=5, Lmax=5) #rCutSoft = rCutHard - 3.0
print("soap for each atom in structure")
y = soapPy.get_soap_structure(atoms, rCutHard=8.0, NradBas=5, Lmax=5) #rCutSoft = rCutHard - 3.0

# SOAP solution: x
#print("soap size: ", np.shape(x))
np.savetxt('test.txt',x)
np.savetxt('structure_test.txt',y)

