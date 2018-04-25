import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import soapPy
import ase
import numpy as np
from numpy import genfromtxt

# example structure
atoms = ase.io.read("Cu_110.pdb")
atoms_c = atoms.copy()

print("local soap based on given positions")
Hpos = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0],[10.0,11.0,12.0]]
#Hpos = genfromtxt('H.dat').tolist()
x = soapPy.get_periodic_soap_locals(atoms_c, Hpos, rCutHard=8.0, NradBas=5, Lmax=5) #rCutSoft = rCutHard - 3.0
print("soap for each atom in structure")
y = soapPy.get_periodic_soap_structure(atoms, rCutHard=8.0, NradBas=5, Lmax=5)




# SOAP solution: x
#print("soap size: ", np.shape(x))
np.savetxt('periodic_test.txt',x)
np.savetxt('periodic_structure_test.txt',y)


