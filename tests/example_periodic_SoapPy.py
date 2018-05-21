import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import soaplite
import genBasis
import ase
import numpy as np

# example structure
atoms = ase.io.read("../Structs/Cu_110.pdb")
#atoms_c = atoms.copy()

myAlphas, myBetas = genBasis.getBasisFunc(10.0, 5) # input: (rCut, NradBas)
print("local soap based on given positions")
Hpos = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0],[10.0,11.0,12.0]]
x = soaplite.get_periodic_soap_locals(atoms, Hpos, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=5) 
print("soap for each atom in structure")
y = soaplite.get_periodic_soap_structure(atoms, myAlphas, myBetas, 10.0, NradBas=5, Lmax=5)




# SOAP solution: x
#print("soap size: ", np.shape(x))
np.savetxt('periodic_test.txt',x)
np.savetxt('periodic_structure_test.txt',y)


