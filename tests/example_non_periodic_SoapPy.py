import soaplite
import datetime
import ase
import numpy as np
from soaplite import genBasis
from numpy import genfromtxt
import time

#-------------- Define Structure ------------------------------------
#atoms = ase.Atoms('CH4', positions=[[-1.5, 0.0, 0.0],[ 1.5, 0.0, 0.0],[0.0, 1.5, 0.0],[ -1.5, 0.0, 0.0]])
#atoms = ase.io.read("../Structs/mos2_51.xyz")
atoms = ase.io.read("../Structs/au40cu40.xyz")

#-------------- Define Position of Local Environment ----------------
#Hpos = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0],[10.0,11.0,12.0]]
#Hpos = genfromtxt('../Structs/mos2H.xyz').tolist()
Hpos = genfromtxt('../Structs/au40cu40H.dat').tolist()

#-------------- set Basis Function (rCut--soft, N_max) Environment ----------------
myAlphas, myBetas = genBasis.getBasisFunc(10.0, 5) # input:(rCut, NradBas)

#-------------- run local chemical environments on desired points ----------------
start = time.time()
x = soaplite.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=5,crossOver=True) 
endTime = time.time()
totalTime = endTime - start
print("Soap ran in seconds:", totalTime)
np.savetxt('au40cu40H.txt',x)

#-------------- run local chemical environments on each atom ----------------
start = time.time()
y = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=5,crossOver=True) 
endTime = time.time()
totalTime = endTime - start

print("Soap ran in seconds:", totalTime)
np.savetxt('structure_test.txt',y)
#------------ Check Results in *txt files --------------------------------
#------------ Make sure that rCut and N_max are the same across all calculations 
