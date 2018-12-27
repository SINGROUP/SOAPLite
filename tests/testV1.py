import soaplite
import datetime
import ase
from ase import io
import numpy as np
from soaplite import genBasis
from numpy import genfromtxt

print(soaplite)

Hpos = [[1.0,2.0,3.0]]
myAlphas, myBetas = genBasis.getBasisFunc(10.0, 5) # input:(rCut, NradBas)

atoms = ase.Atoms('H',positions=[[0.1,1,1]])

x = soaplite.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCut=10.0, nMax=5, Lmax=6,crossOver=True, eta=1.5)
print(x)
x = soaplite.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCut=10.0, nMax=5, Lmax=6,crossOver=False, eta=1.5)

atoms = ase.io.read("Structs/au40cu40.xyz")
x = soaplite.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCut=10.0, nMax=5, Lmax=6,crossOver=True, eta=1.5)
print(x)
x = soaplite.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCut=10.0, nMax=5, Lmax=6,crossOver=False, eta=1.5)
print(x)

print(x)
