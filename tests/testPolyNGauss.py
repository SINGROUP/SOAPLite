import soaplite
import datetime
import ase
import numpy as np
from soaplite import genBasis
from numpy import genfromtxt

Hpos =[[0,0,0]]
atoms = ase.Atoms('H',positions=[[0.0,0,0]])
x = soaplite.get_soap_locals_gauss(atoms, Hpos,  rCut=5.0, nMax=5, Lmax=5,  all_atomtypes=[], eta=1.5)
print(x)
x = soaplite.get_soap_locals_poly(atoms, Hpos, rCut=5.0, nMax=5, Lmax=5, all_atomtypes=[], eta=1.5)
print(x)
np.savetxt('test3.xyz',x)
