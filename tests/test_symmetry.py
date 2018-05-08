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
import time
from ase.structure import molecule
from ase.collections import g2
from scipy.spatial.distance import pdist
# GLOBALS

IS_PASS = True

#-------------- Define Structure ------------------------------------
#atoms = ase.Atoms('CH4', positions=[[-1.5, 0.0, 0.0],[ 1.5, 0.0, 0.0],[0.0, 1.5, 0.0],[ -1.5, 0.0, 0.0]])
#atoms = ase.io.read("../Structs/mos2_51.xyz")
atoms = ase.io.read("../Structs/au40cu40.xyz")

#-------------- Define Position of Local Environment ----------------
#Hpos = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0],[10.0,11.0,12.0]]
#Hpos = genfromtxt('../Structs/mos2H.xyz').tolist()
Hpos = genfromtxt('../Structs/au40cu40H.dat').tolist()

#-------------- set Basis Function (rCut--soft, N_max) Environment ----------------
myAlphas, myBetas = genBasis.getBasisFunc(10.0, 5)


#-------------- run local chemical environments on each atom ----------------
features = soapPy.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 
orig_atoms = atoms.copy()
# rotation check
for rotation in ['x', 'y', 'z']:
    print("rotating in", rotation)
    atoms.rotate(45, rotation)
    rot_features = soapPy.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

    deviation = np.max(features- rot_features)
    print("maximal numerical deviation:", deviation)
    if deviation > 10e-8:
        IS_PASS = False


# translation check
atoms = orig_atoms.copy()
for translation in [[1.0, 1.0, 1.0], [-5.0, 5.0, -5.0], [1.0, 1.0, -10.0],]:
    print("translating towards", translation)
    atoms.translate(translation)
    trans_features = soapPy.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

    deviation = np.max(features- trans_features)
    print("maximal numerical deviation:", deviation)
    if deviation > 10e-9:
        IS_PASS = False

#symmetry check
print("check with molecules")
print(g2.names)
atoms = molecule('SiH4')
features = soapPy.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

#deviation = pdist(features[:4])
print(deviation)
for idx in np.arange(3):
    deviation = np.max(features[idx] - features[idx + 1])
    print("maximal numerical deviation:", deviation)
    if deviation > 10e-8:
        IS_PASS = False


#------------ Check if test passed --------------------------------
if IS_PASS:
    print("test passed")
else:
    print("WARNING, test failed")
