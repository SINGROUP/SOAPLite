import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import soaplite
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
print("hello")


#-------------- run local chemical environments on each atom ----------------
features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 
orig_atoms = atoms.copy()
# rotation check
for rotation in ['x', 'y', 'z']:
    print("rotating in", rotation)
    atoms.rotate(45, rotation)
    rot_features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

    deviation = np.max(np.abs(features- rot_features))
    print("maximal numerical deviation:", deviation)
    if deviation > 10e-8:
        IS_PASS = False


# translation check
atoms = orig_atoms.copy()
for translation in [[1.0, 1.0, 1.0], [-5.0, 5.0, -5.0], [1.0, 1.0, -10.0],]:
    print("translating towards", translation)
    atoms.translate(translation)
    trans_features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

    deviation = np.max(np.abs(features- trans_features))
    print("maximal numerical deviation:", deviation)
    if deviation > 10e-9:
        IS_PASS = False

for i in g2.names:
#Rotational Check for many molecules
    atoms = molecule(i)
    print(i)
    features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 
    for rotation in ['x', 'y', 'z']:
        print("rotating in", rotation)
        atoms.rotate(45, rotation)
        rot_features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

        deviation = np.max(np.abs(features- rot_features))
        print("maximal numerical deviation:", deviation)
        if deviation > 10e-8:
            IS_PASS = False


#Translation check for many molecules
    atoms = molecule(i)
    for translation in [[1.0, 1.0, 1.0], [-5.0, 5.0, -5.0], [1.0, 1.0, -10.0],]:
        print("translating towards", translation)
        atoms.translate(translation)
        trans_features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

        deviation = np.max(np.abs(features- trans_features))
        print("maximal numerical deviation:", deviation)
        if deviation > 10e-9:
            IS_PASS = False

#symmetry check
print("check with molecules")
print(g2.names)
atoms = molecule('SiH4')
features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, NradBas=5, Lmax=9,crossOver=True) 

#deviation = pdist(features[:4])
print("deviation: ",deviation)
for idx in np.arange(3):
    deviation = np.max(np.abs(features[idx] - features[idx + 1]))
    print("maximal numerical deviation:", deviation)
    if deviation > 10e-8:
        IS_PASS = False


#------------ Check if test passed --------------------------------
if IS_PASS:
    print("test passed")
else:
    print("WARNING, test failed")
