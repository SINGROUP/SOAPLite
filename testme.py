import soaplite

from ase import Atoms

import numpy as np

import math

from soaplite import genBasis



orig_atoms = Atoms(

    cell=[

        [1.0, 0.0, 0.0],

        [0.0, 1.0, 0.0],

        [0.0, 0.0, 1.0]

    ],

    positions=[

        [0, 0, 0],

        [0.95, 0, 0],

        [0.95*(1+math.cos(76/180*math.pi)), 0.95*math.sin(76/180*math.pi), 0.0]

    ],

    symbols=["H", "O", "H"],

)



# Test the GTO basis

myAlphas, myBetas = genBasis.getBasisFunc(10.0, 5)

atoms = orig_atoms.copy()

features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, nMax=5, Lmax=9, crossOver=True)

for rotation in ['x', 'y', 'z']:

    # print("rotating in", rotation)

    atoms.rotate(45, rotation)

    rot_features = soaplite.get_soap_structure(atoms, myAlphas, myBetas, rCut=10.0, nMax=5, Lmax=9, crossOver=True)



    deviation = np.max(np.abs(features - rot_features))

    print("maximal numerical deviation:", deviation)

    if deviation > 10e-8:

        IS_PASS = False



# Test the poly basis

atoms = orig_atoms.copy()

features = soaplite.get_soap_structure_poly(atoms, rCut=10.0, nMax=5, Lmax=9)

for rotation in ['x', 'y', 'z']:

    # print("rotating in", rotation)

    atoms.rotate(45, rotation)

    rot_features = soaplite.get_soap_structure_poly(atoms, rCut=10.0, nMax=5, Lmax=9)



    deviation = np.max(np.abs(features - rot_features))

    print("maximal numerical deviation:", deviation)

    if deviation > 10e-8:

        IS_PASS = False
