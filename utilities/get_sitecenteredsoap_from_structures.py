import numpy as np
import pandas as pd
import ase.io
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.spatial import ConvexHull
import argparse, time
import soapPy, genBasis

### GLOBALS ###

NN_CUTOFF = 2.0
feature_lst = []

# helper functions

def get_distance_matrix(obj):
    last_atom = obj.positions[-1]
    distmat = np.linalg.norm(last_atom - obj.positions[0:-1], axis = 1)
    return distmat

def get_nearest_neighbour(obj, distmat):
    nn_id = np.argmin(distmat)
    distance = distmat[nn_id]
    nn_type = obj.get_chemical_symbols()[nn_id]
    return distance, nn_type, nn_id

def get_type_info(obj, idx = True):
    nn_type = obj.get_atomic_numbers()[idx]

    atomtype_set = set(nn_type)
    num_atomtypes = len(atomtype_set)

    atomtype_lst = np.sort(list(atomtype_set))
    n_atoms_per_type_lst = []
    for atomtype in atomtype_lst:
        condition = obj.get_atomic_numbers() == atomtype
        nn_ids = np.zeros(len(condition))
        nn_ids[idx] = 1
        
        n_onetype = len(obj.get_atomic_numbers()[(condition) & (nn_ids == 1)])
        n_atoms_per_type_lst.append(n_onetype)
    return nn_type, atomtype_lst, n_atoms_per_type_lst 

def get_first_shell(obj, shelldist, distmat):
    first_shell= np.where(distmat < shelldist)
    first_shell = first_shell[0]
    #print(first_shell)
    nn_dist_lst = distmat[first_shell]
    ids = np.argsort(nn_dist_lst)
    #print(ids, ids.shape, ids[0], type(ids[0]))
    nn_dist_lst = nn_dist_lst[ids] #np.sort(nn_dist_lst)
    first_shell = first_shell[ids]

    nn_type, atomtype_lst, n_atoms_per_type_lst = get_type_info(obj, idx = first_shell)
    return first_shell, nn_dist_lst, nn_type, n_atoms_per_type_lst

def get_nnsoap(obj, first_shell, alphas, betas, rcut=6, nmax=10, lmax=9, all_atomtypes=[]):
    """Takes cluster structure and nearest neighbour information of a datapoint,
    Returns concatenated soap vectors for each nearest
    neighbour (up to 3). Top, bridge, hollow fill the initial
    zero soap vector from left to right.
    """

    soap_vector = []
    nnn = len(first_shell)
    for tbh in range(0,3):
        try:
            atom_idx = first_shell[tbh]
        except:
            soap_vector.append(soap_zero)
        else:
            Hpos = []
            print(atom_idx)
            pos = obj.get_positions()[atom_idx]
            Hpos.append(pos)
            x = soapPy.get_soap_locals(obj, Hpos, myAlphas, myBetas, rCut=rcut, NradBas=nmax, Lmax=lmax,crossOver=False, all_atomtypes=all_atomtypes)
            soap_zero = np.zeros(x.shape)
            soap_vector.append(x)
    print(len(soap_vector), soap_vector[0].shape, soap_vector[1].shape, soap_vector[2].shape)
    print("exemplary soapvalues",soap_vector[0][0,1], soap_vector[1][0,1], soap_vector[2][0,1])
    soap_array = np.hstack(soap_vector)
    return soap_array

def get_sitecenteredsoap(obj, first_shell, alphas, betas, rcut=6, nmax=10, lmax=9, all_atomtypes=[]):
    """Takes cluster structure and nearest neighbour information of a datapoint,
    Returns concatenated soap vectors for each nearest
    neighbour (up to 3). Top, bridge, hollow fill the initial
    zero soap vector from left to right.
    """

    soap_vector = []
    nnn = len(first_shell)
    center_of_atoms = 1.0 / nnn * np.mean(obj.get_positions()[first_shell], axis = 0)
    #print("center of atoms", center_of_atoms)
    Hpos = [center_of_atoms]
    soap_vector = soapPy.get_soap_locals(obj, Hpos, myAlphas, myBetas, rCut=rcut, NradBas=nmax, Lmax=lmax,crossOver=False, all_atomtypes=all_atomtypes)
    #print(len(soap_vector), soap_vector.shape)
    #print("exemplary soapvalues",soap_vector[0,0], soap_vector[0,1], soap_vector[0,2])
    return soap_vector

### INPUT ###
parser = argparse.ArgumentParser(description='Give xyzfile')
parser.add_argument('arguments', metavar='args', type=str, nargs='+',
                                   help='xyzfile in trajectory format')
args = parser.parse_args()
print("Passed arguments:", args.arguments)
if len(args.arguments) < 1:
    print('Not enough arguments')
    exit(1)  
if len(args.arguments) > 1:
    print('Too many arguments')
    exit(1)

trajfile = args.arguments[0]
rootname = trajfile.replace(".xyz", "") 

traj = ase.io.read(trajfile, index= ":")

print("number of snapshots:", len(traj))

feature_mat = np.zeros((len(traj), 21))
NMAX = 10
LMAX = 9
RCUT = 8.0
myAlphas, myBetas = genBasis.getBasisFunc(RCUT, NMAX)
n_datapoints = len(traj)


for i, atoms in zip(np.arange(len(traj)),traj):
    print(i)
    cluster = atoms[:-1]
    last_atom = atoms[-1]
    distmat = get_distance_matrix(atoms)

    # get number of nearest neighbours, classify as top, bridge, hollow, 4-fold hollow
    first_shell, nn_dist_lst, nn_type, n_atoms_per_type_lst = get_first_shell(atoms, NN_CUTOFF, distmat)

    nnn = len(nn_dist_lst)
    print("number of nearest neighbours", nnn)
    soap_array = get_sitecenteredsoap(cluster, first_shell, myAlphas, myBetas, RCUT, NMAX, LMAX, all_atomtypes=[29,79])
    #print(soap_array.shape)

    if i == 0:
        n_features = soap_array.shape[1]
        print("soap first", soap_array.shape)
        print(n_datapoints, n_features)
        soapmatrix = np.zeros((n_datapoints, n_features))


    print("Processing " + str(atoms.info), end="\r")
    soapmatrix[i,:] = soap_array = get_sitecenteredsoap(cluster, first_shell, myAlphas, myBetas, RCUT, NMAX, LMAX, all_atomtypes=[29,79])
    print("")

# infos
print("shape", soapmatrix.shape)
np.save(rootname + "_citecenteredsoap_" + str(NMAX) + "L" + str(LMAX) + "R" + str(RCUT) + ".npy", soapmatrix)
