### DEFINE ###
from ctypes import *
import os, argparse
import numpy as np
import genBasis 
import ase, ase.io
import os

def _format_ase2clusgeo(obj):
    """ Takes an ase Atoms object and returns numpy arrays and integers
    which are read by the internal clusgeo. Apos is currently a flattened
    out numpy array
    """
    #atoms metadata
    totalAN = len(obj)

    atomtype_set = set(obj.get_atomic_numbers())
    num_atomtypes = len(atomtype_set)

    atomtype_lst = np.sort(list(atomtype_set))
    n_atoms_per_type_lst = []
    pos_lst = []
    for atomtype in atomtype_lst:
        condition = obj.get_atomic_numbers() == atomtype
        pos_onetype = obj.get_positions()[condition]
        n_onetype = pos_onetype.shape[0]

        # store data in lists
        pos_lst.append(pos_onetype)
        n_atoms_per_type_lst.append(n_onetype)

    typeNs = n_atoms_per_type_lst
    Ntypes = len(n_atoms_per_type_lst)
    atomtype_lst
    Apos = np.concatenate(pos_lst).ravel()
    return Apos, typeNs, Ntypes, atomtype_lst, totalAN

def _get_supercell(obj, rCutHard=8.0):
    """ Takes atoms object (with a defined cell) and a radial cutoff.
    Returns a supercell centered around the original cell
    generously extended to contain all spheres with the given radial
    cutoff centered around the original atoms. 
    """

    xyz_arr = np.abs(np.diag(obj.get_cell()))
    cell_images = np.ceil(rCutHard/xyz_arr)
    nx = int(cell_images[0])
    ny = int(cell_images[1])
    nz = int(cell_images[2])

    suce = obj * (1 + 2*nx, 1+ 2*ny,1+2*nz)
    shift = obj.get_cell()

    shifted_suce = suce.copy()
    shifted_suce.translate(-shift[0]*nx -shift[1]*ny - shift[1]*nz)

    return suce


def get_soap_locals(obj, Hpos, alp, bet, rCutHard=8.0, NradBas=5, Lmax=5):
    print("XDXX")
    assert Lmax <= 9, "l cannot exceed 9. Lmax={}".format(Lmax) 
    assert Lmax >= 0, "l cannot be negative.Lmax={}".format(Lmax) 
    assert rCutHard < 10.0001 , "hard redius cuttof cannot be larger than 10 Angs. rCut={}".format(rCutHard) 
    assert rCutHard > 4.999 , "hard redius cuttof cannot be lower than 5 Ang. rCut={}".format(rCutHard)
    assert NradBas >= 2 , "number of basis functions cannot be lower than 2. NradBas={}".format(NradBas)
    assert NradBas <= 10 , "number of basis functions cannot exceed 10. NradBas={}".format(NradBas)

    # get clusgeo internal format for c-code
    Apos, typeNs, py_Ntypes, atomtype_lst, totalAN = _format_ase2clusgeo(obj)
    # flatten Hpos array
    Hpos = np.array(Hpos)
    py_Hsize = Hpos.shape[0]
    Hpos = Hpos.flatten()

    print("XXXA")
    # convert int to c_int
    alphas = (c_double*len(alp))(*alp)
    betas = (c_double*len(bet))(*bet)
    lMax = c_int(Lmax)
    Hsize = c_int(py_Hsize)
    Ntypes = c_int(py_Ntypes)
    totalAN = c_int(totalAN)
    rCutHard = c_double(rCutHard)
    Nsize = c_int(NradBas)
    #convert int array to c_int array
    typeNs = (c_int * len(typeNs))(*typeNs)

    #print(l, Hsize, Ntypes, totalAN, typeNs, Nsize)
    # convert to c_double arrays
    #Apos
    axyz = (c_double * len(Apos))(*Apos.tolist())
    #Hpos
    hxyz = (c_double * len(Hpos))(*Hpos.tolist())
    print("RXXX")

    ### START SOAP###
    path_to_so = os.path.dirname(os.path.abspath(__file__))
    libsoap = CDLL(path_to_so + '/src/libsoapPy.so')
    libsoap.soap.argtypes = [POINTER (c_double),POINTER (c_double), POINTER (c_double),POINTER (c_double),POINTER (c_double), POINTER (c_int),c_double,c_int,c_int,c_int,c_int,c_int]
    libsoap.soap.restype = POINTER (c_double)
    # double* c, double* Apos,double* Hpos,int* typeNs,
    # int totalAN,int Ntypes,int Nsize, int l, int Hsize);
    c = (c_double*(NradBas*NradBas*(Lmax+1)*py_Ntypes*py_Hsize))()
    c = libsoap.soap( c, axyz, hxyz, alphas, betas, typeNs, rCutHard, totalAN, Ntypes, Nsize, lMax, Hsize)
    print("XXXX")
    #   return c;
    return np.ctypeslib.as_array( c, shape=(py_Hsize,NradBas*NradBas*(Lmax+1)*py_Ntypes))

def get_soap_structure(obj, alp, bet, rCutHard=8.0, NradBas=5, Lmax=5):
    Apos, typeNs, py_Ntypes, atomtype_lst, totalAN = _format_ase2clusgeo(obj)
    Hpos = Apos.copy().reshape((-1,3))
    arrsoap = get_soap_locals(obj, alp, bet, Hpos, rCutHard, NradBas, Lmax)
    return arrsoap

def get_periodic_soap_locals(obj, Hpos, alp, bet, rCutHard=8.0, NradBas=5, Lmax=5):
    # get supercells
    suce = _get_supercell(obj, rCutHard=rCutHard)

    arrsoap = get_soap_locals(suce, Hpos, alp, bet, rCutHard=rCutHard, 
        NradBas=NradBas, Lmax=Lmax)

    return arrsoap

def get_periodic_soap_structure(obj, rCutHard=8.0, NradBas=5, Lmax=5):
    Apos, typeNs, py_Ntypes, atomtype_lst, totalAN = _format_ase2clusgeo(obj)
    Hpos = Apos.copy().reshape((-1,3)) 
    suce = _get_supercell(obj, rCutHard=rCutHard)

    arrsoap = get_soap_locals(suce, Hpos, alp, bet, rCutHard, NradBas, Lmax)
    return arrsoap
