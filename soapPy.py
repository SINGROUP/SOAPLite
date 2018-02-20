### DEFINE ###
from ctypes import *
import os, argparse
import numpy as np
import ase, ase.io

def format_ase2clusgeo(obj):
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


def soap(obj, Hpos, NradBas=4, Lmax=5):
    # get clusgeo internal format for c-code
    Apos, typeNs, py_Ntypes, atomtype_lst, totalAN = format_ase2clusgeo(obj)
    # flatten Hpos array
    Hpos = np.array(Hpos)
    py_Hsize = Hpos.shape[0]
    Hpos = Hpos.flatten()

    # convert int to c_int
    l = c_int(Lmax)
    Hsize = c_int(py_Hsize)
    Ntypes = c_int(py_Ntypes)
    totalAN = c_int(totalAN)
    Nsize = c_int(NradBas)
    #convert int array to c_int array
    typeNs = (c_int * len(typeNs))(*typeNs)

    print(l, Hsize, Ntypes, totalAN, typeNs, Nsize)
    # convert to c_double arrays
    #Apos
    axyz = (c_double * len(Apos))(*Apos.tolist())
    #Hpos
    hxyz = (c_double * len(Hpos))(*Hpos.tolist())

    ### START SOAP###
    libsoap = CDLL('./libsoapPy.so')
    libsoap.soap.argtypes = [POINTER (c_double),POINTER (c_double), POINTER (c_double), POINTER (c_int),c_int,c_int,c_int,c_int,c_int]
    libsoap.soap.restype = POINTER (c_double)
    # double* c, double* Apos,double* Hpos,int* typeNs,
    # int totalAN,int Ntypes,int Nsize, int l, int Hsize);
    c = (c_double*(NradBas*NradBas*(Lmax+1)*py_Ntypes*py_Hsize))()
    c = libsoap.soap(c, axyz,hxyz,typeNs,totalAN, Ntypes, Nsize, l, Hsize)
    #   return c;
    return np.ctypeslib.as_array(c, shape=(py_Hsize,NradBas*NradBas*(Lmax+1)*py_Ntypes))
