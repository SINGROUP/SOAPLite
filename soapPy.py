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

    typeNs = atomtype_lst
    Ntypes = len(n_atoms_per_type_lst)
    Nsize = n_atoms_per_type_lst
    Apos = np.concatenate(pos_lst).ravel()
    return Apos, typeNs, Ntypes, Nsize

### INPUT ###

parser = argparse.ArgumentParser(description='Give xyzfile')
parser.add_argument('arguments', metavar='args', type=str, nargs='+', help='xyzfile')
args = parser.parse_args()
print("Passed arguments:", args.arguments)
if len(args.arguments) < 1:
    print('Not enough arguments')
    exit(1)
if len(args.arguments) > 1:
    print('Too many arguments')
    exit(1)

inpfile = args.arguments[0]
rootname = inpfile.replace(".xyz", "")

### PROCESS ###
structure = ase.io.read(inpfile, index = -1)

Apos, typeNs, Ntypes, Nsize = format_ase2clusgeo(structure)
print("Apos")
print(Apos)
print("typeNs")
print(typeNs)
print("Ntypes")
print(Ntypes)
print("Nsize")
print(Nsize)

### START SOAP###
#libsoap = CDLL('./libsoapPy.so')
#libsoap.soap.argtypes = [POINTER (c_double),POINTER (c_double), POINTER (c_double), POINTER (c_int),c_int,c_int,c_int,c_int,c_int]

#libsoap.soap.restype = POINTER (c_double)

def soap(Apos, Hpos, NradBas=4, Lmax=5):

    types = set()
    for i in Apos: types.add(i[0])
    Ntypes = len(types)

    c = (c_double*(NradBas*NradBas*(Lmax+1)*Ntypes*len(Hpos)))
    axyz = (c_double*(len(Apos)*3))()
    hxyz = (c_double*(len(Hpos)*3))()

    for i in range(len(Apos)):
        axyz[3*i + 0] = c_double(Apos[i][1])
        axyz[3*i + 1] = c_double(Apos[i][2])
        axyz[3*i + 2] = c_double(Apos[i][3])
    for i in range(len(Hpos)):
        hxyz[3*i + 0] = c_double(Hpos[i][1])
        hxyz[3*i + 1] = c_double(Hpos[i][2])
        hxyz[3*i + 2] = c_double(Hpos[i][3])
    #   c = libsoap.soap(c, axyz,hxyz,totalAN, Ntypes, Nsize, l, Hsize)
    #   soap(axyz,typeNs.data_as(c_double_p),totalAN, Ntypes, Nsize, l, Hsize))
    # transform into c_type types
    # axyz
    # hxyz
    """
    typeNs
    totalAN =
    Ntypes
    Nsize
    l
    Hsize
    """

    soap(axyz,hxyz,typeNs,totalAN, Ntypes, Nsize, l, Hsize)
    numArrya = array(c)
#   return c;
    return numpy.ctypeslib.as_array(c, shape=(Hsize,NradBas*NradBas*(lmax+1)*Ntypes))
