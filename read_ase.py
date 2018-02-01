#!/usr/bin/env python3
# reads an xyzfile
# outputs 3 files:
#   metadatafile with: TOTAL_N_ATOMS N_TYPES
#   typecountfile with: N_ATOMS N_ATOMS ....
#       for each atom type in ascending order
#   posfile with:
#   ATOMIC_NUMBER_1
#   X Y Z
#   ...
#   ATOMIC_NUMBER_2
#   X Y Z
#   ...
#   ATOMIC_NUMBER_3
#   X Y Z
#   ...

### DEFINE ###
import os, argparse
import numpy as np
import ase, ase.io

### INPUT ###

parser = argparse.ArgumentParser(description='Give xyzfile')
parser.add_argument('arguments', metavar='args', type=str, nargs='+',
                                   help='xyzfile')
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

###
#atoms metadata
totnum = len(structure)
atomtype_set = set(structure.get_atomic_numbers())
num_atomtypes = len(atomtype_set)

atomtype_lst = np.sort(list(atomtype_set))
n_atoms_per_type_lst = []
pos_lst = []
for atomtype in atomtype_lst:
    condition = structure.get_atomic_numbers() == atomtype
    pos_onetype = structure.get_positions()[condition]
    n_onetype = pos_onetype.shape[0]

    # store data in lists
    pos_lst.append(pos_onetype)
    n_atoms_per_type_lst.append(n_onetype)


### OUTPUT ###

print("")

# write to files
print("writing files")
metadata = [totnum, num_atomtypes]

print("writing metadata.dat")
with open("metadata.dat", 'w') as f:
    for item in metadata:
        f.write("%s\n" % item)

print("writing atomtypecount.dat")
with open("atomtypecount.dat", 'w') as f:
    for item in n_atoms_per_type_lst:
        f.write("%s\n" % item)

print("writing type_pos.dat")
with open("type_pos.dat", 'w') as f:
    f.write("")

for atomtype, pos_onetype in zip(atomtype_lst, pos_lst):
    with open("type_pos.dat", 'a') as f:
        f.write(str(atomtype) + "\n")

    with open("type_pos.dat", 'ab') as f:
        np.savetxt(f, pos_onetype, fmt='%.8f')

#    with open("type_pos.dat", 'a') as f:
#        f.write("\n")

print("Done")
