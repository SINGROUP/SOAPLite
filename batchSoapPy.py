import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import soapPy
import genBasis
import time
import datetime
import ase
import numpy as np
import argparse

def get_lastatom_soap(atoms, cutoff, myAlphas, myBetas):
    hardcutoff = cutoff + 3
    lastatom = atoms[-1]
    Hpos = [lastatom.position]
    structure = atoms[:-1]
    x = soapPy.get_soap_locals(structure, Hpos, myAlphas, myBetas, rCutHard=hardcutoff, NradBas=10, Lmax=9)
    return x

def create(atoms_list, cutoff = 0, affix="",):
    """Takes a trajectory xyz file and writes soap features
    """
    hardcutoff = cutoff + 3
    myAlphas, myBetas = genBasis.getBasisFunc(hardcutoff, 10)
    # get information about feature length
    n_datapoints = len(atoms_list)
    atoms = atoms_list[0]
    x = get_lastatom_soap(atoms, cutoff, myAlphas, myBetas)
    n_features = x.shape[1]
    print("soap first", x.shape)
    print(n_datapoints, n_features)
    soapmatrix = np.zeros((n_datapoints, n_features))

    i = -1
    for atoms in atoms_list:
        i +=1
        #atoms
        print("Processing " + str(atoms.info),
         " Run time: " + str(time.time()-t0_total),
        end="\r")
        soapmatrix[i,:] = get_lastatom_soap(atoms, cutoff, myAlphas, myBetas)
    print("")

    # infos
    print("shape", soapmatrix.shape)

    # add cutoff info
    if cutoff != 0:
        cutoffstr = "_" + str(cutoff)
    else:
        cutoffstr = ""

    # write descriptor or predictor
    np.save("soap" + cutoffstr + affix + ".npy", soapmatrix)
    return None

##########################################################################

if __name__ == '__main__':
    t0_total = time.time()

    ### INPUT ###
    parser = argparse.ArgumentParser(description=
        'Give input filename trajectory')
    parser.add_argument('arguments', metavar='args', type=str, nargs='+',
                                   help='[filename]')
    args = parser.parse_args()
    print("Passed arguments:", args.arguments)
    if len(args.arguments) < 1:
        print('Not enough arguments')
        exit(1)  
    infilename = args.arguments[0]
    ### PROCESS ###

    atoms_list = ase.io.read(infilename, ':')
    t0 = time.time()
    create(atoms_list, cutoff = 8.0, affix = "")
    t1 = time.time()
    dt = t1 - t0
    t1_total = time.time()
    print("Total run time:", str(t1_total - t0_total))

"""
print("local soap based on given positions")
myAlphas, myBetas = genBasis.getBasisFunc(9.0, 10)
a = datetime.datetime.now()
x = soapPy.get_soap_locals(atoms, Hpos, myAlphas, myBetas, rCutHard=9.0, NradBas=10, Lmax=9) #rCutSoft = rCutHard - 3.0
b = datetime.datetime.now()
c = b - a
print("xyz",c.microseconds)
print("soap for each atom in structure")
#y = soapPy.get_soap_structure(atoms, myAlphas, myBetas, rCutHard=9.0, NradBas=10, Lmax=9) #rCutSoft = rCutHard - 3.0

print("XXX")
# SOAP solution: x
#print("soap size: ", np.shape(x))
np.savetxt('test.txt',x)
np.savetxt('structure_test.txt',y)
"""
