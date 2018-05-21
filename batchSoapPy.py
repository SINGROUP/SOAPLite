import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import soaplite
import genBasis
import time
import datetime
import ase
import numpy as np
import argparse

def get_lastatom_soap(atoms, cutoff, myAlphas, myBetas, i, j, all_atomtypes=[]):
    lastatom = atoms[-1]
    Hpos = [lastatom.position]
    structure = atoms[:-1]
    x = soapPy.get_soap_locals(structure, Hpos, myAlphas, myBetas, rCut=cutoff, NradBas=i, Lmax=j, crossOver=True)
    return x

def create(atoms_list,N, L, cutoff = 0, all_atomtypes=[]):
    """Takes a trajectory xyz file and writes soap features
    """
    myAlphas, myBetas = genBasis.getBasisFunc(cutoff, N)
    # get information about feature length
    n_datapoints = len(atoms_list)
    atoms = atoms_list[0]
    x = get_lastatom_soap(atoms, cutoff, myAlphas, myBetas,N,L, all_atomtypes=all_atomtypes)
    n_features = x.shape[1]
    print("soap first", x.shape)
    print(n_datapoints, n_features)
    soapmatrix = np.zeros((n_datapoints, n_features))

    i = -1
    for atoms in atoms_list:
        i +=1
        #atoms
        print("Processing " + str(atoms.info)," Run time: " + str(time.time()-t0_total), end="\r")
        soapmatrix[i,:] = get_lastatom_soap(atoms, cutoff, myAlphas, myBetas, N, L, all_atomtypes=all_atomtypes)
    print("")

    # infos
    print("shape", soapmatrix.shape)
    return soapmatrix

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
    rootname = infilename.replace(".xyz", "") 
    ### PROCESS ###

    atoms_list = ase.io.read(infilename, ':')
    t0 = time.time()
    for N in range(2,3):
      for L in range(6,7):
        for cutoff in range(6,7):
            print("N:", N)
#            all_atomtypes = [29,79]
            soapmatrix = create(atoms_list, N, L, cutoff*1.0)
            # write descriptor or predictor
            np.save(rootname +  "_soap_N" + str(N) + "L"+ str(L) +"R" + str(cutoff) + ".npy", soapmatrix)
            print("soaps saved.")
    t1 = time.time()
    dt = t1 - t0
    t1_total = time.time()
    print("Total run time:", str(t1_total - t0_total))

