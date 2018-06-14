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
x = []
#-------------- Define Structure ------------------------------------
for i in range(2,15):
  myAlphas, myBetas = genBasis.getBasisFunc(5.0, i) # input:(rCut, NradBas)
  atoms1 = ase.io.read("../Structs/h2o.xyz")
  atoms2 = ase.io.read('../Structs/h2oDiff.xyz')

#-------------- set Basis Function (rCut--soft, N_max) Environment ----------------

#-------------- run local chemical environments on each atom ----------------
#start = time.time()
  y1 = soaplite.get_soap_structure(atoms1, myAlphas, myBetas, rCut=5.0, NradBas=i, Lmax=9,crossOver=True) 
  y2 = soaplite.get_soap_structure(atoms2, myAlphas, myBetas, rCut=5.0, NradBas=i, Lmax=9,crossOver=True) 
#endTime = time.time()
#totalTime = endTime - start
  print np.sum((y1 - y2)**2, axis=1)[2] 
#print("Soap ran in seconds:", totalTime)
#np.savetxt('structure_test.txt',y)
#------------ Check Results in *txt files --------------------------------
#------------ Make sure that rCut and N_max are the same across all calculations 
