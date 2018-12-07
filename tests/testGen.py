import soaplite.getBasis
import soaplite
import numpy as np
import matplotlib.pyplot as plt
import ase

alphas = np.loadtxt("alphasPy.dat");

#Define cutoff and desired basis functions (normalized or not normalized)
NMax, rx, gss = soaplite.getBasis.getGns(
      5.0, 5, # rCut=0.5 is only for the integration cutoff, you have to design the basis function so that it cuts off before rCut.
      [lambda x: x*x*np.exp(-alphas[0]*x*x), # Here, input your desired basis function
       lambda x: x*x*np.exp(-alphas[1]*x*x), # It will Normalize them automatically.
       lambda x: x*x*np.exp(-alphas[2]*x*x), # It can be as many basis functions as desired
       lambda x: x*x*np.exp(-alphas[3]*x*x), # but with minimum 2
       lambda x: x*x*np.exp(-alphas[4]*x*x)] # 
      );

#Define xyz structure
atoms = ase.io.read("Structs/au40cu40.xyz")

#Define positions
Hpos = np.genfromtxt('Structs/au40cu40H.dat') # Must be a list of a list

#Run soap calculations
x = soaplite.get_soap_locals_general(atoms, Hpos, rx, gss,  rCut=5.0, nMax=NMax, Lmax=2, all_atomtypes=[], eta=1.0)
np.savetxt("testme.txt",x)

atoms = ase.io.read("Structs/au40cu40.xyz")
Hpos=[[0,0,0]]
x = soaplite.get_soap_locals_general(atoms, Hpos, rx, gss,  rCut=5.0, nMax=NMax, Lmax=5, all_atomtypes=[], eta=1.0)

