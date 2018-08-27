import soaplite.getBasis
import soaplite
import numpy as np
import matplotlib.pyplot as plt
import ase


NMax, rx, gss = soaplite.getBasis.getGns(
      [lambda x: np.exp(-x), lambda x: np.exp(-0.5*x), lambda x: np.exp(-0.25*x), lambda x: np.exp(-0.125*x)
   ]);

#print(nMax)
#print(rx)
#print(gss)

for i in range(4):
  plt.plot(gss[i])


atoms = ase.io.read("tests/Structs/au40cu40.xyz")
Hpos=[[0,0,0]]
x = soaplite.get_soap_locals_general(atoms, Hpos, rx, gss, gaussAlpha=1.0, rCut=5.0, nMax=NMax, Lmax=5, all_atomtypes=[])

np.savetxt("testme.txt",x)
