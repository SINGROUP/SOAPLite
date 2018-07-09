import numpy as np
from soaplite import getBasisFunc, get_periodic_soap_structure
from ase.build import bulk
from scipy.spatial.distance import cdist
import ase.io
IS_PASS = True

atoms = bulk('NaCl', 'rocksalt', a=5.64)
suceatoms = atoms * (2,2,2)

n_max = 5
l_max = 5
r_cut = 10.0
my_alphas, my_betas = getBasisFunc(r_cut, n_max)

x = get_periodic_soap_structure(
    suceatoms,
    my_alphas,
    my_betas,
    rCut=r_cut,
    NradBas=n_max,
    Lmax=l_max,
    crossOver=True
)

#ase.io.write("rocksalt_suce.xyz", suceatoms)

print(x.shape)
print(suceatoms.get_atomic_numbers())


print("sodium soap")
na = x[0::2]
print(na.shape)
deviation = np.max(cdist(na, na), )


print("maximum dissimilarity", deviation)

if deviation > 10e-8:
    IS_PASS = False


if IS_PASS:
    print("test passed")
else:
    print("WARNING, test failed")
