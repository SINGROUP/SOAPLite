import numpy as np 
import ase
from scipy.spatial.distance import cdist
import soaplite
from ase.build import molecule

def compute_envkernel(local_a,local_b):
    return cdist(local_a, 
            local_b, 
            lambda u, v: np.exp(-np.sqrt(((u-v)**2).sum())),
            )



def rematch_entry(envkernel, gamma = 0.1, threshold = 1e-6):
    """
    Compute the global similarity between two structures A and B.
    It uses the Sinkhorn algorithm as reported in:
    Phys. Chem. Chem. Phys., 2016, 18, p. 13768
    Args:
        envkernel: NxM matrix of structure A with 
            N and structure B with M atoms
        gamma: parameter to control between best match gamma = 0
            and average kernel gamma = inf.
    """
    n, m = envkernel.shape
    K = np.exp(-(1 - envkernel) / gamma)

    # initialisation
    u = np.ones((n,)) / n
    v = np.ones((m,)) / m

    en = np.ones((n,)) / float(n)
    em = np.ones((m,)) / float(m)

    Kp = (1 / en).reshape(-1, 1) * K

    # converge balancing vectors u and v
    itercount = 0
    error = 1
    while (error > threshold):
        uprev = u
        vprev = v
        v = np.divide(em, np.dot(K.T, u))
        u = np.divide(en, np.dot(K, v))

        # determine error every now and then
        if itercount % 5:
            error = np.sum((u - uprev) ** 2) / np.sum((u) ** 2) + np.sum((v - vprev) ** 2) / np.sum((v) ** 2)
        itercount += 1

    # using Tr(X.T Y) = Sum[ij](Xij * Yij)
    # P.T * C
    # P_ij = u_i * v_j * K_ij
    pity = np.multiply( np.multiply(K, u.reshape((-1,1))) , v)

    glosim = np.sum( np.multiply( pity, envkernel))

    return glosim

def test_glosim_molecules():
    is_pass = True

    # check if the same molecules give global similarity of around 1
    from ase.collections import g2
    myAlphas, myBetas = soaplite.genBasis.getBasisFunc(10.0, 5)
    for molname in g2.names:
        atoms = molecule(molname)
        local_a = soaplite.get_soap_structure(atoms, myAlphas, 
            myBetas, rCut=10.0, NradBas=5, Lmax=5,crossOver=True)

        envkernel = compute_envkernel(local_a, local_a)
        glosim = rematch_entry(envkernel, gamma = 0.01)

        if (0.99 < glosim < 1.01):
            pass
        else:
            is_pass = False

    # check randomly a few combinations of molecules, just for no errors
    all_atomtypes = [1,6,7,8]
    for molname1 in g2.names:
        for molname2 in g2.names:
            if np.random.rand(1,1) < 0.99:
                continue
            atoms1 = molecule(molname1)        
            atoms2 = molecule(molname2)

            local_a = soaplite.get_soap_structure(atoms1, myAlphas, 
                myBetas, rCut=10.0, NradBas=5, Lmax=5,crossOver=True, 
                all_atomtypes = all_atomtypes)
            local_b = soaplite.get_soap_structure(atoms2, myAlphas, 
                myBetas, rCut=10.0, NradBas=5, Lmax=5,crossOver=True,
                all_atomtypes = all_atomtypes)

            envkernel = compute_envkernel(local_a, local_b)
            glosim = rematch_entry(envkernel, gamma = 0.01)
            #print(glosim)

    return is_pass


if __name__ == "__main__":
    print("Start")
    myAlphas, myBetas = soaplite.genBasis.getBasisFunc(10.0, 5) # input:(rCut, NradBas)

    A = molecule('H2O')
    B = molecule('H2O2')
    local_a = soaplite.get_soap_structure(A, myAlphas, 
            myBetas, rCut=10.0, NradBas=5, Lmax=5,crossOver=True) 
    local_b = soaplite.get_soap_structure(B, myAlphas, 
            myBetas, rCut=10.0, NradBas=5, Lmax=5,crossOver=True) 

    envkernel = compute_envkernel(local_a, local_b)

    glosim = rematch_entry(envkernel, gamma = 0.01)

    print("global similarity")
    print(glosim)

    print("same")
    envkernel = compute_envkernel(local_a, local_a)
    glosim = rematch_entry(envkernel, gamma =0.01)
    print(glosim)

    print("unity envkernel")

    envkernel = np.ones((5,5)) * 1.0
    glosim = rematch_entry(envkernel, gamma =0.01)
    print(glosim)

    print("test glosim on different molecules")
    pass_test = test_glosim_molecules()
    print(pass_test)