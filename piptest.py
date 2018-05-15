from soaplite import get_soap_structure
from soaplite import getBasisFunc
from ase.build import molecule

atoms = molecule('CO2')

a,b = getBasisFunc(5.0, 10)
x = get_soap_structure(atoms, a, b)


print(x.shape)
print(x)
