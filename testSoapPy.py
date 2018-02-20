import soapPy
import ase
import numpy as np

# example structure
#atoms = ase.io.read("teststr.xyz")

h2o = ase.Atoms('H2O', positions=[[0.0, 0.1, 0.2 ], [3.2,4.5,5.5],[5.0, 6.0, 7.0]])

#Apos = [[[0.0, 0.1, 0.2 ], [3.2,4.5,5.5]],[[5.0, 6.0, 7.0]]]
Hpos = [[2.3, 4.1 ,3.5],[2.2,3.3,4.4]]

x = soapPy.soap(h2o, Hpos, 4, 5)

print("SOAP solution: ")
print(x)
print("soap size: ", np.shape(x))
print("type: ", type(x))


"""
# old input from soapPy.py, in case it is still needed for debugging
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

Apos, typeNs, Ntypes, Nsize, totalAN = format_ase2clusgeo(structure)
print("Apos")
print(Apos)
print("typeNs")
print(typeNs)
print("Ntypes")
print(Ntypes)
print("Nsize")
print(Nsize)

Hpos = [[2.3, 4.1 ,3.5],[2.2,3.3,4.4]]
soap(structure, Hpos)
"""
