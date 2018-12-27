<p align="center">
  <img src="logoSoapLite.png" height="250">
</p>

# SOAPLite

Smooth Overlap of Atomic Positions (SOAP) is an algorithm used for accurately
classifying and machine learning chemical environments [1,2]. For a detailed
documentation, please read soapDoc.pdf in this repository and visit [DScribe](https://github.com/SINGROUP/dscribe).

## Getting Started

This is a low level, lightweight and fast implementation of SOAP for machine learning in quantum chemistry and materials physics. Inputting structure and SOAP parameters, SOAPLite will spit out the SOAP spectra of local points in space. For a higher level inplementation with kernel methods and Neural Networks, please use [DScribe](https://github.com/SINGROUP/dscribe) instead.

Here is an example of the python interface:
```python
from soaplite import getBasisFunc, get_soap_locals
from ase.build import molecule

#-------------- Define structure -----------------------------------------------
atoms = molecule("H2O")

#-------------- Define positions of desired local environments ----------------
hpos = [
    [0, 1, 2],
    [2, 3, 4]
]

#------------------ Basis function settings (rCut, N_max) ----------------------
n_max = 5
l_max = 5
r_cut = 10.0
my_alphas, my_betas = getBasisFunc(r_cut, n_max)

#--------- Get local chemical environments for each defined position -----------
x = get_soap_locals(
    atoms,
    hpos,
    my_alphas,
    my_betas,
    rCut=r_cut,
    NradBas=n_max,
    Lmax=l_max,
    crossOver=True
)

print(x)
print(x.shape)
```

### Installation

We provide a python interface to the code with precompiled C-extension. This
precompiled version should work with linux-based machines, and can be installed
with:
```
pip install soaplite
```

Or by cloning this repository, you can install it by
```
pip3 install .
```
in the SOAPLite directory.

## Authors

* **Eiaki V. Morooka** - [Aki78]( https://github.com/Aki78)
* **Marc Jäger** - [marchunter](https://github.com/marchunter)
* **Lauri Himanen** - [lauri-codes](https://github.com/lauri-codes)

See also the list of contributors who participated in this project.
* **Yu Ninomiya** - [Yu](http://www.sp.u-tokai.ac.jp/~bentz/Members.html)
* **Filippo Federici** - [fullmetalfelix](https://github.com/fullmetalfelix)
* **Yashasvi S. Ranawat** - [yashasvi-ranawat](https://github.com/yashasvi-ranawat)
* **Adam Foster** - [suurimonster](https://github.com/suurimonster)


## License

This project is licensed under the GNU LESSER GENERAL PUBLIC LICENSE - see the [LICENSE.md](LICENSE.md) file for details

## References
If you use this software, please cite

* [1] On representing chemical environments  - Albert P. Bartók, Risi Kondor, Gábor Csányi [paper](https://arxiv.org/abs/1209.3140)
* [2] Comparing molecules and solids across structural and alchemical space -  Sandip De, Albert P. Bartók, Gábor Cásnyi, and Michele Ceriotti [paper](https://arxiv.org/pdf/1601.04077.pdf)
* Machine learning hydrogen adsorption on nanoclusters through structural descriptors - Marc O. J. Jäger, Eiaki V. Morooka, Filippo Federici Canova, Lauri Himanen & Adam S. Foster   [paper](https://www.nature.com/articles/s41524-018-0096-5)

The theory is based on the first two papers, and the implementation is based on the third paper.
