# SOAP LIGHT 

Smooth Overlap of Atomic Positions (SOAP) is an algorithm used for accurately classifying and machine learning chemical environments[1,2].


## Getting Started

This is a very light weight and fast SOAP implementation for machine learning in quantum chenistry and materials physics. Once you give SOAP the .xyz in ASE format,  list of positions, radius cutoff, number of basis functions and l (number of spherical harmonics), soap will return a numpy matrix of the power spectrum for each point. Each row corresponds to each specified point, and each column corresponds the the spectrum value.


### Prerequisites

Numpy and Scipy and ASE are required. To install them by typing in the terminal: 

```
sudo pip install numpy
```
```
sudo pip install scipy
```
and
```
sudo pip install ase
```
If you don't have super user access, install them by
```
pip install numpy --user
```
```
pip install scipy --user
```
and
```
pip install ase --user
```
You also need a gcc compiler.


### Installing

Simply type into the terminal:
```
make
```
if you have gcc compilers.

## Running the tests

Run 
```
python testSoapPy.py
```
This should generate a test.txt file which contains the SOAP spectrum for every H.dat point of a au40cu40 molecule (in au40cu40.xyz).

## Possible Applications 

By taking the differences of the soap spectrum, we can compare the differences of the chemical environment. For example, if a point P1 gave a power
spectrum S1 and at point P2 gave  a spectrum S2, the difference of the chemical environment will be |S2 - S1| where || denotes the Euclidean distance.
We can use this differences to classify similar/different chemical environments.

The power spectrum can also be used as an input for a neural network, kernel ridge regression or other machine learning algorithms.

## Authors

* **Eiaki V. Morooka** - [Aki78]( https://github.com/Aki78)
* **Marc Jäger** - [marchunter](https://github.com/marchunter)

See also the list of contributors who participated in this project.
* **Yu Ninomiya** - [Yu](http://www.sp.u-tokai.ac.jp/~bentz/Members.html)
* **Filippo Federici** - [fullmetalfelix](https://github.com/fullmetalfelix)
* **Yashasvi S. Ranawat** - [yashasvi-ranawat](https://github.com/yashasvi-ranawat)
* **Lauri Himanen** - [lauri-codes](https://github.com/lauri-codes)
* **Adam Foster** - [suurimonster](https://github.com/suurimonster)


## License

This project is licensed under the GNU LESSER GENERAL PUBLIC LICENSE - see the [LICENSE.md](LICENSE.md) file for details

## References
* [1] On representing chemical environments  - Albert P. Bartók, Risi Kondor, Gábor Csányi [paper](https://arxiv.org/abs/1209.3140)
* [2] Comparing molecules and solids across structural and alchemical space -  Sandip De, Albert P. Bartók, Gábor Cásnyi, and Michele Ceriotti [paper](https://arxiv.org/pdf/1601.04077.pdf)

