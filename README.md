# SOAP LIGHT 

One Paragraph of project description goes here

## Getting Started

SOAP is an algorithm used for accurately classifying and machine learning chemical environments[].
This is a very light weight and fast SOAP implementation for machine learning in quantum chenistry and materials physics. Once you give SOAP the .xyz in ASE format,  list of positions, radius cutoff, number of basis functions and l (number of spherical harmonics), soap will return a numpy matrix of the power spectra for each point. Each row corresponds to each specified point, and each column corresponds the the spectrum value.


### Prerequisites

Numpy and Scipy are required. To install them by typing in the terminal: 

```
sudo pip install numpy
```
and
```
sudo pip install scipy
```
If you don't have super user access, install them by
```
pip install numpy --user
```
and
```
pip install scipy --user
```
You also need a gcc and g++ compiler.


### Installing

Simply type into the terminal:
```
make
```
if you have gcc and g++ compilers.

## Running the tests

Run 
```
python testSoapPy.py
```
This should generate a test.txt file which contains the SOAP spectra for every H.dat point of a au40cu40 molecule (in au40cu40.xyz).


## Authors

* **Eiaki V. Morooka** - [Aki78]( https://github.com/Aki78)
* **Marc Jäger** - [marchunter](https://github.com/marchunter)

See also the list of contributors who participated in this project.
* **Yu Ninomiya** - [lab](http://www.sp.u-tokai.ac.jp/~bentz/Members.html)
* **Filippo Federici** - [fullmetalfelix](https://github.com/fullmetalfelix)
* **Yashasvi S. Ranawat** - [yashasvi-ranawat](https://github.com/yashasvi-ranawat)

## License

This project is licensed under the GNU LESSER GENERAL PUBLIC LICENSE - see the [LICENSE.md](LICENSE.md) file for details


