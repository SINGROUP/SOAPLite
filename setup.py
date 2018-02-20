from distutils.core import setup, Extension
import numpy

module = Extension('soapAnalPy', sources = ['soapAnalPy.cpp'], include_dirs=[numpy.get_include()])

setup(name = 'soapAnalPy', version = '1.0', descriptoin = 'my own module', ext_modules = [module])
