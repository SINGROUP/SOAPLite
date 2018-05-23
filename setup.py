from setuptools import setup, find_packages, Distribution

class BinaryDistribution(Distribution):
  def has_ext_modules(foo):
    return True

if __name__=="__main__":
  setup(name="soaplite",
      url="https://github.com/SINGROUP/SOAPLite",
      version="0.9.6",
      description=("fast lightweight smooth overlap atomic position (SOAP) calculator. see github.com/SINGROUP/SOAPLite for detailed documentations."), author="Eiaki V. Morooka", author_email="eiaki.morooka@aalto.fi",
      packages = find_packages(),
      install_requires =["numpy",
      "scipy",
      "ase"], python_requires='>=2.2, <4', keywords="smooth overlap of atomic positions materials science machine learning soap descriptor fast analytic quantum chemistry local environment materials physics symmetry reduction",
      license="LGPLv3", classifiers=['Topic :: Scientific/Engineering :: Physics', 'Operating System :: POSIX :: Linux' ,'Topic :: Scientific/Engineering :: Chemistry','Topic :: Scientific/Engineering :: Artificial Intelligence','License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)', 'Development Status :: 4 - Beta', 'Intended Audience :: Science/Research','Intended Audience :: Religion', 'Intended Audience :: Education','Intended Audience :: Developers','Programming Language :: Python','Programming Language :: C' ],
      package_data={
        'soaplite':['src/libsoapPy2.so','src/libsoapPy3.so','src/libsoapPy4.so','src/libsoapPy5.so','src/libsoapPy6.so','src/libsoapPy.so']})
      #      ,distclass=BinaryDistribution)
