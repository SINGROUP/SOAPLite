from setuptools import setup, find_packages, Extension

ext_list = []
for extname, soname in zip(
    [
        "src/soapAnalFullPySigma.c",
        "src/soapGTO.c",
        "src/soapGeneral.c",
    ],
    [
        "lib.libsoapPySig",
        "lib.libsoapGTO",
        "lib.libsoapGeneral",
    ]):
    ext_list.append(Extension(soname,
        [extname],
        include_dirs=["src"],
        extra_compile_args=["-O3", "-std=c99"]
    ))


extensions = ext_list

if __name__ == "__main__":
    setup(
        name="soaplite",
        url="https://github.com/SINGROUP/SOAPLite",
        version="1.0.3",
        description=("fast lightweight smooth overlap atomic position (SOAP) calculator. see github.com/SINGROUP/SOAPLite for detailed documentations."), author="Eiaki V. Morooka", author_email="eiaki.morooka@aalto.fi",
        packages=find_packages(),
        install_requires=[
            "numpy",
            "scipy",
            "future",
            "ase"
        ],
        python_requires='>=2.2, <4',
        keywords="smooth overlap of atomic positions materials science machine learning soap descriptor fast analytic quantum chemistry local environment materials physics symmetry reduction",
        license="LGPLv3",
        classifiers=[
            'Topic :: Scientific/Engineering :: Physics',
            'Operating System :: POSIX :: Linux',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Artificial Intelligence',
            'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
            'Development Status :: 4 - Beta', 'Intended Audience :: Science/Research',
            'Intended Audience :: Religion',
            'Intended Audience :: Education',
            'Intended Audience :: Developers',
            'Programming Language :: Python',
            'Programming Language :: C'
        ],
        ext_modules=extensions,
        package_data={
            'soaplite': [
                'tests/convergenceTest.py',
                'tests/example_non_periodic_SoapPy.py',
                'tests/example_periodic_SoapPy.py',
                'tests/test_symmetry.py',
                'utilities/batchSoapPy.py',
                'tests/Structs/au40cu40H.dat',
                'tests/Structs/au40cu40.xyz',
                'tests/Structs/Cu_110.pdb',
                'tests/Structs/h2oDiff.xyz',
                'tests/Structs/h2o.xyz',
                'tests/Structs/mos2_51.xyz',
                'tests/Structs/mos2H.dat',
            ]
        }
    )
