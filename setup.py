from setuptools import setup, find_packages, Extension

ext_list = []
for extname, soname in zip(
    [
        "src/soapAnalFullPy.c",
        "src/soapAnalFullCro2Py.c",
        "src/soapAnalFullCro3Py.c",
        "src/soapAnalFullCro4Py.c",
        "src/soapAnalFullCro5Py.c",
        "src/soapAnalFullCro6Py.c",
        "src/soapAnalFullPySigma.c",
        "src/soapAnalFullCro2PySigma.c",
        "src/soapAnalFullCro3PySigma.c",
        "src/soapAnalFullCro4PySigma.c",
        "src/soapAnalFullCro5PySigma.c",
        "src/soapAnalFullCro6PySigma.c",
        "src/soapGeneral.c",
    ],
    [
        "lib.libsoapPy",
        "lib.libsoapPy2",
        "lib.libsoapPy3",
        "lib.libsoapPy4",
        "lib.libsoapPy5",
        "lib.libsoapPy6",
        "lib.libsoapPysig",
        "lib.libsoapPy2sig",
        "lib.libsoapPy3sig",
        "lib.libsoapPy4sig",
        "lib.libsoapPy5sig",
        "lib.libsoapPy6sig",
        "lib.libsoapGeneral",
    ]):
    ext_list.append(Extension(soname,
        [extname],
        include_dirs=["src"],
        # libraries=["m"],
        extra_compile_args=["-O3", "-std=c99"]
    ))


extensions = ext_list

if __name__ == "__main__":
    setup(
        name="soaplite",
        url="https://github.com/SINGROUP/SOAPLite",
        version="0.14.5",
        description=("fast lightweight smooth overlap atomic position (SOAP) calculator. see github.com/SINGROUP/SOAPLite for detailed documentations."), author="Eiaki V. Morooka", author_email="eiaki.morooka@aalto.fi",
        packages=find_packages(),
        install_requires=[
            "numpy",
            "scipy",
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
