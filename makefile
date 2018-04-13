serial:
	gcc -fPIC -O -shared -o src/libsoapPy.so src/soapAnalFullCroPy.c -lm
omp:
	gcc -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalFullCroPy.c -lm -fopenmp
