omp:
	gcc -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalFullPy.c -fopenmp
serial:
	gcc -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalFullPy.c 
