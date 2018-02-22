omp:
	g++ -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalPy.cpp -fopenmp
serial:
	g++ -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalPy.cpp -fopenmp
