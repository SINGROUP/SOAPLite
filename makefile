all:
	g++ -O3  soapAnalOptOpt.cpp 
omp:
	g++ -O3  soapAnalOptOpt.cpp -fopenmp
