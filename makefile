all:
	g++ -O3  soapAnal.cpp 
omp:
	g++ -O3  soapAnal.cpp -fopenmp
