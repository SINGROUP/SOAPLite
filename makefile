all:
	g++ -fPIC -O3 -shared -o libsoapPy.so soapAnalPy.cpp
