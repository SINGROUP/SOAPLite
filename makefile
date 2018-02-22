all:
	g++ -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalPy.cpp
