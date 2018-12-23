all:
	gcc -fPIC -O3 -shared -o src/libsoapGTO.so src/soapGTO.c -lm 
	gcc -fPIC -O3 -shared -o src/libsoapPysig.so src/soapAnalFullPySigma.c -lm 
