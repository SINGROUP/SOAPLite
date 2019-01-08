all:
	gcc -fPIC -O3 -shared -o src/libsoapGTO.so src/soapGTO.c -lm 
	gcc -fPIC -O3 -shared -o src/libsoapPysig.so src/soapAnalFullPySigma.c -lm 

half:
	gcc -fPIC -O3 -shared -o src/libsoapPy2.so src/soapAnalFull2Half.c -lm
