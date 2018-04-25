serial:
	gcc -fPIC -O -shared -o src/libsoapPy.so src/soapAnalFullCroPy.c -lm
noCross:
	gcc -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalFullPy.c -lm 
