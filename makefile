all:
	gcc -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalFullPy.c -lm 
	gcc -fPIC -O3 -shared -o src/libsoapPy2.so src/soapAnalFullCro2Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy3.so src/soapAnalFullCro3Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy4.so src/soapAnalFullCro4Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy5.so src/soapAnalFullCro5Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy6.so src/soapAnalFullCro6Py.c -lm

half:
	gcc -fPIC -O3 -shared -o src/libsoapPy2.so src/soapAnalFull2Half.c -lm
