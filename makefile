all:
	gcc -fPIC -O3 -shared -o src/libsoapGTO.so src/soapGTO.c -lm 
	gcc -fPIC -O3 -shared -o src/libsoapPy.so src/soapAnalFullPy.c -lm 
	gcc -fPIC -O3 -shared -o src/libsoapPy2.so src/soapAnalFullCro2Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy3.so src/soapAnalFullCro3Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy4.so src/soapAnalFullCro4Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy5.so src/soapAnalFullCro5Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy6.so src/soapAnalFullCro6Py.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPysig.so src/soapAnalFullPySigma.c -lm 
	gcc -fPIC -O3 -shared -o src/libsoapPy2sig.so src/soapAnalFullCro2PySigma.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy3sig.so src/soapAnalFullCro3PySigma.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy4sig.so src/soapAnalFullCro4PySigma.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy5sig.so src/soapAnalFullCro5PySigma.c -lm
	gcc -fPIC -O3 -shared -o src/libsoapPy6sig.so src/soapAnalFullCro6PySigma.c -lm

half:
	gcc -fPIC -O3 -shared -o src/libsoapPy2.so src/soapAnalFull2Half.c -lm
