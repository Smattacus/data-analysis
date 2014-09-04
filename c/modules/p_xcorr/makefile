rebuild: p_xcorr.c p_xcorr.h
		rm p_xcorr.o libxcorr.so
			gcc -Wall -fPIC -fopenmp -c p_xcorr.c
				gcc -shared -fopenmp -o libxcorr.so p_xcorr.o

debug:
	gcc -Wall -fPIC -g -c p_xcorr.c
		gcc -shared -g -o libxcorr.so p_xcorr.o

build: p_xcorr.c p_xcorr.h
		gcc -Wall -fPIC -c p_xcorr.c
			gcc -shared -o libxcorr.so p_xcorr.o
buildP: p_xcorr.c p_xcorr.h
		gcc -Wall -fPIC -fopenmp -c p_xcorr.c
		gcc -shared -fopenmp -o libxcorr.so p_xcorr.o


clean: 
		rm p_xcorr.o libxcorr.so

buildO1: p_xcorr.c p_xcorr.h
	gcc -Wall -fPIC -fopenmp -O1 -c p_xcorr.c
	gcc -shared -fopenmp -O1 -o libxcorr.so p_xcorr.o

buildO2: p_xcorr.c p_xcorr.h
	gcc -Wall -fPIC -fopenmp -O2 -c p_xcorr.c
	gcc -shared -fopenmp -O2 -o libxcorr.so p_xcorr.o

buildO3: p_xcorr.c p_xcorr.h
	gcc -Wall -fPIC -fopenmp -O3 -c p_xcorr.c
	gcc -shared -fopenmp -O3 -o libxcorr.so p_xcorr.o
