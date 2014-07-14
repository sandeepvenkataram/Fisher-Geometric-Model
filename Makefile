CC=g++
CFLAGS=-lgsl -lgslcblas
DEPS = src/allele.h src/environment.h src/modelFunctions.h src/population.h src/randomv.h src/runningStat.h 

a:	src/tester.cpp src/allele.cpp src/environment.cpp src/modelFunctions.cpp src/population.cpp src/randomv.cpp src/runningStat.cpp src/stability.cpp
	$(CC) -O3 -o bin/a.out src/tester.cpp src/allele.cpp src/environment.cpp src/modelFunctions.cpp src/population.cpp src/randomv.cpp src/runningStat.cpp src/stability.cpp $(CFLAGS)
geometry: src/geometry.cpp src/randomv.cpp src/modelFunctions.cpp src/environment.cpp
	$(CC) -O3 -o bin/geometry.exe src/geometry.cpp src/randomv.cpp src/modelFunctions.cpp src/environment.cpp $(CFLAGS)
clean: 
	rm -f bin/a.out bin/geometry.exe