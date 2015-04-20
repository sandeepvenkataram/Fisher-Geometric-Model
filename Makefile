CC=g++
CFLAGS=-lgsl -lgslcblas
DEPS = src/allele.h src/environment.h src/modelFunctions.h src/population.h src/randomv.h src/runningStat.h 

a:	src/tester.cpp src/allele.cpp src/environment.cpp src/modelFunctions.cpp src/population.cpp src/randomv.cpp src/runningStat.cpp src/stability.cpp
	$(CC) -O3 -o bin/a.out src/tester.cpp src/allele.cpp src/environment.cpp src/modelFunctions.cpp src/population.cpp src/randomv.cpp src/runningStat.cpp src/stability.cpp $(CFLAGS)
stability: src/computeStability.cpp src/stability.cpp
	$(CC) -O3 -o bin/stability src/computeStability.cpp src/stability.cpp $(CFLAGS)
weinreichBackwardPredictability: src/environment.cpp src/modelFunctions.cpp src/randomv.cpp src/computeHeterozygousWeinreichPaths.cpp
	$(CC) -O3 -std=c++0x -ffloat-store -o bin/computeHeterozygousWeinreichPaths src/environment.cpp src/modelFunctions.cpp src/randomv.cpp src/computeHeterozygousWeinreichPaths.cpp $(CFLAGS) -I boost_1_48_0/boost
clean: 
	rm -f bin/a.out bin/stability
