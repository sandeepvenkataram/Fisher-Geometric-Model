CC=g++
CFLAGS=-lgsl -lgslcblas
DEPS = allele.h environment.h modelFunctions.h population.h randomv.h runningStat.h 

a:	tester.cpp allele.cpp environment.cpp modelFunctions.cpp population.cpp randomv.cpp runningStat.cpp stability.cpp
	$(CC) -O3 -o a.out tester.cpp allele.cpp environment.cpp modelFunctions.cpp population.cpp randomv.cpp runningStat.cpp stability.cpp $(CFLAGS)
geometry: geometry.cpp randomv.cpp modelFunctions.cpp environment.cpp
	$(CC) -O3 -o geometry.exe geometry.cpp randomv.cpp modelFunctions.cpp environment.cpp $(CFLAGS)
clean: 
	rm -f a.out balanced.exe