# Evolutionary predictability with historical reconstruction

The repository hosts an implementation of forward Wright-Fisher simulations under Fisher's geometric model (Fisher 1931).

It uses code modified from Sellis et al 2011 (PNAS) to allow for arbitrary numbers of dimensions. This code is published in conjunction with the publication Venkataram, Sellis and Petrov (in preparation).

Source code is hosted in the `src/` folder, and can be compiled if the generated binaries in the `bin/` folder are not compatible with your system.

Example input files are hosted in the `example/` folder, which also describe the various parameter options and the format of the files.


*Build Instructions*

Extract the archive and change directory into the root directory. This software requires the gcc C++ compiler, along with the GSL library

`make` - makes the program FGM, which contains the simulation software. Requires the gcc compiler

`make stability` - makes the program stability, which calculates whether a set of alleles can be maintained as a stable polymorphism using the method of Kimura 1969

`make clean` - removes all built files

*Execution Instructions*

**FGM**

execution of the program:

./bin/FGM initPopulation.dat initEnvironment.dat par.dat output

output 

	.table	one line per allele for each generation, gives allele frequencies
	
	.ts	   one line per genaration, gives mean fitness per generation etc
	
	.edges	mutation events parent,offspring and time
	
	.b		summary statistics for each run
	

**Stability**

This program computes whether a set of alleles can be maintained as a stable polymorphism according to the method of Kimura et al 1969

Command line arguments for input: 
  
d - the number of alleles to be tested, an integer value
	
A further d*(d+1)/2 arguments defining the fitness of all possible genotypes for the d alleles. 
	
The order of these fitness value arguments should essentially be a row-major representation of the upper-right half of the dxd fitness matrx
	
e.g. if there are 3 alleles A, B and C, the fitness values should be given in the order of AA, AB, AC, BB, BC, CC

