This is an implementation of forward Wright-Fisher simulations under Fisher's geometric model.

It uses code modified from Sellis et al 2011 (PNAS) to allow for arbitrary numbers of dimensions. This code is published in conjunction with the publication Venkataram, Sellis and Petrov (in preparation).

Source code is hosted in the src/ folder, and can be compiled if the generated binaries in the bin/ folder are not compatible with your system

Example input files are hosted in the example/ folder, which also describe the various parameter options and the format of the files


<b>Build Instructions</b>

make - makes the program a.out, which contains the simulation software. Requires the gcc compiler

make clean - removes all built files
