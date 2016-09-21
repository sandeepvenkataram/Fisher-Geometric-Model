#include <vector>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include "stability.h"
using namespace std;

/*
	Computes whether a set of alleles can be maintained as a stable polymorphism according to the method of Kimura et al 1969
	Command line arguments for input: 
	d - the number of alleles to be tested, an integer value
	a further (d*(d+1)/2) arguments defining the fitness of all possible genotypes for the d alleles. 
	
	the order of these fitness value arguments should essentially be a row-major representation of the upper-right half of the dxd fitness matrx
	e.g. if there are 3 alleles A, B and C, the fitness values should be given in the order of AA, AB, AC, BB, BC, CC
*/

int main(int argc, char* argv[]){
	if(argc < 3){
		cerr<<"Insufficient arguments!";
		exit(1);
	}
	int d = atoi(argv[1]);
	
	if((argc - 2) != (d*(d+1)/2)){
		cerr<<"Insufficient values for "<<d<<" alleles!";
		exit(1);
	}
	
	/*if(d==1){
		cout<<1<<endl;
		exit(0);
	}*/
	
	vector<double> vec;
	
	for(int i= 2; i<argc; i++){
		vec.push_back(atof(argv[i]));
	}
	
	stability stab(vec,d);
	int stable = stab.stable();
	cout<<stable<<"\t"<<stab.getMeanFitness()<<"\t";
	gsl_matrix *eqFreqs = stab.getEqFreqs();
	for(int i=0; i<d; i++){
		cout<<gsl_matrix_get(eqFreqs,i,0)<<"\t";
	}
	cout<<endl;

}
