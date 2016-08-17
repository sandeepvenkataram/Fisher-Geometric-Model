#include <vector>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include "stability.h"
using namespace std;
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
