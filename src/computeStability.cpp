#include <vector>
#include <iostream>
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
	
	if(d==1){
		cout<<1<<endl;
		exit(0);
	}
	
	vector<double> vec;
	
	for(int i= 2; i<argc; i++){
		vec.push_back(atof(argv[i]));
	}
	
	stability stab(vec,d);
	int stable = stab.stable();
	cout<<stable<<endl;

}