#include <cmath>
#include <vector>
#include "modelFunctions.h"
#include "environment.h"
#include "allele.h"

//! default constructor initialize everything
allele::allele(){
  p        = 0.;
  rv       = vector<double> (1,0);
  id       = 0;
  idp      = -1; //unknown ancestry
  bd       = 0;
  fixed    = false;
  phe      = 0;
  distance = 0;
}

//! constructor with coordinates
/*! initialize everything and set fixed true only if its frequency is 1
*/
allele::allele(double a, vector <double> b, modelFunctions &myModelRef){
  p        = a;
  rv       = b;
  id       = 0;
  idp      = -1;
  bd       = 0;
  phe      = 0;
  distance = 0;
  if (p == 1){
    fixed = true;
  }else{
    fixed = false;
  }

}

double allele::getW(modelFunctions &myModelRef, environment &envRef){
  double w = envRef.fW(rv, myModelRef);
  return w;
}

/* unused
void allele::incrBalanced(int a){
  if (balanced.size() < a+1){
    balanced.push_back(0);
  }
  balanced.at(a)++;
}
*/

vector<double> allele::getCartesian(modelFunctions &myModelRef){
 return myModelRef.polarToCartesian(rv);
}
