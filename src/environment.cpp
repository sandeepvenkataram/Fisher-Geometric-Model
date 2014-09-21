#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_math.h>
#include "randomv.h"
#include "modelFunctions.h"
#include "environment.h"

environment::environment(){
  optimum.push_back(0);
  optimum.push_back(0);
  timeOfChange = 0;
}

environment::environment(vector<double> opt){
  optimum = opt;
  timeOfChange = 0;
}

void environment::change(randomv &rnd, modelFunctions &myModelRef, int t){ //this is broken with the current implementation as fitness is only calculated when the new allele is created.
  return;
  double r;
  double theta;
  vector<double> move;
  if(getPeriodic() == true){
    /*
      Simple harmonic motion
      f:   frequency
      phi: phase
      d:   displacement
    */
    double d = 1; //in order b:[0,2]
    theta = fixedTheta; //move in a straight line
    double td  = (double) t;
    r = sin(2*M_PI*f*td + phi) + d;
    move.push_back(r);
    move.push_back(theta);
    optimum = move;    
  }else{
    r = myModelRef.getStep(rnd);
    theta = myModelRef.fTheta(rnd);
    move.push_back(r);
    move.push_back(theta);
    optimum = myModelRef.add(move, optimum);    
  }
  timeOfChange = t;
}

double environment::fW(vector<double> rv, modelFunctions &myModelRef){
  /*
    'Gaussian function'
    f(x)=a*exp(-(x-b)^d/(2*c^2))
    b: mean
    a: maximum
    c: shape
    b = xd
    a = 1
    c = 1
    d = 2
  */
  double w = 0;
  //find distance from optimum
 // cerr<<"env.fW";
  for(int i=0; i<rv.size();i++){
//	cerr<<" "<<rv[i];
  }

  double x = myModelRef.vdistance(rv, optimum,0);
  //apply fitness function
  if (myModelRef.getFitnessFunction() == 'g'){
    //gaussian
    w = myModelRef.getA()*exp((-1*pow(x,myModelRef.getD()))/(2*pow(myModelRef.getC(),2)));
  }else if(myModelRef.getFitnessFunction() == 'q'){ // quadratic
    w = 1-pow(x,2)/myModelRef.getA();
    if(w<0){
      w=0;
    }
  }else if (myModelRef.getFitnessFunction() == 'l'){
    //linear
    if(x>=myModelRef.getMaxR()){
      w = 0;
    }else{
      w = myModelRef.getA()*(1-(x/myModelRef.getMaxR()));
    }
  }else{
    cout<<"unknown fitness function: "<<myModelRef.getFitnessFunction() <<endl;
    exit(1);
    
  }
// cerr<<" "<<w<<endl;
  return w;
}

double environment::diploidFitness(vector<double> a, vector<double> b, modelFunctions &m){
  //add vectors
  vector<double> newRv = m.add(a, b);
  //divide by 2 as we are assuming phenotypes are completely co-dominant
  newRv.front()/=2.;
  //find fitness (distance from optimum is computed within fW)
  double wab = fW(newRv, m);
  return wab;
}
