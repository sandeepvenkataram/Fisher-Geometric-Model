#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include "randomv.h"
#include "modelFunctions.h"

using namespace std;

modelFunctions::modelFunctions(randomv &r){
  //for Gaussian
  b = 1; //mean
  c = 1; //shape
  a = 1; //maximum for s
  d = 2;
}

double modelFunctions::fR(randomv myR){ // sample random magnitude for mutation vector based on parameters
  double newR;
  switch (mutationFunction){
  case 'g': //GAMMA
    { 
      double shape = parameters.front();
      newR = myR.sampleGamma(shape,1);
    }
    break;
  case 'u': //UNIFORM
    newR = myR.sampleUniform()*parameters.front();
    break;
  case 'e': //EXPONENTIAL
    newR = myR.sampleExponential(parameters.front());
    break;
  case 'n': //NORMAL
    {
      newR = myR.sampleNormal(parameters.at(0),parameters.at(1));
      if (newR<0){
	newR = -newR;
      }
    }
    break;
  default:
    cout<<"unknown function: "<<mutationFunction<<endl;
    exit(1);
  }
  return newR;
}

double modelFunctions::fTheta(randomv myR){ // sample random direction in radians based on uniform distribution
  double newTheta = myR.sampleUniform();
  return newTheta*2*M_PI - M_PI; // (-PI,PI]
}

vector<double> modelFunctions::getMutationVector(randomv myR, int numDimensions){ // sample mutation vector
	double fr = fR(myR);
	double direction[numDimensions];
	myR.sampleDirVector(numDimensions, direction);
	std::vector<double> result (direction, direction+sizeof(direction)/sizeof(double));
	result = cartesianToPolar(result);
	result[0] =fr;	
	return result;	
}

vector<double> modelFunctions::cartesianToPolar(vector<double> cart){ // actually converting it to hyperspherical coordinates for more than 2 dimensions, but naming due to legacy code
	double r = 0;
	int numDim = cart.size();
	vector<double> result;
	for(int i=0; i<cart.size();i++){ // initialize vector, and compute vector magnitude
		r+=cart[i]*cart[i];
		result.push_back(0);
	}
	r = pow(r,0.5);
	result[0]=r;
	if(r==0){ // if vector has 0 magnitude, set all angle values to 0
		for(int i=1; i<cart.size();i++){
			result[i]=0;
		}
		return result;
	}
	for(int i=1; i<cart.size();i++){ 
		double val = cart[i]*cart[i];
		for(int j=1;j<=i;j++){
			result[j]+=val;
		}
	}
	for(int j=1;j<cart.size()-1;j++){ // for all but the last angle
		if(result[j]==0){
			result[j]=0;
		}else{
			result[j]=M_PI/2-atan(cart[j-1]/pow(result[j],0.5));
		}
	}
	if(cart[numDim-1]==0){ // if the last dimension is 0, the last angle is 0
		result[numDim-1]=0;
	}else{ // otherwise, it is this mess
		result[numDim-1]=M_PI-2*atan((cart[numDim-2]+pow(cart[numDim-1]*cart[numDim-1]+cart[numDim-2]*cart[numDim-2],0.5))/cart[numDim-1]);
	}
	return result;
	
}

vector<double> modelFunctions::polarToCartesian(vector<double> polar){ // convert hyperspherical coordinates into cartesian
	double r = polar.front();
	
	vector<double> Cartesian;
	for(int i=1;i<=polar.size();i++){
		double res = r;
		for(int j=1;j<=i;j++){
			if(j==i&&j<polar.size()){
				res=res*cos(polar[j]);
			}else if(j<polar.size()){
				res=res*sin(polar[j]);
			}else{
				
			}
		}
		Cartesian.push_back(res);
	}
	return Cartesian;
}

vector<double> modelFunctions::add(vector<double> a, vector<double> b){
  //variables with convenient names
  vector<double> cartA = polarToCartesian(a);
  vector<double> cartB = polarToCartesian(b);
  vector<double> sum;
  for(int i=0; i<cartA.size();i++){
	sum.push_back(cartA[i]+cartB[i]);
  }
 // for(int i=0; i<cartA.size();i++){
//	cout<<"add: "<<a[i]<<" "<<cartA[i]<<" "<<b[i]<<" "<<cartB[i]<<" "<<sum[i]<<endl;
 // }
  return cartesianToPolar(sum);
}

double modelFunctions::tau(double p, double N){
  double t = 0;
  t = 4*N*(-p*log(p)/(1-p));
  return t;
}

double modelFunctions::getStep(randomv myR){
  double step;
  if (walk == false){
    step = myR.sampleNormal(0,sigma);
  }else{
    step = myR.samplePareto(alpha, beta);
  }
  return step;
}

vector<double> modelFunctions::getMeanPhenotype(vector< vector<double> > alleles, vector<double> freqs){
	vector<double> resultPhenotype;
	if(alleles.size()!=freqs.size()){
		return resultPhenotype;
	}
	for(int i=0; i<getNumDimensions();i++){ // initialize vector, and compute vector magnitude
		resultPhenotype.push_back(0);
	}
	
	for(int i=0; i<alleles.size(); i++){
		vector<double> cartVec = polarToCartesian(alleles[i]);
		for(int j=0; j<cartVec.size(); j++){
			resultPhenotype[j]+=cartVec[j]*freqs[i];
		}
	}
	return cartesianToPolar(resultPhenotype);
}

double modelFunctions::vdistance(vector<double> a, vector<double> b, int print){ // distance between two polar coordinates
  vector<double> cartA = polarToCartesian(a);
  vector<double> cartB = polarToCartesian(b);
  double result;
  if(print!=0){
	cout<<"Starting vdistance computation!"<<endl;
  }
  for(int i=0; i<cartA.size();i++){
	  if(print!=0){
	  cout<<"vdistance: "<<a[i]<<" "<<cartA[i]<<" "<<b[i]<<" "<<cartB[i]<<endl;
	  }
	result+=pow(cartA[i]-cartB[i],2);
  }
  if(print!=0){
  cout<<"End Distance is:\t"<<pow(result,0.5)<<endl;
  }
  return pow(result,0.5);
  
}

int modelFunctions::getStates(double wa, double wb, double wab){
  int state = -1;
  /*
    1 Underdominance
    2 Recesive
    3 Dominant
    4 Overdominance
    5 Codominant
    6 Complete dominance
    7 Complete recesiveness
  */
  if (wab>wa && wab>wb){
    state = 4;
  }else if (wa > wb){
    double h = (wa-wab)/(wa-wb);
    if(h>1){
      state = 1;
    }else if(h>0.5 && h <1){
      state = 3;
    }else if (h == 0.5){
      state = 5;
    }else if(h>0 && h < 0.5){
      state = 2;
    }else{
      state = 6;
    }
  }else if (wb > wa){
    double h = (wb-wab)/(wb-wa);
    if(h>1){
      state = 1;
    }else if(h>0.5 && h <1){
      state = 2;
    }else if (h == 0.5){
      state = 5;
    }else if(h>0 && h < 0.5){
      state = 3;
    }else{
      state = 7; //Complete recesiveness
    }
    
  }else{
    state = 5; // waa == wbb codominance
  }
  return state;
}
