#include <iostream>
#include <vector>
#include <gsl/gsl_math.h>
#include "randomv.h"
#include "modelFunctions.h"
#include "environment.h"
// version = 4;

/*
g++ geometry.cpp randomv.cpp modelFunctions.cpp environment.cpp -o geometry.exe -lgsl -lgslcblas

usage:
./geometry.exe [DETAILED] [rinit] > geom.out &

 where 
 DETAILED   can be 0 or 1. Default is 0
            if 1 the output will be: coordinates (r, theta) - heterozygote relative fitnes - state - probability of not being extinct - relative fitness of homozygote
            if 0:                    Pbalanced - max. mutation step/distance of A from optimum

calculates the probability of balanced alleles

with uniform mutation function
#b: ./geometry.exe > data/uniformpbal.dat &
uniform with rinit = 1
./geometry.exe 1 2 > data/ugeomclose.dat &
./geometry.exe 1 7 > data/ugeomfar.dat

with exponential mutaton function and rinit = 2
#a: ./geometry.exe 1 2 > data/geom.dat &
#b: ./geometry.exe > data/pbal.dat &

pdf("manuscript/figures/geom.pdf")
#png("manuscript/figures/presentation/geom.png")
see figures.R
dev.off()

pdf("manuscript/figures/pbal.pdf")
#png("manuscript/figures/presentation/pbal.png")
library(ggplot2)
b<-read.table("data/pbal.dat")
qplot(b$V1,b$V2,xlab="max. mutation step / distance of A from optimum", ylab = expression(P[balanced]))
dev.off()

*/
using namespace std;
int main(int argc, char* argv[]){
  
  bool detailed = false; //default values
  double r      = 1;
  if (argc == 3){
    detailed = argv[1];
    r = atof(argv[2]);
  }
  bool uSampling = false; //sampling method form mutation function

  //input A(ra,rtheta) and Optimum (ro,thetao)
  double theta   = 0;
  double ro      = 0;
  double tho     = 0;
  int    repeat  = 100000;
  vector<double> a;
  a.push_back(r);
  a.push_back(theta);
  //instantiate classes
  randomv myR;
  randomv &ra = myR;
  modelFunctions myModel(ra);
  environment myEnvironment(ro, tho);
  environment &e = myEnvironment;
  myModel.setFitnessFunction('g');
  myModel.setA(1);
  myModel.setB(0);
  myModel.setC(1);
  myModel.setD(2);
  myModel.setMutationFunction('u');
  double maxStep = 5; //maximum size of mutation step for uniform sampling
  for(double c = 0; c<maxStep*r; c+=0.01){
    vector<double> parameters;  //new mutation function parameters
    if (detailed == true){
      c = maxStep*r;
      parameters.push_back(5*r); //Mutation function ~ Exp(1)
    }else{
      parameters.push_back(c);
    }
    myModel.setParameters(parameters);
    double sumBal = 0;
    double sumTot = 0;
    //loop repeat times
    for(int i = 0; i < repeat; i++){
      //  mutate A according to some function
      double mr = 0;
      double mth = 0;
      if(uSampling == true){
	//method 1 uniform
	double xrange = c;
	double yrange = c;
	double x =  (ra.sampleUniform() - 0.5) * 2 * xrange;
	double y =  (ra.sampleUniform() - 0.5) * 2 * yrange;
	// find r and theta
	mr = sqrt(y*y+x*x); //keep only r > 0
	mth = atan2(y,x); //[-pi,pi]
      }else{
	//method 2 angular
	mr = myModel.fR(ra);
	mth = myModel.fTheta(ra);
      }      
      vector<double> m; //mutator vector
      m.push_back(mr);
      m.push_back(mth);
      vector<double> b = myModel.add(a,m);
      //  find in which region is B (i, ii:balanced, or iii)
      //     find waa, wab, wbb
      double waa = e.fW(a,myModel);
      double wbb = e.fW(b,myModel);
      double wab = e.diploidFitness(a,b,myModel);
      double hs = wab/waa - 1;
      double s = wbb/waa - 1;
      int state = -1;
      if(hs>0){
	sumTot += hs;
      }
      if (wab>waa && wab>wbb){
	state = 4;
	sumBal += hs;
      }
      if (detailed == true){
	// 1 Underdominance
	// 2 Recesive
	// 3 Dominant
	// 4 Overdominance
	// 5 Codominant
	// 6 Complete dominance
	// 7 Complete recesiveness
	if (wab>waa && wab>wbb){
	  state = 4;
	  sumBal += hs;
	}else if (waa > wbb){
	  double h = (waa-wab)/(waa-wbb);
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
	}else if (wbb > waa){
	  double h = (wbb-wab)/(wbb-waa);
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

	}else if(waa == wbb){
	  state = 5; // waa == wbb codominance
	}
	//in boundary regions state = -1
	double pi = 2.*hs; //probability of new allele survival (assuming infinite N and that selection acts on the heterozygote)
	cout <<b.front()<<' '<<b.back()<<' '<<(hs>0?hs:0)<<' '<<state<<' '<<pi<<' '<<s<<' '<<mr<<endl;
      }
    }
    if(detailed ==false){
      cout <<c/r<<' '<<sumBal/sumTot<<endl;
    }
  }
  return EXIT_SUCCESS;
}


