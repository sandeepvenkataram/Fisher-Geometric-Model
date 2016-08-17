#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

/*
	Defines the environment that the population is evolving on, e.g. the position of the phenotypic optimum
	The code to allow for periodic fluctuation in the optimum's position has not been extensively tested
*/

#include <vector>
#include "randomv.h"
#include "modelFunctions.h"

using namespace std;
class environment{
public:
  environment(void);
  environment(vector<double> opt);
  void setOptimum( vector<double> x ){ optimum = x; }
  void setPeriodic( bool x ){ periodic = x; }
  vector<double> getOptimum( void ){ return optimum; }
  bool getPeriodic( void ){ return periodic; }
  int getTimeOfChange(void){return timeOfChange; }
  void change(randomv &r, modelFunctions &myModelRef, int t);
  
  double fW(vector<double> rv, modelFunctions &myModelRef);
  double diploidFitness(vector<double> a, vector<double> b, modelFunctions &m); //find diploidFitness from coordinates of each allele
  ///  double diploidFitness(vector<double> a, vector<double> b, randomv &r, modelFunctions &m); //find diploidFitness from
  void setUsePredefinedMutations(bool x){usePredefinedMutations = x;}
  bool getUsePredefinedMutations(void){return usePredefinedMutations;}
  
  void setPredefinedMutations(vector< vector<double> > &x){ predefinedMutations = x;}
  vector< vector<double> > getPredefinedMutations(void){ return predefinedMutations;}
  //for Harmonic
  void   setF(double x){ f = x; }
  double getF(void){ return f; }
  void   setPhi(double x){ phi = x; }
  double getPhi(void){ return phi; }
  void   setFixedTheta(double x){ fixedTheta = x; }
  double getFixedTheta(void){ return fixedTheta; }
  
private:
  //for the harmonic motion of the mean
  double f;               //! Frequency
  double phi;             //! Phase
  double fixedTheta;      //! Fixed angle
  vector<double> optimum; //! Coordinates for the optimum
  bool periodic;          //! True if the environmental change is periodic
  int timeOfChange;       //! Most recent time the optimum changed
  bool usePredefinedMutations;
  vector< vector<double> > predefinedMutations;
  
};
#endif
