#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

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
  ///  double diploidFitness(vector<double> a, vector<double> b, randomv &r, modelFunctions &m); //find diploidFitness from
  double diploidFitness(vector<double> a, vector<double> b, modelFunctions &m); //find diploidFitness from coordinates

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
};
#endif
