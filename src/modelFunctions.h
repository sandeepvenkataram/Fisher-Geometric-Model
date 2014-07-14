#ifndef MODELFUNCTIONS_H
#define MODELFUNCTIONS_H

#include <vector>
#include "randomv.h"

//! Class containing the model functions for organisms and the environment
//! version = 1;
using namespace std;
class modelFunctions {
 public:
  modelFunctions(randomv &r); //! constructor
  double fR(randomv myR);     
  double fTheta(randomv myR);
  vector<double> add(vector<double> a, vector<double> b);
  double vdistance(vector<double> a, vector<double> b, int print);
  double tau(double p, double N);
  //fitness function
  char   getFitnessFunction(void)  {return fitnessFunction;}
  void   setFitnessFunction(char x){fitnessFunction = x;}
  // for linear
  double getMaxR(void)    { return maxR;}
  void   setMaxR(double x){ maxR = x; }
  // for Gaussian
  double getA(void)    { return a; }
  void   setA(double x){ a = x; }
  double getB(void)    { return b; }
  void   setB(double x){ b = x; }
  double getC(void)    { return c; }
  void   setC(double x){ c = x; }
  double getD(void)    { return d; }
  void   setD(double x){ d = x; }
  //parameters for new mutation functions
  vector<double> getParameters(void)   { return parameters; }
  void   setParameters(vector<double> x) { parameters = x; }
  void   setMutationFunction(char a)     { mutationFunction = a; }
  char   getMutationFunction(void)       { return mutationFunction; }
  //environmental
  void   setWalk(  bool x )   { walk = x; }
  void   setSigma( double x ) { sigma = x; }
  void   setAlpha( double x ) { alpha = x; }
  void   setBeta(  double x ) { beta = x; }
  bool   getWalk(  void )     { return walk; }
  double getSigma( void )     { return sigma; }
  double getAlpha( void )     { return alpha; }
  double getBeta(  void )     { return beta; }
  double getStep(randomv myR);
  int    getStates(double waa, double wbb, double wab);
  int getNumDimensions(void){return numDimensions;}
  void setNumDimensions(int x){numDimensions=x;}
  vector<double> polarToCartesian(vector<double> polar);
  vector<double> cartesianToPolar(vector<double> cart);
  vector<double> getMutationVector(randomv myR, int numDimensions);
  //global printing
  //  bool printTable;   //table of balanced pairs
 private:
  char   fitnessFunction; //! linear or Gaussian fitness function
  double maxR;            //! maximum distance from optimum
  //@{ 
  /**Parameters of Gaussian fitness function
     \f$f(x)=a*e^{(\frac-(x-b)^d/(2*c^2))}\f$
   */
  double b; //! mean
  double a; //! maximum
  double c; //! shape
  double d; //! shape of curve
  //}@
  //for new mutation functions
  vector<double> parameters;
  char mutationFunction;
  //environment
  bool walk; // false: Wiener process (Brownian motion)
             // true: Levy flight
  double sigma; // normal step for Wiener 
  double alpha; // pareto step for Levy flight
  double beta;  //  >>
  int numDimensions;
  
};
#endif
