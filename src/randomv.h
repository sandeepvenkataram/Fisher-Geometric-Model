#ifndef RANDOMV_H
#define RANDOMV_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//! random variable class
/*! an object oriented interface for sampling from distributions using gsl
 */
class randomv{
 public:
  randomv(void);                                  //! default constructor
  double sampleUniform(void);                     //! sample from the uniform distribution [0,1)
  int sampleUniformInt(int n);                    //! sample from the uniform distribution [0,n-1)
  unsigned int sampleBinomial(double p, int N);   //! sample from the binomial distribution
  double sampleNormal(double mu, double sigma);   //! sample from a normal distribution with mean mu and standard deviation sigma
  double sampleNormal(double sigma);              //! sample from a normal distribution with mean 0 and standard deviation sigma
  double sampleExponential(double mu);            //! sample from an exponential distribution with mean mu
  double sampleGamma(double alpha, double beta);  //! sample from a gamma distribution with parameters alpha and beta
  void sampleMultinomial(size_t K, unsigned int N, const double p[], unsigned int n[]); //sample from a multinomial distribution
  double samplePareto(double alpha, double beta); //! sample from a pareto distribution with parameters alpha and beta \fp(x) dx = (a/b) / (x/b)^{a+1} dx

  double samplePoisson(double lambda);
  void sampleDirVector(size_t n, double * x);
 private:
  const gsl_rng_type *Gt; 
  gsl_rng *Gr;
};
#endif
