#ifndef STABILITY_H
#define STABILITY_H
#include <gsl/gsl_matrix.h>
#include <vector>
using namespace std;

class stability{
 public:
  stability(vector<double> ws, int n);
  int stable(void);
 private:
  int negative_definite(gsl_matrix *M, int n);
  int delta_condition(gsl_matrix *M, int n);
  gsl_matrix *A;
  int n; //! dimensions
};
#endif
