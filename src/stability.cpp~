#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include "stability.h"
using namespace std;

stability::stability(vector<double> ws, int d){
  n = d;
  A = gsl_matrix_alloc(n,n);
  int counter = 0;
  for (int i = 0; i < n; i++){
    for (int j = i; j < n; j++){
      gsl_matrix_set(A,i,j, ws.at(counter));
      if(i != j){
	gsl_matrix_set(A,j,i, ws.at(counter)); //symmetric	
      }
      counter++;
    }
  }  
}

int stability::negative_definite(gsl_matrix *M, int n)
{
  int neg_def=1;

  gsl_vector *eval = gsl_vector_alloc(n);     
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(n);       
  gsl_eigen_symm(M,eval,w);     
  gsl_eigen_symm_free(w);  
               
  for (int i=0; i<n; i++)
    {
      double eval_i  = gsl_vector_get(eval,i);   
      if(eval_i>=0) { neg_def=0; }
    }

  gsl_vector_free(eval);
 
  return neg_def;
}

int stability::delta_condition(gsl_matrix *M, int n)
{
  int condition=1;

  int s;

  for (int r=0; r<n; r++)
    {
      // calculate determinant delta_r from M by substituting all elements of the r-th column by 1 

      gsl_matrix *M_r = gsl_matrix_alloc(n,n);

      for(int i=0; i<n; i++) 
	{ 
	  for(int j=0; j<n; j++)
	    {
	      if(j!=r) { gsl_matrix_set(M_r,i,j,gsl_matrix_get(M,i,j)); }
	      else     { gsl_matrix_set(M_r,i,j,1.0); }
	    }
	}

      // calculate sign of determinant of M_r via LU decomposition

      gsl_permutation *p = gsl_permutation_alloc(n);     
      gsl_linalg_LU_decomp(M_r,p,&s);
      double sgn = gsl_linalg_LU_sgndet(M_r,s);

      // check condition

      if(pow(-1,n-1)*sgn<=0) { condition=0; }

      gsl_permutation_free(p);
      gsl_matrix_free(M_r);
    }

  return condition;
}

int stability::stable(void){
  // quadratic form T from Eq. (2.4)

  gsl_matrix *T = gsl_matrix_alloc(n,n);
  
  for(int i=0; i<n; i++)
    {
      for(int j=0; j<n; j++)
	{
	  double element =  gsl_matrix_get(A,i,j)-gsl_matrix_get(A,i,n-1)-gsl_matrix_get(A,j,n-1)-gsl_matrix_get(A,n-1,n-1);

	  gsl_matrix_set(T,i,j,element);
	}
    }
  
  // check conditions
  if ((negative_definite(T,n) == 1) &&
      (delta_condition(A,n) == 1)){
    return 1;
  }else{
    return 0;
  }
  return -1;
}
