#ifndef ALLELE_H
#define ALLELE_H
#include <iostream>
#include <vector>
#include "modelFunctions.h"
#include "environment.h"

//! allele class
/*! It represents a single locus haplotype
 */
class allele{
 public:
  allele(void); //! default constructor
  allele(double a, vector<double> b, modelFunctions &myModelRef); //! constructor with parameters
  void           setP(double a){ if(a>=0){p = a;}else{cerr<<"Warning: negative frequency!";exit(0);} }
  void           setRv(vector<double> &a){ rv = a; }
  void           setId(int a){ id = a; }
  void           setIdp(int a){ idp = a; }
  void           setBd(int a){ bd = a; }
  /* unused
   *  void           setBalanced(vector<int> a){ balanced = a; }
   *  void           incrBalanced(int a);
   */
  void           setFixed(bool a){ fixed = a; }
  void           setPhe(double x){ phe = x; }
  void           setDistance(double x){ distance = x; }
  void		 setfW(double x){ fW = x; }
  double         getP(void){ return p; }
  vector<double> getRv(void){ return rv; }
  double         getW(modelFunctions &myModelRef, environment &envRef);
  int            getId(void){ return id; }
  int            getIdp(void){ return idp; }
  int            getBd(void){ return bd; }
  //  vector<int>    getBalanced(void){ return balanced; }
  bool           getFixed(void){ return fixed; }
  double         getPhe(void){ return phe; }
  double         getDistance(void){ return distance; }
  double         getfW(void){ return fW; }
  vector<double> getCartesian(modelFunctions &myModelRef); //!return cartesian coordinates in vector
  void addSNP(int x){ SNPs.push_back(x); }
  int length(void) {return SNPs.size(); } 
  void setSNPs( vector<int> x){ SNPs = x; }
  vector<int> getSNPs(void){ return SNPs; }

 private:
  double         p;        //! Frequency
  vector<double> rv;       //! Vector of polar coordinates (r,theta1, theta2, theta3, theta4 etc.)
  int            id;       //! Identification number for tracking alleles and debuging
  int            idp;      //! Id of parental allele
  int            bd;       //! Allele's birthday (generation that it appeared)
  //  vector<int>    balanced; /// Number of generations the allele is under balancing selection
                           /** multiple elements represent generations where the allele
			    is balancing in two or more equilibria*/
  bool           fixed;    //! Whas the allele ever fixed? (ancestral)
  double         phe;      //! Phenotypic effect (Wmutated-Wmutant)
  double         distance; //! Eucledian distance between mutated and mutant
  double	fW;	   //! fitness of homozygouts allele
  vector<int> SNPs; 
};
#endif
