#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <gsl/gsl_sort.h>
#include <iterator>
#include "randomv.h"
#include "modelFunctions.h"
#include "environment.h"
#include "allele.h"
#include "population.h"
#include "stability.h"

using namespace std;
population::population(randomv &r, modelFunctions &myModelRef, vector<allele> &init, int pl, environment &envRef, ofstream &fsr){
  simTime         = 0;
  alleles         = init;
  ploidy          = pl;
  wm              = updatewm(myModelRef, envRef); //initialize mean fitness
  alleleCounter   = alleles.size();
  wildW           = 0;
  maxWij          = 0;
  //  balancedPairs   = 0;
  phm.push_back(0);
  phm.push_back(0);
  for(int i=0; i<myModelRef.getNumDimensions();i++){
	phenotypeVar.push_back(0);
  }
 // cout<<"pop1"<<endl;
  updatePhm(myModelRef);
 // cout<<"pop2"<<endl;
  areBalanced(myModelRef, envRef);
 // cout<<"pop3"<<endl;
}

double population::FindNeutralSphere(modelFunctions &myModelRef){
  r0 = myModelRef.getC()*sqrt(2./(double)N);
  return r0;
}

void population::evolve(randomv &r, modelFunctions &myModelRef, environment &envRef, ofstream &fse){
  /*  if(ploidy == 2){
    balancedPairs = 0;
    }*/
  //mutate
 // cout<<"evolve1"<<endl;
  mutate(r, myModelRef, envRef, fse);
//  cout<<"evolve2"<<endl;
  propagate(r, myModelRef, envRef);
 // cout<<"evolve3"<<endl;
  wm = updatewm(myModelRef, envRef);
 // cout<<"evolve4"<<endl;
  if (alleles.size() == 1){
    alleles.front().setFixed(true);
  }
  simTime++;
}

void population::printStatus(modelFunctions &myModelRef, ofstream &fp_out, ofstream &fts, environment &envRef, ofstream &fsr){
  /* unused
   *  //make a histogram-vector for currently balanced alleles
   *  vector<int> h; 
   */
  int counter = 0;
  updatePhm(myModelRef);
  for (vector<allele>::iterator it = alleles.begin();
       it != alleles.end(); ++it){
    /* unused
     *    if(ploidy == 2){
     *      int multiBalanced = curBalanced.at(counter);
     *      if( multiBalanced >= h.size()){
     *	h.resize(multiBalanced + 1,0);
     *      }else{
     *	h.at(multiBalanced)++;
     *      }
     *    }
     */
    double wi = (*it).getfW();
    //double si;
    //    if (wi < wildW){
    //      si = 1-wi/wildW;
    //    }else if (abs(wi - wildW) < numeric_limits<double>::min()){
    //      si = 1;
    //    }else{
    //      cerr<<"error 4"<<endl;
    //      exit(1);
    //    }
    double si = 0;
    if (abs(wi - maxGenW) < numeric_limits<double>::min()){
      si = 1;
    }else{
      si = 1-wi/maxGenW;
    }
    vector<double> rvi = (*it).getRv();
    double pi = (*it).getP();
    int idi = (*it).getId();
    int age = simTime - (*it).getBd();
    fp_out << simTime <<' ';         // 1. time (generations)
    fp_out <<pi<< ' ';               // 2. frequency (p)
    for(int i=0; i<rvi.size(); i++){
	fp_out <<rvi[i]<<' ';       // 3. polar coordinates r
    }
    fp_out <<idi<<' ';               // 5. allele id
    fp_out <<wm<<' ';                // 6. mean fitness (wm)
    fp_out <<wi<<' ';                // 6. fitness (w)
    fp_out <<age<<' ';               // 7. allele age
    fp_out <<si<<' ';                // 8. selection coefficient (s)
    fp_out <<(*it).getFixed()<<' ';  // 9. is the allele ancestral (ever been fixed)
    //    if (ploidy == 2 && curBalanced.at(counter) != 0 ){ //10. if balanced
    //      fp_out << 1<<' ';
      /*
       *      vector<int> balancedi = (*it).getBalanced();
       *      for(vector<int>::iterator jt = balancedi.begin();
       *	  jt != balancedi.end(); ++jt){
       *	fp_out << *jt/age<<' ';
       *	// 11,12... single,double...balanced ratio
       *      }
       */
    //    }else{
    //      fp_out << -1<<' '; // not balanced or haploid
    //    }
    fp_out <<(*it).getPhe()<<' ';         //10. phe
    fp_out <<(*it).getDistance()<<' ';    //11. distance
    fp_out <<(*it).getIdp()<<' ';         //12. parental allele
    fp_out<<(*it).length();

    fp_out <<endl;      //new line
    counter++;
  }
  
  fts <<simTime<<' ';                           // 1. time (generations)
  fts <<wm<<' ';                                // 2. mean fitness
  fts <<wildW<<' ';                             // 3. maximum allele's fitness
  fts <<alleles.size()<<' ';                    // 4. number of different alleles
  for(int i=0; i<envRef.getOptimum().size(); i++){
	fts <<envRef.getOptimum()[i]<<' ';       // 3. polar coordinates r
  }
  fts <<var<<' ';                               // 7. fitness variance
  for(int i=0; i<phm.size(); i++){
	fts <<phm[i]<<' ';       // 3. polar coordinates r
  }
  for(int i=0; i<phenotypeVar.size(); i++){
	fts <<phenotypeVar[i]<<' ';       // 3. polar coordinates r
  }
  fts <<maxGenW<<' ';                           //12. maximum genotype fitness
  fts <<simTime - envRef.getTimeOfChange()<<' '; //13. time since last change in the environment's optimum
  if (ploidy == 2){
    fts <<polymorphisms<<' ';                   //14. number of polymorphic alleles (with pi > d && pi < d)
    fts <<stableEquilibrium<< ' ';              //15. stable equilibrium for all alleles (with pi > d && taking 1/N into account)
    // compare maximum fitness difference between homozygotes - heterozygotes and 1/N
    double maxWdiff = 0; //maximum difference between wij and max(wii,wjj)
    //for each genotype calculate hetA and hetB and its frequency
    double sumA = 0;
    double sumB = 0;
    for (vector<allele>::iterator it = alleles.begin();	 it !=alleles.end(); ++it){
	//fts<<envRef.fW((*it).getRv(), myModelRef)<<" "<<(*it).getRv()[0]<<" "<<myModelRef.vdistance((*it).getRv(), envRef.getOptimum(),1)<<" ";
      for (vector<allele>::iterator jt = it; jt !=alleles.end(); ++jt){
	int iid = (*it).getId();
	vector<double> rvi = (*it).getRv();
	//(*it).setfW(envRef.fW(rvi,myModelRef)); //!recalculate the fitness every generation in case the optimum moved
	double wi = (*it).getfW();
	double pi = (*it).getP();
 	int jid = (*jt).getId();
	vector<double> rvj = (*jt).getRv();
	//(*jt).setfW(envRef.fW(rvj,myModelRef)); //!recalculate the fitness every generation in case the optimum moved
	double wj = (*jt).getfW();
	double pj = (*jt).getP();
	double wij = envRef.diploidFitness(rvi, rvj, myModelRef);
	double pij = pi*pj;
	double hetA = (wij - (wi+wj)/2)/wij;
	double diff = wij - max(wi,wj);
	double hetB = diff/wij;
	if (iid != jid){
	  sumA += 2 * pij * hetA;
	  sumB += 2 * pij * hetB;
	}else{
	  sumA += pij * hetA;
	  sumB += pij * hetB;
	}
	if(maxWdiff < diff){
	  maxWdiff = diff/maxGenW;
	}
      }
    }
    //Normalized average heterosis in current generation
    fts<<sumA<<' ';      //16. Heterosis A (wab - (waa+wbb)/2)/wab
    fts<<sumB<<' ';      //17. Heterosis B (wab - max(waa,wbb))/wab
    fts<<maxWdiff;       //18. maximum difference
  }
  fts <<endl;
}

int population::printStates(modelFunctions &myModelRef, environment &envRef, ofstream &fsr){
  for (vector<allele>::iterator it = alleles.begin();
       it != alleles.end(); ++it){
    int idi = (*it).getId();
    double pi = (*it).getP();
    vector<double> rvi = (*it).getRv();
    double wi = (*it).getfW();
    for (vector<allele>::iterator jt = it+1;
	 jt != alleles.end(); ++jt){
      int idj = (*jt).getId();
      double pj = (*jt).getP();
      vector<double> rvj = (*jt).getRv();
      double wj = (*jt).getfW();
      double wij = envRef.diploidFitness(rvi, rvj, myModelRef);
      int state = myModelRef.getStates(wi, wj, wij);
    
      fsr<<simTime<<' '<<idi<<' '<<idj<<' '<<state<<' '<<wi<<' '<<wj<<' '<<wij<<endl;
    }
  }
  return 1;
}
/*
int population::getSample(randomv &r, modelFunctions
			  &myModelRef, unsigned int sampleSize, ofstream &fsb, ofstream &fsr){

  double threshold = 0.1; //lower frequencies will be ignored
                          //should be moved to parameter file
  int alleleNumber = alleles.size();
  if (alleleNumber == 1){ //if allele fixed no sampling necessary
    return 0;
  }
  int counter = 0;
  //make arrays and initialize
  vector<double> p;
  vector<double> age;
  vector<int> bArray;
  vector<int> ids;
  double sumFreq = 0 ; //sum of the frequencies of the alleles to be sampled
  bool flag = false; //no balanced found
  
  for (vector<allele>::iterator it = alleles.begin();
       it != alleles.end(); ++it){
    double pi = (*it).getP();
    if (pi > threshold){ //ignore low frequencies
      p.push_back(pi);
      age.push_back(simTime - (*it).getBd());
      ids.push_back((*it).getId());
      sumFreq += pi;
      if (curBalanced[counter] != 0){
	bArray.push_back((*it).getBalanced().at(0)); //only record single not multiple 
	flag = true; //at least one balanced allele
      }else{
	bArray.push_back(-1);
      }
    }
    counter++;
  }
  if (flag == false){// if no balanced alleles
    //    return 0; keep track of how many generations
  }

  //normalize frequencies to add up to 1 (for the multinomial sampling)
  double pArray[p.size()];
  counter = 0;
  for (vector<double>::iterator it = p.begin();
       it !=p.end(); ++it){
    pArray[counter] = *it / sumFreq;
    counter++;
  }
  //make sample array
  unsigned int* sample = new unsigned int [p.size()];
  for (int i = 0; i< p.size(); i++){
    sample[i]=0;
  }
  r.sampleMultinomial(p.size(), sampleSize, pArray, sample);

  for (int i = 0; i < p.size(); i++){
    double sampledFrequency = (double) sample[i]/ (double) sampleSize;
    int degree = getDegree(ids.at(i));
    if (bArray[i] == -1){
      //print in one file
      fsr <<simTime<<' '<<p[i]<<' '<<sampledFrequency<<' '<<age[i]<<' '<<myModelRef.tau(p[i], N)<<' '<<ids.at(i)<<' '<<degree<<endl;
    }else{
      double balancingRatio = (double) bArray[i] / age[i] ;
      //print in other file
      fsb <<simTime<<' '<<p[i]<<' '<<sampledFrequency<<' '<<age[i]<<' '<<myModelRef.tau(p[i], N)<<' '<<ids.at(i)<<' '<< degree<<' '<<balancingRatio <<endl;
    }
  }
  return 1;
}
*/
int population::isFixed(int id){
  int present = 0;
  if(alleles.size() == 1){ //if someone is fixed
    if ( alleles.front().getId() == id){
      //      cout <<id<<" fixed 1"<<endl;
      return 1; //if the specific one is fixed
    }else{
      //      cout<<alleles.front().getId()<<" fixed 2"<<endl;
      return 2; //someone else is fixed
    }
  }else{
    return 0; //no one is fixed
  }
}

void population::mutate(randomv &r, modelFunctions &myModelRef, environment &envRef, ofstream &fse){
  //find how many mutations, if any
  int mutationNumber = r.sampleBinomial(u,ploidy*N);///haploid or diploid
  if (mutationNumber == 0){
    return;
  }
  //fse<<u<<" "<<ploidy<<" "<<N<<" "<<mutationNumber<<" "<<alleles.size()<<endl;
  double newAlleleFrequency = 1./(ploidy*N); //haploid or diploid
  int alleleNumber = alleles.size();
  //find which alleles will be mutated and how many times (for multiple hits)
  int i = 0;
  double Pmu[alleleNumber];
  for (vector<allele>::iterator it = alleles.begin();
       it != alleles.end();++it){
    Pmu[i] = (*it).getP();
    i++;
  }
  unsigned int* mutations = new unsigned int [alleleNumber];
  r.sampleMultinomial(alleles.size(),mutationNumber,Pmu,mutations);
  vector<allele> mutants;
  int counter = 0;
  for (vector<allele>::iterator it = alleles.begin();it != alleles.end();++it){
    if(mutations[counter]==0){
      counter++;
      continue;
    }
    double pi = (*it).getP();
    vector<double> rv = (*it).getRv();
    int id = (*it).getId();
    //create new mutation vector
    for(int i = 0; i<mutations[counter]; i++){
	//fse<<id<<" "<<counter<<" "<<mutations[counter]<<endl;
      vector<double> mutation = myModelRef.getMutationVector(r, myModelRef.getNumDimensions());
      vector<double> newRv = myModelRef.add(mutation, rv);
      allele newAllele(newAlleleFrequency,newRv, myModelRef);
      newAllele.setId(alleleCounter);
      newAllele.setIdp(id);
      newAllele.setBd(simTime);
      double wmutated = envRef.fW(newRv,myModelRef);
     // (*it).setfW(envRef.fW(rvi,myModelRef)); //!recalculate the fitness every generation in case the optimum moved
      double wmutant = (*it).getfW();
      newAllele.setfW(wmutated);
      newAllele.setPhe(wmutated - wmutant); // new -old
      newAllele.setDistance(mutation[0]);
      //if repeats >> 1 can produce very large output
      fse<<id<<' '<<alleleCounter<<' '<<simTime;
      for(int j=0; j<mutation.size();++j){
	      fse<<' '<<mutation[j];
      }
      for(int j=0; j<rv.size();++j){
	      fse<<' '<<rv[j];
      }
      for(int j=0; j<newRv.size();++j){
	      fse<<' '<<newRv[j];
      }
      
      fse<<' '<<wmutant<<' '<<wmutated<<endl; //parent
      //fse<<id<<' '<<alleleCounter<<' '<<simTime<<endl; //parent offspring time
      newAllele.addSNP(alleleCounter);
      alleleCounter++;
      //reduce the frequency of the mutated one
      if (pi > newAlleleFrequency){
	(*it).setP(pi-newAlleleFrequency);
      }else{
	(*it).setP(0);
      }
      mutants.push_back(newAllele);  
    }
    
      counter++;
  }
  delete [] mutations;
  for (vector<allele>::iterator it = mutants.begin();
       it != mutants.end();++it){
    alleles.push_back(*it);
  }

  // vector<int> mutants;
  // for (int c = 0; c < mutationNumber; c++){
  //   mutants.push_back(r.sampleUniformInt(alleleNumber));
  // }
  // //An allele can mutate more than once
  // //make mutations
  // for (int i = 0; i < mutants.size(); i++){
  // //   int location = mutants.at(i);
  // //   double pi = alleles.at(location).getP();
  //   int id = alleles.at(location).getId();
  //   vector<double> rv =  alleles.at(location).getRv();
  //   //create new mutation vector
  //   vector<double> mutation;
  //   double newR = myModelRef.fR(r);
  //   double newTheta = myModelRef.fTheta(r);
  //   mutation.push_back(newR);
  //   mutation.push_back(newTheta);
  //   vector<double> newRv = myModelRef.add(mutation, rv);
  //   allele newAllele(newAlleleFrequency,newRv, myModelRef);
  //   newAllele.setId(alleleCounter);
  //   newAllele.setIdp(id);
  //   newAllele.setBd(simTime);
  //   double wmutated = envRef.fW(newRv,myModelRef);
  //   double wmutant = envRef.fW(rv,myModelRef);
  //   newAllele.setPhe(wmutated - wmutant); // new -old
  //   newAllele.setDistance(myModelRef.vdistance(newRv,rv));
  //   //if repeats >> 1 can produce very large output
  //   fse<<id<<' '<<alleleCounter<<' '<<simTime<<endl; //parent offspring time
  //   alleleCounter++;
  //   //reduce the frequency of the mutated one
  //   if (pi > newAlleleFrequency){
  //     alleles.at(location).setP(pi-newAlleleFrequency);
  //   }else{
  //     alleles.at(location).setP(0);
  //   }
  //   alleles.push_back(newAllele);
  // }
}

void population::propagate(randomv &r, modelFunctions &myModelRef, environment &envRef){
  ///  vector<allele> temp;  //temporary copy (memory or speed?)
  int counter = 0;
  //build array of pd's 
  int alleleNumber = alleles.size();
  double pds[alleleNumber];
  //build array of offspring numbers
  unsigned int* newGeneration = new unsigned int [alleleNumber];
 // cout<<"propagate1"<<endl;
  //initialize
  for (int i = 0; i < alleleNumber; i++){
    pds[i] = 0;
    newGeneration[i] = 0;
  }
 // cout<<"propagate2 "<<alleles.size()<<endl;
  wildW = 0;
  //find number of descendants for each genotype
  for (vector<allele>::iterator it = alleles.begin();
       it != alleles.end();++it){
    double pi = (*it).getP();   
    vector<double> rvi = (*it).getRv();
   // for(int i=0; i<rvi.size();i++){
//	    cout<<"allele1 "<<rvi[i]<<endl;
   // }
    double pd = 0;
    int iid = (*it).getId();
   // (*it).setfW(envRef.fW(rvi,myModelRef)); //!recalculate the fitness every generation in case the optimum moved
    double wi = (*it).getfW();
    if(wi>wildW){
      wildW = wi;
    }
   // cout<<pi<<" "<<pd<<" "<<iid<<" "<<wi<<endl;
    //haploid or diploid
    if (ploidy == 2){
      for(vector<allele>::iterator jt = alleles.begin();
	  jt!= alleles.end();++jt){
	//      cout<<"Inner Loop: "<<endl;
	double pj = (*jt).getP();
	
	vector<double> rvj = (*jt).getRv();
     //  for(int i=0; i<rvj.size();i++){
//	    cout<<"allele2 "<<rvj[i]<<endl;
   // }
	int jid = (*jt).getId();
	//find fitness of diploid genotype
	double wij = envRef.diploidFitness(rvi, rvj, myModelRef);
	//compute new frequency after selection
	/*
	  wij = wji
	  if i != j pd += ((1/2)*pi*pj*wij+(1/2)*pj*pi*wji)/wij = pi*pj*wij/wm
	  if i == j pd += pi*pj*wij /wm
	*/
	pd += pi*pj*wij;
      }
    }else if (ploidy == 1){
      pd += pi*wi;
    }else {
      cerr <<"error 1"<<endl;
      exit(1);
    }
    if (wm != 0){
      pd = pd/wm; // divide once
    }else{
      cerr<<"error 2"<<endl;
      exit(1);
    }
    (*it).setP(pd);
    pds[counter] = pd;
    counter++;
  }
 // cout<<"propagate3"<<endl;
  counter = 0;
  //multinomial sampling to find out how many offspring
  r.sampleMultinomial(alleleNumber, ploidy*N, pds, newGeneration); 
  //create new generation with appropriate frequences
 // cout<<"propagate4"<<endl;
  for (vector<allele>::iterator it = alleles.begin();     
       it != alleles.end();){
    if (newGeneration[counter] == 0){ //allele extinction
      alleles.erase(it);
    }else{
      (*it).setP((double) newGeneration[counter]/ (double) (ploidy*N) ); //fixed population size
      ++it;
    }
    counter++;
  }
  delete[] newGeneration;
}

double population::updatewm(modelFunctions &myModelRef, environment &envRef){
  double sum = 0;
  double varsumSqr = 0; // for calculating the variance
  maxGenW = 0; //recalculate at each generation
  for (vector<allele>::iterator it = alleles.begin(); it !=alleles.end(); ++it){
    double pi = (*it).getP();
    vector<double> rvi = (*it).getRv();
    double wi = (*it).getfW();
    int iid = (*it).getId();
    if (ploidy == 2){ 	//haploid or diploid
      for (vector<allele>::iterator jt = it;  jt !=alleles.end(); ++jt){
	double pj = (*jt).getP();
	vector<double> rvj = (*jt).getRv();
	int jid = (*jt).getId();
	double wij = envRef.diploidFitness(rvi, rvj, myModelRef);
	if(wij > maxWij){
	  maxWij = wij;
	}
	if(wij > maxGenW){
	  maxGenW = wij;
	}
	if (iid != jid){
	  varsumSqr += 2*pi*pj*wij*wij;
	  sum += 2*pi*pj*wij; //== pi*pj*wij + pj*pi*wji, wij=wji
	}else{
	  varsumSqr += pi*pj*wij*wij;
	  sum += pi*pj*wij; //== pi*pj*wij, wij = wji
	}
      }
    }else if (ploidy == 1){
      if(wi > maxGenW){
	maxGenW = wi;
      }
      varsumSqr += pi*wi*wi;
      sum += pi*wi;
    }else {
      cerr<<"error 3"<<endl;
      exit(1);
    }
  }
  if(alleles.size() == 1){
    var = 0;
  }else{
    var = varsumSqr - sum*sum;
  }
  return sum;
}

int population::isAdapted(modelFunctions &myModelRef, environment &envRef){
  // find if 90% of the population is closer than r0
  //double loop through alleles
  double sum = 0; //proportion of individuals closer to the optimum than r0
  for (vector<allele>::iterator it = alleles.begin();
       it !=alleles.end(); ++it){
    double pi = (*it).getP();
    vector<double> rvi = (*it).getRv();
    int iid = (*it).getId();
    if (ploidy == 2){ 	//haploid or diploid
      for (vector<allele>::iterator jt = it;
	   jt !=alleles.end(); ++jt){
	vector<double> rvj = (*jt).getRv();
	vector<double> hetRv = myModelRef.add(rvi, rvj);
	hetRv.front()/=2.;
	if(hetRv.front()<r0){
	  double pj = (*jt).getP();
	  int jid = (*jt).getId();
	  if (iid != jid){
	    sum += 2*pi*pj;
	  }else{
	    sum += pi*pj;
	  }
	}
      }
    }else if (ploidy == 1){
      if(myModelRef.vdistance(rvi, envRef.getOptimum(),0) < r0){
	sum+=pi;
      }
    }
  }
  //  return sum;
  if(sum>0.9){
    return 1;
  }else{
    return 0;
  }
  //add p^2 +2pqs
  //r0;
  //c*sqrt(2./(double)N);
}

int population::ispseudoFixed(double threshold){
  /*
  //FREQUENCY
  for (vector<allele>::iterator it = alleles.begin();
       it !=alleles.end(); ++it){
    //if some allele, besides the initial one, reaches to majority frequency
    if((*it).getP()>(threshold) && (*it).getId() != 0){
	return 1;
    }
  }
  return 0;
  */
  /*
  //VARIANCE !
  if (var < 0.01 && simTime > 5){
    return 1;
  }else{
    return 0;
  } 
  */

  /*
  //GENETIC LOAD 
  double load = (wildW - wm)/wildW;
  if (load < 0.1 && simTime > 2){
    return 1;
  }else{
    return 0;
  } 
  */
  /*
  //MAXIMUM INDIVIDUAL FITNESS ! 
  if (wildW > 1.99){
    return 1;
  }else{
    return 0;
  } 
  */
  //MEAN FITNESS !
  if (wm > threshold){
    return 1;
  }else{
    return 0;
  }
  /*  
  //FIND MINIMUM VALUE ! Not Verified for haploid !

  double mymin = alleles.front().getRv().front();
    for (vector<allele>::iterator it = alleles.begin();
       it !=alleles.end(); ++it){
      if (mymin > (*it).getRv().front())
	{
	  mymin = (*it).getRv().front();
	}
      if (mymin < rmin && (*it).getP() > 0.99){
	return 1;
      }
    }
    return 0;
  */

 /*

 //   FIND MEDIAN VALUE ! Not Verified for haploid !

  //iterate over alleles and make array of r's
  int size = alleles.size();
  double tempar[size];
  int i = 0 ;
  for (vector<allele>::iterator it = alleles.begin();
       it !=alleles.end(); ++it){
    tempar[i]=(*it).getRv().front();
    i++;
  }
  gsl_sort(tempar,1,size);
  //find median value
  double median;
  int middle = size/2;
  if(size%2 ==0){
    median = (tempar[middle]+tempar[middle+1])/2;
  }else{
    median = tempar[middle];
  }
  //compare and return
  if(median<rmin){
    return 1;
  }else{
    return 0;
  }
  */
}

bool population::areBalanced(modelFunctions &myModelRef, environment &envRef){
  //check if polymorphism is stable
  double maxWDiff = 0; //maximum difference between fitness of heterozygote and homozygote with max fitness max(wij-max(wi,wj))
  vector<double> wVector;
  polymorphisms = 0;
  vector<double> frequencies;
  //  vector<int> ages;
  for (vector<allele>::iterator it = alleles.begin();
       it !=alleles.end(); ++it){ 
    vector<double> rvi = (*it).getRv();
    double pi = (*it).getP();
    double wi = (*it).getfW();
    double d = 0.05;
    //    int age = simTime - (*it).getBd();
    if (pi > d && pi < (1-d)){
      frequencies.push_back(pi);
      //      ages.push_back(age);
      polymorphisms++;
      for(vector<allele>::iterator jt = it;
	  jt !=alleles.end(); ++jt){
	double pj = (*jt).getP();
	if(pj > d && pi < (1-d)){
	  vector<double> rvj = (*jt).getRv();
	  double wj = (*jt).getfW();
	  double wij = envRef.diploidFitness(rvi, rvj, myModelRef);
	  wVector.push_back(wij);
	  double diff = wij-max(wi,wj);
	  if(maxWDiff < diff){
	    maxWDiff = diff;
	  }
	}
      }
    }
  }
  if ((polymorphisms > 1) && (1/N < maxWDiff/maxGenW)){
    stability stbf(wVector,polymorphisms);
    if (stbf.stable() == 1){
      stableEquilibrium = 1;
      //      for (int i = 0; i<ages.size(); i++){ 
      //	fsr<<simTime<<" "<<ages.at(i)<<" "<<frequencies.at(i)<<endl;
      //      }
    }else{
      stableEquilibrium = 0;
    }
  }else{
    stableEquilibrium = 0;
  }
  return stableEquilibrium;
}

/*
void population::addb (int a, int b){

  //if table size is zero
  // add first two lines (A B
  //                      B A)

  vector<int> alpha(1,a);
  alpha.push_back(b);
  vector<int> beta(1,b);
  beta.push_back(a);
  if(table.size() == 0){
    table.push_back(alpha);
    table.push_back(beta);
  }

  bool foundA = false;
  bool foundB = false;
  //loop through vector (it)
  for(vector<vector <int> >::iterator it = table.begin();
      it != table.end(); ++it){
    bool unique = true;
    if (a == (*it).front()){
      foundA = true;
      for(vector <int>::iterator jt = (*it).begin() + 1;
	  jt != (*it).end(); ++jt){
	if(b == *jt){
	  unique = false;
	  break;
	}
      }
      if(unique == true){
	(*it).push_back(b);
      }
    }else if (b == (*it).front()){
      foundB = true;
      for(vector <int>::iterator jt = (*it).begin() + 1;
	  jt != (*it).end(); ++jt){
	if(a == *jt){
	  unique = false;
	  break;
	}
      }
      if(unique == true){
	(*it).push_back(a);
      }
    }
  }
  if (foundA == false){
    table.push_back(alpha);
  }
  if (foundB == false){
    table.push_back(beta);
  }
}

void population::printTable(){
  for(vector<vector <int> >::iterator it = table.begin();
      it != table.end(); ++it){
    for(vector <int>::iterator jt = (*it).begin();
	jt != (*it).end(); ++jt){
      cout << *jt<< ' ';
    }
    cout<<endl;
  }
  cout<<endl;
}
*/
/*
int population::getDegree(int id){
    for(vector<vector <int> >::iterator it = table.begin();
      it != table.end(); ++it){
      if ((*it).front() == id){
	return (*it).size();
      }
    }
    return 0;
}
*/
/*
double population::getBalancedRatio(void){
  double br = 0;
  int totalAlleles = alleles.size();
  if(totalAlleles != 0){
    br = (double) 2*balancedPairs/((double)totalAlleles*(totalAlleles-1));
  }
  return br;
}
*/
int population::updatePhm(modelFunctions &myModelRef){
	vector<double> sum;
	vector<double> varsumSqr; // for calculating the variance
	for(int i=0; i<myModelRef.getNumDimensions();i++){
		varsumSqr.push_back(0);
		sum.push_back(0);
	}
//	cout<<"Phm1"<<endl;
	for (vector<allele>::iterator it = alleles.begin(); it !=alleles.end(); ++it){
		double pi = (*it).getP();
		vector<double> xyi = (*it).getCartesian(myModelRef);
		if (ploidy == 2){ 	//haploid or diploid
			int iid = (*it).getId();
			for (vector<allele>::iterator jt = it; jt !=alleles.end(); ++jt){
				double pj = (*jt).getP();
				vector<double> xyj = (*jt).getCartesian(myModelRef);
				double pipj = pi*pj;
				int jid = (*jt).getId();
				if (iid != jid){
					for(int i=0; i<myModelRef.getNumDimensions();i++){
						varsumSqr[i] += 2*pipj*(xyi[i] + xyj[i])*(xyi[i] + xyj[i]);
						sum[i] += 2*pipj*(xyi[i] + xyj[i]);
					}
				}else{
					for(int i=0; i<myModelRef.getNumDimensions();i++){
						varsumSqr[i] += pipj*xyi[i]*xyi[i];
						sum[i] += pipj*xyi[i];
					}
				}
			}
		}
		else if (ploidy == 1){
			for(int i=0; i<myModelRef.getNumDimensions();i++){
				varsumSqr[i] += pi*xyi[i]*xyi[i];
				sum[i] += pi*xyi[i];
			}
		}else {
			cerr<<"error 3"<<endl;
			exit(1);
		}	
	}
	//cout<<"Phm2"<<endl;
	for(int i=0; i<myModelRef.getNumDimensions();i++){
		phenotypeVar[i]=varsumSqr[i]-sum[i]*sum[i];
	}
//	cout<<"Phm3"<<endl;
	phm=myModelRef.cartesianToPolar(sum);
//	cout<<"Phm4"<<endl;
	
	return 1;
}