#include <iomanip>
#include <iostream>
#include <math.h>
#include <regex>
#include <sstream>
#include <string>
#include <time.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <algorithm> 
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <fpu_control.h>
#include "environment.h"
#include "modelFunctions.h"
#include "randomv.h"


using namespace std;
using namespace boost;

template<typename Type>
string joinArray(vector<Type> vec, char delim)
{
   std::stringstream ss;
	for(size_t i = 0; i < vec.size(); ++i)
	{
	if(i != 0)
		ss << delim;
	ss << vec[i];
	}
	return ss.str();
}



double invasionProbDiploid(double s1, double h1, double s2, double h2, double h3, double q, int N, double expectedEqFreq, int numSims){
	int t = 10000;
	
	double r = 1./(2*N);
	double p = 1-q-r;
	//cout<<s1<<"\t"<<h1<<"\t"<<s2<<"\t"<<h2<<"\t"<<h3<<"\t"<<p<<"\t"<<q<<"\t"<<r<<endl;
	double established =0;
	double currentMeanFitness = p*p+2*p*q*(1+s1*h1)+q*q*(1+s1);
	double hetFit = ((1+s2*h2)*p+(1+(s1+s2)/2*h3)*q)/currentMeanFitness-1;
	
	double establishmentFreq = 1./(2*N*hetFit);
	
	if(establishmentFreq<0){
		//cout<<0<<endl;
		return 0;
	}
	
	const gsl_rng_type *Gt; 
	gsl_rng *Gr;
	gsl_rng_env_setup();
	gsl_rng_default_seed = (unsigned long)time(0)*getpid();
	Gt = gsl_rng_default;
	Gr = gsl_rng_alloc (Gt);
	
	
	
	for(int sim=0;sim<numSims;sim++){
		
		double genotypeFreqs [6] = {p*p,2*p*q, q*q, 2*p*r, 2*q*r, r*r};
		
		
		for(int gen=0; gen<t; ++gen){
			//selection
			double meanFitness = genotypeFreqs[0]+genotypeFreqs[1]*(1+s1*h1)+genotypeFreqs[2]*(1+s1)+genotypeFreqs[3]*(1+s2*h2)+genotypeFreqs[4]*(1+(s1+s2)/2*h3) + genotypeFreqs[5]*(1+s2);
			genotypeFreqs[0]*=1./meanFitness;
			genotypeFreqs[1]*=(1+s1*h1)/meanFitness;
			genotypeFreqs[2]*=(1+s1)/meanFitness;
			genotypeFreqs[3]*=(1+s2*h2)/meanFitness;
			genotypeFreqs[4]*=(1+(s1+s2)/2*h3)/meanFitness;
			genotypeFreqs[5]*=(1+s2)/meanFitness;
			
			//drift
			unsigned int* genotypes = new unsigned int [6];
			gsl_ran_multinomial(Gr, 6, N, genotypeFreqs, genotypes);
			for(int idex=0;idex<6;idex++){
				genotypeFreqs[idex]= genotypes[idex]*1./N;				
			}
			
			
			//mating
			double newP = genotypeFreqs[0]+genotypeFreqs[1]/2 + genotypeFreqs[3]/2;
			double newQ = genotypeFreqs[2]+genotypeFreqs[1]/2 + genotypeFreqs[4]/2;
			double newR = 1-newP-newQ;
			if(newR==0){
				//cout<<"LOST!\t"<<gen<<"\t"<<genotypeFreqs[0]<<"\t"<<genotypeFreqs[1]<<"\t"<<genotypeFreqs[2]<<"\t"<<genotypeFreqs[3]<<"\t"<<genotypeFreqs[4]<<"\t"<<genotypeFreqs[5]<<endl;
				break;
			}
			if(newR>= 0.9 * expectedEqFreq){
				//cout<<"ESTABLISHED!\t"<<gen<<"\t"<<genotypeFreqs[0]<<"\t"<<genotypeFreqs[1]<<"\t"<<genotypeFreqs[2]<<"\t"<<genotypeFreqs[3]<<"\t"<<genotypeFreqs[4]<<"\t"<<genotypeFreqs[5]<<endl;
				established++;
				break;
			}
			genotypeFreqs[0]=newP*newP;
			genotypeFreqs[1]=2*newP*newQ;
			genotypeFreqs[2]=newQ*newQ;
			genotypeFreqs[3]=2*newP*newR;
			genotypeFreqs[4]=2*newQ*newR;
			genotypeFreqs[5]=newR*newR;
			delete[] genotypes;
			
		}
		
		
	}
	return established/numSims;
}


vector<double> computeMaxMeanFitness(double homOne, double het, double homTwo){
	vector<double> result;
	result.push_back(0.0);
	result.push_back(0.0);
	
	double s = homTwo/homOne-1;
	double h=0;
	if(s!=0){
		h = (het/homOne-1)/s;
	}
	
	//cout<<"computing max mean fitness with values "<<s<<" "<<h<<endl;
	
	if(s>0){
		if(h>0 && h<1){ //incomplete dominance, p is lost
			result[0] = homTwo;
		}
		if(h<0){ //underdominance, p is fixed
			result[0] = homOne;
			result[1] = 1;
		}
		if(h>1){ //overdominance
			double balancedFreq = (h-1)/(2*h-1);
			//print "BALANCED:\tbalancedFreq\n";
			result[0] = pow(balancedFreq,2)*homOne + 2*(1-balancedFreq)*(balancedFreq)*het + pow((1-balancedFreq),2) * homTwo;
			result[1] = balancedFreq;
		}
	}else{
		if(h>0){ //mutation is deleterious in all cases, p is fixed
			result[0] = homOne;
			result[1] = 1;
		}
		if(h<0){ //overdominance, p is balanced
			double balancedFreq = (h-1)/(2*h-1);
			//print "BALANCED:\tbalancedFreq\n";
			result[0] = pow(balancedFreq,2)*homOne + 2*(1-balancedFreq)*(balancedFreq)*het + pow((1-balancedFreq),2) * homTwo;
			result[1] = balancedFreq;
		}
	}
	return result;
}


double effectiveNumPaths;
unordered_map<string, double> invasionProbMap;
unordered_map<string, int> validEndStates;
vector< vector<double> > mutVectors;
vector<double> optimum;
vector<double> initialPosition;
int loopCount;
int numPaths;
int validPaths;
double sumProb;
int N, d, pathLength;
vector<double> probArray;



void computeHeterozygousWeinreichPaths(vector<int> oneAlleles, vector<int> twoAlleles, vector<double> oneAlleleVector, double oneFreq, vector<double> twoAlleleVector, double twoFreq, vector<string> observedPoints, vector<double> observedMuts, vector<double> currentFitness, int depth, double probability, int numOverdominant, vector<int> alleleMutated, string previousStates, vector<string> meanPopulationPhenotype, modelFunctions  &myModelRef, environment	 &envRef){
	//cout<<"\n\nInside function to compute het weinreich paths at depth "<<depth<<" and probability "<<probability<<" num overdom: "<<numOverdominant<<" previous states: "<<previousStates<<endl;
	//cout<<"one alleles: "<<joinArray(oneAlleles,'_')<<"\t"<<joinArray(oneAlleleVector,'_')<<"\t"<<oneFreq<<endl;
	//cout<<"two alleles: "<<joinArray(twoAlleles,'_')<<"\t"<<joinArray(twoAlleleVector,'_')<<"\t"<<twoFreq<<"\n"<<endl;
	//cout<<"observed points: "<<joinArray(observedPoints,':')<<endl;
	//cout<<"observed muts: "<<joinArray(observedMuts,'_')<<endl;
	//cout<<"currentFitness: "<<joinArray(currentFitness,'_')<<endl;
	//cout<<"alleleMutated: "<<joinArray(alleleMutated,'_')<<endl;
	//cout<<"meanPopulationPhenotype: "<<joinArray(meanPopulationPhenotype,':')<<endl;
	
	
	vector<double> initHetPhenotype = myModelRef.add(oneAlleleVector, twoAlleleVector);
	initHetPhenotype[0]/=2;
	double initHomOneFreq = pow(oneFreq,2);
	double initHomTwoFreq = pow(twoFreq,2);
	double initHetFreq = 2 * oneFreq * twoFreq;
	vector<vector<double> > initAlleleVec;
	initAlleleVec.push_back(oneAlleleVector);
	initAlleleVec.push_back(twoAlleleVector);
	initAlleleVec.push_back(initHetPhenotype);
	vector<double> initAlleleFreqVec;
	initAlleleFreqVec.push_back(initHomOneFreq);
	initAlleleFreqVec.push_back(initHomTwoFreq);
	initAlleleFreqVec.push_back(initHetFreq);
	vector<double> initPopMeanPhenotype = myModelRef.getMeanPhenotype(initAlleleVec,initAlleleFreqVec);
	meanPopulationPhenotype.push_back(joinArray(initPopMeanPhenotype,'_'));
	if(probability < pow(10,-5) || oneFreq == 0){
		//cout<<"probability is "<<probability<<" and frequency is "<<oneFreq<<" so we are quitting!"<<endl;
		return;
	}
	if(mutVectors.size() == oneAlleles.size() || mutVectors.size() == twoAlleles.size()){
		numPaths++;
		vector<double> mutString;
		
		
		
		if(mutVectors.size() == oneAlleles.size()){
			for(int i=0; i<oneAlleles.size(); i++){
				mutString.push_back(mutVectors[oneAlleles[i]][0]);
			}
			cout<<"Path Found:\t"<<joinArray(observedPoints,':')<<"\t"<<probability<<"\t"<<joinArray(currentFitness,'_')<<"\t"<<joinArray(mutString,'_')<<"\t"<<numOverdominant<<"\t"<<depth;
			cout<<"\t"<<" "+joinArray(oneAlleles,'_')+" "<<"\t"<<'_'+joinArray(mutString,'_')<<"\t"<<joinArray(twoAlleles,'_')<<"\t"<<oneFreq<<"\t"<<twoFreq<<"\t"<<joinArray(alleleMutated,'_')<<"\t"<<':'+joinArray(meanPopulationPhenotype,':')<<endl;
		}else{
			for(int i=0; i<twoAlleles.size(); i++){
				mutString.push_back(mutVectors[twoAlleles[i]][0]);
			}
			cout<<"Path Found:\t"<<joinArray(observedPoints,':')<<"\t"<<probability<<"\t"<<joinArray(currentFitness,'_')<<"\t"<<joinArray(mutString,'_')<<"\t"<<numOverdominant<<"\t"<<depth;
			cout<<"\t"<<" "+joinArray(twoAlleles,'_')+" "<<"\t"<<'_'+joinArray(mutString,'_')<<"\t"<<joinArray(oneAlleles,'_')<<"\t"<<twoFreq<<"\t"<<oneFreq<<"\t"<<joinArray(alleleMutated,'_')<<"\t"<<':'+joinArray(meanPopulationPhenotype,':')<<endl;
		}
		probArray.push_back(probability);
		sumProb+=probability;
		vector<int> validEndR;
		if(oneFreq > 0){
			validEndR.push_back(int(oneAlleleVector[0]*10000));
		}
		
		if(twoFreq > 0){
			validEndR.push_back(int(twoAlleleVector[0]*10000));
		}
		sort (validEndR.begin(), validEndR.end()); 
		validEndStates[joinArray(validEndR,'_')] = 1;
		return;
	}
	
	unordered_map<int, double> invasionProbHashOne;
	unordered_map<int, double> invasionProbHashTwo;
	double totalInvasionProb = 0;
	
	for(int i=0; i<mutVectors.size(); i++){
		vector<double> mutThetas = mutVectors[i];
		double mutR = mutThetas[0];
		//cout<<"current mutR is "<<mutR<<endl;
		if(find(oneAlleles.begin(), oneAlleles.end(),i) == oneAlleles.end()){
			vector<double> newPositionHomOne = myModelRef.add(oneAlleleVector, mutThetas);
			mutThetas[0]/=2;
			vector<double> newPositionHetOneA =  myModelRef.add(oneAlleleVector, mutThetas);
			vector<double> newPositionHetOneB =  myModelRef.add(newPositionHomOne, twoAlleleVector);
			newPositionHetOneB[0]/=2;
			vector<double> currentHetPosition = myModelRef.add(oneAlleleVector, twoAlleleVector);
			currentHetPosition[0]/=2;
			double fitHomOne = envRef.fW(newPositionHomOne,myModelRef);
			//cout<<fitHomOne<<endl;
			double fitStartOne =  envRef.fW(oneAlleleVector,myModelRef);
			//cout<<fitStartOne<<endl;
			double fitStartTwo =  envRef.fW(twoAlleleVector,myModelRef);
			//cout<<fitStartTwo<<endl;
			double fitHetOneA =  envRef.fW(newPositionHetOneA,myModelRef);
			//cout<<fitHetOneA<<endl;
			double fitHetOneB =  envRef.fW(newPositionHetOneB,myModelRef);
			//cout<<fitHetOneB<<endl;
			double fitCurrentHet =  envRef.fW(currentHetPosition,myModelRef);
			//cout<<fitCurrentHet<<endl;
			mutThetas[0]*=2;
			
			vector<double> fitOneVec = computeMaxMeanFitness(fitStartOne, fitHetOneA, fitHomOne);
			vector<double> fitTwoVec = computeMaxMeanFitness(fitStartTwo, fitHetOneB, fitHomOne);
			double eqFreq = 0;
			
			if(fitOneVec[0]>fitTwoVec[0]){
				eqFreq = 1 - fitOneVec[1];
			}else{
				eqFreq = 1 - fitTwoVec[1];
			}
			//cout<<myModelRef.getFitnessFunction()<<" "<<myModelRef.getA()<<" "<<myModelRef.getC()<<" "<<myModelRef.getD()<<" "<<envRef.fW(oneAlleleVector,myModelRef)<<endl;
			//cout<<"current positions are: "<<joinArray(newPositionHomOne,'_')<<"\t"<<joinArray(oneAlleleVector,'_')<<"\t"<<joinArray(twoAlleleVector,'_')<<"\t"<<joinArray(newPositionHetOneA,'_')<<"\t"<<joinArray(newPositionHetOneB,'_')<<"\t"<<joinArray(currentHetPosition,'_')<<endl;
			//cout<<"computed fitnesses are: "<<fitHomOne<<"\t"<<fitStartOne<<"\t"<<fitStartTwo<<"\t"<<fitHetOneA<<"\t"<<fitHetOneB<<"\t"<<fitCurrentHet<<endl;
			fitHomOne/=fitStartOne;
			fitStartTwo/=fitStartOne;
			fitHetOneA/=fitStartOne;
			fitCurrentHet/=fitStartOne;
			fitHetOneB/=fitStartOne;
			fitStartOne=1;
			//cout<<"computed fitnesses post normalization are: "<<fitHomOne<<"\t"<<fitStartOne<<"\t"<<fitStartTwo<<"\t"<<fitHetOneA<<"\t"<<fitHetOneB<<"\t"<<fitCurrentHet<<endl;
			
			double currentMeanFitness = fitStartOne*oneFreq*oneFreq + fitCurrentHet*2*oneFreq*twoFreq + fitStartTwo*twoFreq*twoFreq;
			
			double s1 = fitStartTwo-1;
			double h1 = 0;
			if(s1!=0){h1=(fitCurrentHet-1)/s1;}
			double s2 = fitHomOne-1;
			double h2 = 0;
			if(s2!=0){h2=(fitHetOneA-1)/s2;}
			double h3 = 0;
			if(s2+s1!=0){h3=(fitHetOneB-1)/((s2+s1)/2);}

			double estimatedHetFitness = (fitHetOneA*oneFreq+fitHetOneB*twoFreq)/currentMeanFitness-1;
			
			double invasionProb =0;
			ostringstream probMapStringStream;
			probMapStringStream<<" "<<s1<<" "<<h1<<" "<<s2<<" "<<h2<<" "<<h3<<" "<<twoFreq<<" "<<N<<" ";
			string probMapString (probMapStringStream.str());
			if(estimatedHetFitness<0){
				invasionProb=0;
			}
			else{
				if(invasionProbMap.find(probMapString) != invasionProbMap.end()){
					invasionProb = invasionProbMap[probMapString];
				}else{
					invasionProb = invasionProbDiploid(s1,h1, s2, h2, h3, twoFreq, N, eqFreq, 10000);
					invasionProbMap[probMapString] = invasionProb;
				}
			}
			invasionProbHashOne[i] = invasionProb;
			totalInvasionProb += invasionProb;
			//cout<<"mutating first allele invasion prob is "<<invasionProb<<endl;
			//cout<<"invasionProbMapKey is "<<probMapString<<endl;
		}
		if(twoFreq>0 && find(twoAlleles.begin(),twoAlleles.end(),i) == twoAlleles.end()){ //if we are mutating both alleles and allele 2 has nonzero frequency and does not have this mutation, add mutation in and compute invasion prob
			
			
			vector<double> newPositionHomOne = myModelRef.add(twoAlleleVector, mutThetas);
			mutThetas[0]/=2;
			vector<double> newPositionHetOneA =  myModelRef.add(twoAlleleVector, mutThetas);
			vector<double> newPositionHetOneB =  myModelRef.add(newPositionHomOne, oneAlleleVector);
			newPositionHetOneB[0]/=2;
			vector<double> currentHetPosition = myModelRef.add(twoAlleleVector, oneAlleleVector);
			currentHetPosition[0]/=2;
			double fitHomOne = envRef.fW(newPositionHomOne,myModelRef);
			double fitStartOne =  envRef.fW(oneAlleleVector,myModelRef);
			double fitStartTwo =  envRef.fW(twoAlleleVector,myModelRef);
			double fitHetOneA =  envRef.fW(newPositionHetOneA,myModelRef);
			double fitHetOneB =  envRef.fW(newPositionHetOneB,myModelRef);
			double fitCurrentHet =  envRef.fW(currentHetPosition,myModelRef);
			mutThetas[0]*=2;
			
			vector<double> fitOneVec = computeMaxMeanFitness(fitStartOne, fitHetOneA, fitHomOne);
			vector<double> fitTwoVec = computeMaxMeanFitness(fitStartTwo, fitHetOneB, fitHomOne);
			double eqFreq = 0;
			
			if(fitOneVec[0]>fitTwoVec[0]){
				eqFreq = 1 - fitOneVec[1];
			}else{
				eqFreq = 1 - fitTwoVec[1];
			}
			//cout<<myModelRef.getFitnessFunction()<<" "<<myModelRef.getA()<<" "<<myModelRef.getC()<<" "<<myModelRef.getD()<<" "<<envRef.fW(oneAlleleVector,myModelRef)<<endl;
			//cout<<"current positions are: "<<joinArray(newPositionHomOne,'_')<<"\t"<<joinArray(oneAlleleVector,'_')<<"\t"<<joinArray(twoAlleleVector,'_')<<"\t"<<joinArray(newPositionHetOneA,'_')<<"\t"<<joinArray(newPositionHetOneB,'_')<<"\t"<<joinArray(currentHetPosition,'_')<<endl;
			//cout<<"computed fitnesses are: "<<fitHomOne<<"\t"<<fitStartOne<<"\t"<<fitStartTwo<<"\t"<<fitHetOneA<<"\t"<<fitHetOneB<<"\t"<<fitCurrentHet<<endl;
			
			fitHomOne/=fitStartTwo;
			fitStartOne/=fitStartTwo;
			fitHetOneA/=fitStartTwo;
			fitCurrentHet/=fitStartTwo;
			fitHetOneB/=fitStartTwo;
			fitStartTwo=1;
			//cout<<"computed fitnesses post normalization are: "<<fitHomOne<<"\t"<<fitStartOne<<"\t"<<fitStartTwo<<"\t"<<fitHetOneA<<"\t"<<fitHetOneB<<"\t"<<fitCurrentHet<<endl;
			double currentMeanFitness = fitStartTwo*pow(oneFreq,2) + fitCurrentHet*2*oneFreq*twoFreq + fitStartOne*pow(twoFreq,2);
			
			double s1 = fitStartOne-1;
			double h1 = 0;
			if(s1!=0){h1=(fitCurrentHet-1)/s1;}
			double s2 = fitHomOne-1;
			double h2 = 0;
			if(s2!=0){h2=(fitHetOneB-1)/s2;}
			double h3 = 0;
			if(s2+s1!=0){h3=(fitHetOneA-1)/((s2+s1)/2);}

			double estimatedHetFitness = (fitHetOneB*oneFreq+fitHetOneA*twoFreq)/currentMeanFitness-1;
			ostringstream probMapStringStream;
			probMapStringStream<<" "<<s1<<" "<<h1<<" "<<s2<<" "<<h2<<" "<<h3<<" "<<twoFreq<<" "<<N<<" ";
			string probMapString (probMapStringStream.str());
			
			double invasionProb =0;
			if(estimatedHetFitness<0){
				invasionProb=0;
			}
			else{
				
				if(invasionProbMap.find(probMapString) != invasionProbMap.end()){
					invasionProb = invasionProbMap[probMapString];
				}else{
					invasionProb = invasionProbDiploid(s1,h1, s2, h2, h3, twoFreq, N, eqFreq, 10000);
					invasionProbMap[probMapString] = invasionProb;
				}
			}
			invasionProbHashTwo[i] = invasionProb;
			totalInvasionProb += invasionProb;
			//cout<<"mutating second allele invasion prob is "<<invasionProb<<endl;
			//cout<<"invasionProbMapKey is "<<probMapString<<endl;
			
		}
	}
	
	for(int i=0; i<mutVectors.size(); i++){ //for each mut
		vector<double> mutThetas =  mutVectors[i];
		double mutR = mutThetas[0];
		vector<int> oneAllelesCopy = oneAlleles;
		vector<int> twoAllelesCopy = twoAlleles;
		vector<double> observedMutsCopy = observedMuts;
		vector<double> currentFitnessCopy = currentFitness;
		vector<int> alleleMutatedCopy = alleleMutated;
		
		oneAllelesCopy.push_back(i);
		sort(oneAllelesCopy.begin(), oneAllelesCopy.end());
		sort(twoAllelesCopy.begin(), twoAllelesCopy.end());
		if(find(oneAlleles.begin(), oneAlleles.end(), i)==oneAlleles.end() && joinArray(oneAllelesCopy,'_') != joinArray(twoAllelesCopy,'_')){ //if allele one does not have the mutation, add it in, compute mean fitness and recurse with correct set of alleles
			
			vector<double> newPositionHomOne = myModelRef.add(oneAlleleVector, mutThetas);
			mutThetas[0]/=2;
			vector<double> newPositionHetOneA =  myModelRef.add(oneAlleleVector, mutThetas);
			vector<double> newPositionHetOneB =  myModelRef.add(newPositionHomOne, twoAlleleVector);
			newPositionHetOneB[0]/=2;
			vector<double> currentHetPosition = myModelRef.add(oneAlleleVector, twoAlleleVector);
			currentHetPosition[0]/=2;
			double fitHomOne = envRef.fW(newPositionHomOne,myModelRef);
			double fitStartOne =  envRef.fW(oneAlleleVector,myModelRef);
			double fitStartTwo =  envRef.fW(twoAlleleVector,myModelRef);
			double fitHetOneA =  envRef.fW(newPositionHetOneA,myModelRef);
			double fitHetOneB =  envRef.fW(newPositionHetOneB,myModelRef);
			double fitCurrentHet =  envRef.fW(currentHetPosition,myModelRef);
			mutThetas[0]*=2;
			
			vector<double> fitOneVec = computeMaxMeanFitness(fitStartOne, fitHetOneA, fitHomOne);
			vector<double> fitTwoVec = computeMaxMeanFitness(fitStartTwo, fitHetOneB, fitHomOne);
			
			//cout<<"\n\ndoing recursive steps first!"<<endl;
			//cout<<"one alleles: "<<joinArray(oneAlleles,'_')<<"\t"<<joinArray(oneAlleleVector,'_')<<"\t"<<oneFreq<<endl;
			//cout<<"two alleles: "<<joinArray(twoAlleles,'_')<<"\t"<<joinArray(twoAlleleVector,'_')<<"\t"<<twoFreq<<"\n"<<endl;

			//cout<<"computed fitnesses are: "<<fitHomOne<<"\t"<<fitStartOne<<"\t"<<fitStartTwo<<"\t"<<fitHetOneA<<"\t"<<fitHetOneB<<"\t"<<fitCurrentHet<<endl;
			//cout<<joinArray(fitOneVec,'_')<<"\t"<<joinArray(fitTwoVec,'_')<<endl;
			
			
			double currentMeanFitness = fitStartOne*pow(oneFreq,2) + fitCurrentHet*2*oneFreq*twoFreq + fitStartTwo*pow(twoFreq,2);
			
			
			fitHomOne/=fitStartOne;
			fitStartTwo/=fitStartOne;
			fitHetOneA/=fitStartOne;
			fitCurrentHet/=fitStartOne;
			fitHetOneB/=fitStartOne;
			fitStartOne=1;
			
			double s1 = fitStartTwo-1;
			double h1 = 0;
			if(s1!=0){h1=(fitCurrentHet-1)/s1;}
			double s2 = fitHomOne-1;
			double h2 = 0;
			if(s2!=0){h2=(fitHetOneA-1)/s2;}
			double h3 = 0;
			if(s2+s1!=0){h3=(fitHetOneB-1)/((s2+s1)/2);}
			
			//if(totalInvasionProb==0){
				
			//}else{
			double invasionProb=0;
			if(totalInvasionProb>0){
				invasionProb = invasionProbHashOne[i]/totalInvasionProb;
			}
			//cout<<joinArray(fitOneVec,'_')<<"\t"<<joinArray(fitTwoVec,'_')<<"\t"<<invasionProb<<endl;
			//cout<<"computed fitnesses post normalization are: "<<fitHomOne<<"\t"<<fitStartOne<<"\t"<<fitStartTwo<<"\t"<<fitHetOneA<<"\t"<<fitHetOneB<<"\t"<<fitCurrentHet<<endl;
			//cout<<"invasion prob: "<<invasionProb<<endl;
			if(fitOneVec[0]>=fitTwoVec[0] || twoFreq==0){
				double balancedFreq = fitOneVec[1];
				double balancedStateFound=0;
				if(balancedFreq!=0 && balancedFreq!=1){
					balancedStateFound=1;
				}
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				oneAllelesCopy.push_back(i);
				sort(oneAllelesCopy.begin(), oneAllelesCopy.end());
				sort(twoAllelesCopy.begin(), twoAllelesCopy.end());
				
				double done=0;
				string t1 = "\t"+joinArray(oneAllelesCopy,' ')+':'+joinArray(twoAllelesCopy,' ')+"\t";
				string t2 = "\t"+joinArray(twoAllelesCopy,' ')+':'+joinArray(oneAllelesCopy,' ')+"\t";
				smatch m1, m2;
				if(regex_search(previousStates,m1,regex(t1)) || regex_search(previousStates,m2,regex(t2))){
					continue;
				}
				string previousStatesCopy = previousStates;
				previousStatesCopy+=t1;
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				oneAllelesCopy.push_back(i);
				
				vector<string> observedPointsLocal2 = observedPoints;
				observedPointsLocal2.push_back(joinArray(newPositionHomOne,'_'));
				observedMutsCopy.push_back(mutR);
				currentFitnessCopy.push_back(fitOneVec[0]);
				alleleMutatedCopy.push_back(1);
				
				computeHeterozygousWeinreichPaths(oneAllelesCopy,oneAlleles,newPositionHomOne,1-balancedFreq,oneAlleleVector,balancedFreq,observedPointsLocal2,observedMutsCopy,currentFitnessCopy,depth+1,probability*invasionProb*oneFreq,numOverdominant+balancedStateFound, alleleMutatedCopy, previousStatesCopy, meanPopulationPhenotype, myModelRef, envRef);
				
			}else{
				double balancedFreq = fitTwoVec[1];
				double balancedStateFound=0;
				if(balancedFreq!=0 && balancedFreq!=1){
					balancedStateFound=1;
				}
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				oneAllelesCopy.push_back(i);
				sort(oneAllelesCopy.begin(), oneAllelesCopy.end());
				sort(twoAllelesCopy.begin(), twoAllelesCopy.end());
				
				double done=0;
				string t1 = "\t"+joinArray(oneAllelesCopy,' ')+':'+joinArray(twoAllelesCopy,' ')+"\t";
				string t2 = "\t"+joinArray(twoAllelesCopy,' ')+':'+joinArray(oneAllelesCopy,' ')+"\t";
				smatch m1, m2;
				if(regex_search(previousStates,m1,regex(t1)) || regex_search(previousStates,m2,regex(t2))){
					continue;
				}
				string previousStatesCopy = previousStates;
				previousStatesCopy+=t1;
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				oneAllelesCopy.push_back(i);
				
				vector<string> observedPointsLocal2 = observedPoints;
				observedPointsLocal2.push_back(joinArray(newPositionHomOne,'_'));
				observedMutsCopy.push_back(mutR);
				currentFitnessCopy.push_back(fitTwoVec[0]);
				alleleMutatedCopy.push_back(1);
				
				computeHeterozygousWeinreichPaths(oneAllelesCopy,twoAlleles,newPositionHomOne,1-balancedFreq,twoAlleleVector,balancedFreq,observedPointsLocal2,observedMutsCopy,currentFitnessCopy,depth+1,probability*invasionProb*oneFreq,numOverdominant+balancedStateFound, alleleMutatedCopy, previousStatesCopy, meanPopulationPhenotype, myModelRef, envRef);
			}
			//}
		}
		
		
		observedMutsCopy = observedMuts;
		currentFitnessCopy = currentFitness;
		alleleMutatedCopy = alleleMutated;
		oneAllelesCopy = oneAlleles;
		twoAllelesCopy = twoAlleles;
		twoAllelesCopy.push_back(i);
		sort(oneAllelesCopy.begin(), oneAllelesCopy.end());
		sort(twoAllelesCopy.begin(), twoAllelesCopy.end());
		if(find(twoAlleles.begin(),twoAlleles.end(),i)==twoAlleles.end() && joinArray(twoAllelesCopy,'_') != joinArray(oneAllelesCopy,'_')){ //edit this to make sense!
			vector<double> newPositionHomOne = myModelRef.add(twoAlleleVector, mutThetas);
			mutThetas[0]/=2;
			vector<double> newPositionHetOneA =  myModelRef.add(twoAlleleVector, mutThetas);
			vector<double> newPositionHetOneB =  myModelRef.add(newPositionHomOne, oneAlleleVector);
			newPositionHetOneB[0]/=2;
			vector<double> currentHetPosition = myModelRef.add(twoAlleleVector, oneAlleleVector);
			currentHetPosition[0]/=2;
			double fitHomOne = envRef.fW(newPositionHomOne,myModelRef);
			double fitStartOne =  envRef.fW(oneAlleleVector,myModelRef);
			double fitStartTwo =  envRef.fW(twoAlleleVector,myModelRef);
			double fitHetOneA =  envRef.fW(newPositionHetOneA,myModelRef);
			double fitHetOneB =  envRef.fW(newPositionHetOneB,myModelRef);
			double fitCurrentHet =  envRef.fW(currentHetPosition,myModelRef);
			mutThetas[0]*=2;
			
			vector<double> fitOneVec = computeMaxMeanFitness(fitStartOne, fitHetOneA, fitHomOne);
			vector<double> fitTwoVec = computeMaxMeanFitness(fitStartTwo, fitHetOneB, fitHomOne);
			//cout<<"\n\ndoing recursive steps second!"<<endl;
			//cout<<"one alleles: "<<joinArray(oneAlleles,'_')<<"\t"<<joinArray(oneAlleleVector,'_')<<"\t"<<oneFreq<<endl;
			//cout<<"two alleles: "<<joinArray(twoAlleles,'_')<<"\t"<<joinArray(twoAlleleVector,'_')<<"\t"<<twoFreq<<"\n"<<endl;
			
			
			
			double currentMeanFitness = fitStartTwo*oneFreq*oneFreq + fitCurrentHet*2*oneFreq*twoFreq + fitStartOne*twoFreq*twoFreq;
			
			
			fitHomOne/=fitStartTwo;
			fitStartOne/=fitStartTwo;
			fitHetOneA/=fitStartTwo;
			fitCurrentHet/=fitStartTwo;
			fitHetOneB/=fitStartTwo;
			fitStartTwo=1;
			
			double s1 = fitStartOne-1;
			double h1 = 0;
			if(s1!=0){h1=(fitCurrentHet-1)/s1;}
			double s2 = fitHomOne-1;
			double h2 = 0;
			if(s2!=0){h2=(fitHetOneB-1)/s2;}
			double h3 = 0;
			if(s2+s1!=0){h3=(fitHetOneA-1)/((s2+s1)/2);}
			
			double invasionProb=0;
			if(totalInvasionProb>0){
				invasionProb = invasionProbHashTwo[i]/totalInvasionProb;
			}
			//cout<<joinArray(fitOneVec,'_')<<"\t"<<joinArray(fitTwoVec,'_')<<"\t"<<invasionProb<<endl;
			if(fitOneVec[0]>=fitTwoVec[0]){
				double balancedFreq = fitOneVec[1];
				double balancedStateFound=0;
				if(balancedFreq!=0 && balancedFreq!=1){
					balancedStateFound=1;
				}
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				twoAllelesCopy.push_back(i);
				sort(oneAllelesCopy.begin(), oneAllelesCopy.end());
				sort(twoAllelesCopy.begin(), twoAllelesCopy.end());
				
				double done=0;
				string t1 = "\t"+joinArray(oneAllelesCopy,' ')+':'+joinArray(twoAllelesCopy,' ')+"\t";
				string t2 = "\t"+joinArray(twoAllelesCopy,' ')+':'+joinArray(oneAllelesCopy,' ')+"\t";
				smatch m1, m2;
				if(regex_search(previousStates,m1,regex(t1)) || regex_search(previousStates,m2,regex(t2))){
					continue;
				}
				string previousStatesCopy = previousStates;
				previousStatesCopy+=t1;
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				twoAllelesCopy.push_back(i);
				
				vector<string> observedPointsLocal2 = observedPoints;
				observedPointsLocal2.push_back(joinArray(newPositionHomOne,'_'));
				observedMutsCopy.push_back(mutR);
				currentFitnessCopy.push_back(fitOneVec[0]);
				alleleMutatedCopy.push_back(2);
				
				
				computeHeterozygousWeinreichPaths(twoAllelesCopy,twoAlleles,newPositionHomOne,1-balancedFreq,twoAlleleVector,balancedFreq,observedPointsLocal2,observedMutsCopy,currentFitnessCopy, depth+1,probability*invasionProb*twoFreq,numOverdominant+balancedStateFound, alleleMutatedCopy, previousStatesCopy, meanPopulationPhenotype, myModelRef, envRef);
			}else{
				double balancedFreq = fitTwoVec[1];
				double balancedStateFound=0;
				if(balancedFreq!=0 && balancedFreq!=1){
					balancedStateFound=1;
				}
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				twoAllelesCopy.push_back(i);
				sort(oneAllelesCopy.begin(), oneAllelesCopy.end());
				sort(twoAllelesCopy.begin(), twoAllelesCopy.end());
				
				double done=0;
				string t1 = "\t"+joinArray(oneAllelesCopy,' ')+':'+joinArray(twoAllelesCopy,' ')+"\t";
				string t2 = "\t"+joinArray(twoAllelesCopy,' ')+':'+joinArray(oneAllelesCopy,' ')+"\t";
				smatch m1, m2;
				if(regex_search(previousStates,m1,regex(t1)) || regex_search(previousStates,m2,regex(t2))){
					continue;
				}
				string previousStatesCopy = previousStates;
				previousStatesCopy+=t1;
				
				oneAllelesCopy = oneAlleles;
				twoAllelesCopy = twoAlleles;
				
				twoAllelesCopy.push_back(i);
				
				vector<string> observedPointsLocal2 = observedPoints;
				observedPointsLocal2.push_back(joinArray(newPositionHomOne,'_'));
				observedMutsCopy.push_back(mutR);
				currentFitnessCopy.push_back(fitTwoVec[0]);
				alleleMutatedCopy.push_back(2);
				
				computeHeterozygousWeinreichPaths(twoAllelesCopy,oneAlleles,newPositionHomOne,1-balancedFreq,oneAlleleVector,balancedFreq,observedPointsLocal2,observedMutsCopy,currentFitnessCopy, depth+1,probability*invasionProb*twoFreq,numOverdominant+balancedStateFound, alleleMutatedCopy, previousStatesCopy, meanPopulationPhenotype, myModelRef, envRef);
			}
		}
		
	}
	return;
}

/**ARGS
 * pathFile numDim N pathLength wFuncType gaussA gaussB gaussC gaussD
 */
int main(int argc, char* argv[]){
	if(argc < 10){
		cerr<<"Insufficient arguments!";
		exit(1);
	}
	fpu_control_t fpu_oldcw, fpu_cw;
	_FPU_GETCW(fpu_oldcw); // store old cw
	fpu_cw = (~_FPU_DOUBLE);
	_FPU_SETCW(fpu_cw);
	
	char*   pathFileName = argv[1]; //initial population description
	d = atoi(argv[2]);
	N = atoi(argv[3]);
	pathLength = atoi(argv[4]);
	char wFuncType = *argv[5];
	int gaussA = atoi(argv[6]);
	int gaussB = atoi(argv[7]);
	int gaussC = atoi(argv[8]);
	int gaussD = atoi(argv[9]);
	
	
	optimum.push_back(0);
	initialPosition.push_back(2);
	for(int i=1; i<d; i++){
		optimum.push_back(0);
		initialPosition.push_back(0);
	}
	
	randomv		  myrand;
	randomv		 &myrandRef = myrand;
	modelFunctions   myModel(myrandRef);
	modelFunctions  &myModelRef = myModel;
	environment	  myEnvironment(optimum);
	environment	 &envRef = myEnvironment;
	
	myModelRef.setFitnessFunction(wFuncType);
	myModelRef.setA(gaussA);
	myModelRef.setB(gaussB);
	myModelRef.setC(gaussC);
	myModelRef.setD(gaussD);
	myModelRef.setNumDimensions(d);
	
	string line;
	ifstream pathFile(pathFileName);
	if (pathFile.is_open()){
	  while (! pathFile.eof() ){
		if(mutVectors.size()>=pathLength){
			break;
		}
		getline (pathFile,line);
		//cout<<"current line: "<<line<<endl;
		vector<string> tokens;
		split( tokens, line, is_any_of(" "), token_compress_on );
		//cout<<"tokens: "<<joinArray(tokens,'\t')<<endl;
		//cout<<"size: "<<tokens.size()<<"\t"<<tokens[tokens.size()-1]<<"\t"<<atof(tokens[tokens.size()-2].c_str())<<endl;
		if(tokens.size() < 3 || atof(tokens[tokens.size()-2].c_str())!=1){
			continue;
		}
		
		vector<double> mutVec;
		for(int i=3; i<3+d; i++){
			mutVec.push_back(atof(tokens[i].c_str()));
		}
		//cout<<"mutVec: "<<joinArray(mutVec,'\t')<<endl;
		mutVectors.push_back(mutVec);
	  }
	}
	
	if(mutVectors.size()!=pathLength){
		//cout << mutVectors.size() << " mutations in path but looking for path length of "<<pathLength<<endl;
		exit(1);
	}
	vector<int> oneAlleles, twoAlleles, alleleMutated;
	vector<string> observedPoints, meanPopulationPhenotype;
	vector<double> observedMuts, currentFitness;
	currentFitness.push_back(envRef.fW(initialPosition,myModelRef));
	observedPoints.push_back(joinArray(initialPosition,'_'));
	
	effectiveNumPaths=0;
	loopCount=0;
	numPaths=0;
	validPaths=0;
	sumProb=0;
	
	//cout<<"starting to compute het weinreich paths!"<<endl;
	computeHeterozygousWeinreichPaths(oneAlleles, twoAlleles, initialPosition, 1.0, initialPosition, 0.0, observedPoints, observedMuts,currentFitness, 1, 1,0, alleleMutated, "", meanPopulationPhenotype, myModelRef, envRef);
	
	for(int i=0; i<probArray.size(); i++){
		effectiveNumPaths+=pow((probArray[i]/sumProb),2);
	}

	if(effectiveNumPaths!=0){
		effectiveNumPaths = 1.0/effectiveNumPaths;
	}
	cout<<"Prob array:\t"<<joinArray(probArray,'_')<<endl;
	cout<<"Prob sum:\t"<<sumProb<<endl;
	cout<<"Effective number of paths:\t"<<effectiveNumPaths<<endl;
	
	cout<<"Number Valid End States:\t"<<validEndStates.size()<<endl;
	
	_FPU_SETCW(fpu_oldcw);

}