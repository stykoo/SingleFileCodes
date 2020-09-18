#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
// #include <random>
#include <mkl.h>
#include <mkl_vsl.h>
#include "parameters.h"
#include "observables.h"

#define BATCH_SIZE 256
// See https://software.intel.com/en-us/mkl-developer-reference-c-brng-parameter-definition
#define CUSTOM_RNG VSL_BRNG_SFMT19937
#define MIN(a,b) (((a)<(b))?(a):(b))

struct Observables {
	std::vector< std::vector<long long> > moments;  // All the moments
	std::vector<long long> moments1TP;  // Moments of TP 1
	std::vector<long long> mom2;  // Second moment of each displacement

	long long occPos;  // Occupation of site 1
	long long occNeg;  // Occupation of site -1

	// Observables related to the occupations
	// 0 -> occupations (eta_l)
	// 1 -> occupation * occupation of site X+1 (eta_l * eta_1)
	// 2 -> X * eta_l
	// 3 -> X * eta_l * eta_1 
	std::vector< std::vector<long> > occObs;
};

int runSimulations(const Parameters &p);
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
						   const CoeffsCustom &coeffs,
		                   ObservablesVecCustom &sumObs,
						   const unsigned int seed);
/* void runOneSimulation(const Parameters &p, ObservablesVecCustom &obs,
   					  std::mt19937 &rndGen); */
void runOneSimulation(const Parameters &p, ObservablesVecCustom &obs,
   					  VSLStreamStatePtr stream);
int loadCoeffsCustom(const Parameters &p, CoeffsCustom &coeffs);
int exportObservables(const ObservablesVecCustom &sumObs,
		              const CoeffsCustom &coeffs,
		              const Parameters &p);

#endif
