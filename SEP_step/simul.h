#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
#include <mkl.h>
#include <mkl_vsl.h>
#include "parameters.h"
#include "observables.h"

#define BATCH_SIZE 256
// See https://software.intel.com/en-us/mkl-developer-reference-c-brng-parameter-definition
#define CUSTOM_RNG VSL_BRNG_SFMT19937
#define MIN(a,b) (((a)<(b))?(a):(b))

int runSimulations(const Parameters &p);
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   ObservablesVec &sumObs,
						   const unsigned int seed);
void runOneSimulation(const Parameters &p, ObservablesVec &obs,
   					  VSLStreamStatePtr stream);
int exportObservables(const ObservablesVec &sumObs, const Parameters &p);
int exportProfs(const ObservablesVec &sumObs, const Parameters &p);

#endif
