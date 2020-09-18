#ifndef SIMUL_H
#define SIMUL_H

#include "mkl.h"
#include "mkl_vsl.h"

#ifdef ONLY_CUMS
	#define N_OBS 6
#elif LARGE_NOBS
	#define N_OBS 15
#else
	#define N_OBS 7
#endif

#define BATCH_SIZE 256
// See https://software.intel.com/en-us/mkl-developer-reference-c-brng-parameter-definition
#define CUSTOM_RNG VSL_BRNG_SFMT19937

#define MIN(a,b) (((a)<(b))?(a):(b))

typedef struct {
	long *positions;
	int *occupations;
	long winding;
} State;

/*void run(const long n_sites, const long n_parts, const double prob,
		 const double tfin, const long n_simuls, const long n_times,
		 long **observables, const unsigned long seed, const int determ);*/
void run(const long n_sites, const long n_parts, const double prob,
		 const long n_simuls, const int n_times, const double *tfins, 
		 long **observables, const unsigned long seed, const int determ);

void init_random(State *state, const long n_sites, const long n_parts,
		         VSLStreamStatePtr stream, long *aux_seq, double *aux_ran);
void init_determ(State *state, const long n_sites, const long n_parts);
double evolve(State *state, const long n_sites, const long n_parts,
		      const double prob, const double tini, const double tfin,
			  VSLStreamStatePtr stream);
void updateObs(State *state, const long n_sites, long **observables);

#endif
