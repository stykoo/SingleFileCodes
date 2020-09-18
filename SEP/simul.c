#include <stdio.h>
#include <math.h>
#include "simul.h"

void run(const long n_sites, const long n_parts, const double prob,
		 const long n_simuls, const int n_times, const double *tfins, 
		 long **observables, const unsigned long seed,
		 const int determ) {
	const int n_obs_tot = N_OBS * n_times;
	for (int k = 0 ; k < n_obs_tot ; ++k) {
#ifdef ONLY_CUMS
		observables[k][0] = 0;
#else
		for (long i = 0 ; i < n_sites ; ++i) {
			observables[k][i] = 0;
		}
#endif
	}

	VSLStreamStatePtr stream;
	vslNewStream(&stream, CUSTOM_RNG, seed);

	State state;
	state.positions = malloc(n_parts * sizeof(long));
	state.occupations = malloc(n_sites * sizeof(int));
	long *aux_seq;
	double *aux_ran;
	if (!determ) {
		aux_seq = malloc((n_sites - 1) * sizeof(long));
		aux_ran = malloc((n_parts - 1) * sizeof(double));
	}

	for (long i = 0 ; i < n_simuls ; ++i) {
		// Deterministic / random initialization
		if (determ) {
			init_determ(&state, n_sites, n_parts);
		} else {
			init_random(&state, n_sites, n_parts, stream, aux_seq, aux_ran);
		}
		double ti = 0.0, tf = 0.0;
		for (int k = 0 ; k < n_times ; ++k) {
			tf = tfins[k];
			ti = evolve(&state, n_sites, n_parts, prob, ti, tf, stream);
			// Shift observables depending on index of final time
			updateObs(&state, n_sites, &observables[k * N_OBS]);
		}
	}

	if (!determ) {
		free(aux_seq);
		free(aux_ran);
	}
	free(state.positions);
	free(state.occupations);

	vslDeleteStream(&stream);
}

void init_random(State *state, const long n_sites, const long n_parts,
		         VSLStreamStatePtr stream, long *aux_seq, double *aux_ran) {
	for (long i = 0 ; i < n_sites ; ++i) {
		state->occupations[i] = 0;
	}

	// The TP
	state->positions[0] = 0;
	state->occupations[0] = 1;
	state->winding = 0;

	// 1, 2, ..., n_sites-1
	for (long i = 0 ; i < n_sites - 1 ; ++i) {
		aux_seq[i] = i + 1;
	}

	// Random numbers in [0,1)
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts-1, aux_ran, 0, 1);

	// Partial shuffling
	for (long i = 0 ; i < n_parts - 1 ; ++i) {
		// k is uniform in [i, n_parts-1)
		long k = i + (long) (aux_ran[i] * (double) (n_sites - 1 - i));
		state->positions[i+1] = aux_seq[k];
		state->occupations[aux_seq[k]] = 1;
		// No need to do the full swap
		aux_seq[k] = aux_seq[i];
	}
}

void init_determ(State *state, const long n_sites, const long n_parts) {
	// Different strategies depending whether the system is more or less than
	// half full.
	if (n_sites >= 2 * n_parts) {
		for (long i = 0 ; i < n_sites ; ++i) {
			state->occupations[i] = 0;
		}

		// Place the particles equidistantly
		const long step = n_sites / n_parts; // At least 2
		for (long i = 0 ; i < n_parts ; ++i) {
			long s = i * step;
			state->positions[i] = s;
			state->occupations[s] = 1;
		}

	} else {
		for (long i = 0 ; i < n_sites ; ++i) {
			state->occupations[i] = 1;
		}

		const long n_vacs = n_sites - n_parts;
		const long step = n_sites / n_vacs; // At least 2
		const long hstep = step / 2;
		for (long i = 0 ; i < n_vacs ; ++i) {
			long s = i * step + hstep;
			state->occupations[s] = 0;
		}

		long j = 0;
		for (long i = 0 ; i < n_sites ; ++i) {
			if (state->occupations[i])
				state->positions[j++] = i;
		}
	}
	state->winding = 0;
}

double evolve(State *state, const long n_sites, const long n_parts,
		      const double prob, const double tini, const double tfin,
			  VSLStreamStatePtr stream) {
	const int unbiased = (prob == 0.5); // 0.5 should be represented exactly
	const double biais = prob - 0.5;
	const double lambda = (tfin - tini) * n_parts;
	const long endsite = n_sites - 1;
	int indices[BATCH_SIZE];
	int bern[BATCH_SIZE];
	double us[BATCH_SIZE];

	int n_iters = 0; // Number of time iterations (Poisson distributed)
	viRngPoisson(VSL_RNG_METHOD_POISSON_PTPE, stream, 1, &n_iters, lambda);

	while (n_iters) {
		int n_i = MIN(n_iters, BATCH_SIZE);

		// Generation of the random numbers in batches
		if (unbiased) {
			viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, n_i,
					       bern, 0.5);
		} else {
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_i,
					     us, 0.0, 1.0);
		}
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_i, indices,
					 0, n_parts);
		
		for (int i = 0 ; i < n_i ; ++i) {
			const long part = indices[i];
			const long pos = state->positions[part];

			// Next position
			// Some profiling was done on this
			long pos_next;
			if (unbiased) {
				pos_next = pos + 2 * bern[i] - 1;
			} else {
				double pp = 0.5 + (part == 0) * biais;
				pos_next = pos + 2 * (us[i] < pp) - 1;
			}
			if (pos_next == n_sites)
				pos_next = 0;
			else if (pos_next == -1)
				pos_next = endsite;

			int o = state->occupations[pos_next];
			long dpos = (1 - o) * (pos_next - pos);
			state->positions[part] += dpos; // Move if next site not occupied
			state->occupations[pos] = o; // Occupied only if doesn't move
			state->occupations[pos_next] = 1; // Next site is occupied anyways

			// Add or substract one if the TP did one turn
			long dw = (part == 0) * (1 - o) * (
			           (pos == endsite && pos_next == 0)
					   - (pos == 0 && pos_next == endsite) 
					  );
			state->winding += dw;
			/*if (dw) {
				printf("%ld: %ld, %ld, %ld\n", dw, part, pos, pos_next);
			}*/
		} 

		n_iters -= n_i;
	}
	// This may introduce a small error
	return tfin;
}

/*
double evolve(State *state, const long n_sites, const long n_parts,
	          const double prob, const double tini, const double tfin,
		      VSLStreamStatePtr stream) {
	const double biais = prob - 0.5;
	const double rate = 1.0 / n_parts;
	int indices[BATCH_SIZE];
	double times[BATCH_SIZE], us[BATCH_SIZE];

	double t = tini;

	while (1) {
		// Generation of the random numbers in batches
		vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream, BATCH_SIZE,
					     times, 0, rate);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, us,
					 0.0, 1.0);
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, indices,
					 0, n_parts);
		
		for (long i = 0 ; i < BATCH_SIZE ; ++i) {
			t += times[i];
			if (t >= tfin) {
				return t;
			}

			const long part = indices[i];
			const long pos = state->positions[part];

			// Next position
			// This was found to be the fastest by profiling
			double pp = 0.5 + (part == 0) * biais;
			long pos_next = pos + 2 * (us[i] < pp) - 1;
			if (pos_next == n_sites)
				pos_next = 0;
			else if (pos_next == -1)
				pos_next = n_sites - 1;

			int o = state->occupations[pos_next];
			long dpos = (1 - o) * (pos_next - pos);
			state->positions[part] += dpos; // Move if next site not occupied
			state->occupations[pos] = o; // Occupied only if doesn't move
			state->occupations[pos_next] = 1; // Next site is occupied anyways
		} 
	}
}
*/

void updateObs(State *state, const long n_sites, long **observables) {
	long X = state->positions[0];
	// long n_sites_2 = n_sites / 2;
	// long xPer = (X + n_sites_2) % n_sites - n_sites_2;
	long xPer = X + (state->winding * n_sites);

#ifdef ONLY_CUMS
	long xPer2 = xPer * xPer;
	long xPer4 = xPer2 * xPer2;
	observables[0][0] += xPer;
	observables[1][0] += xPer2;
	observables[2][0] += xPer2 * xPer;
	observables[3][0] += xPer4;
	observables[4][0] += xPer4 * xPer;
	observables[5][0] += xPer4 * xPer2;
#else
	long ep = state->occupations[(X + 1) % n_sites];
	long em = state->occupations[(X + n_sites - 1) % n_sites];
	long xPerEp = xPer * ep;
	long xPerEm = xPer * em;

  #ifdef LARGE_NOBS
	long xPer2 = xPer * xPer;
	long xPer2Ep = xPer2 * ep;
	long xPer2Em = xPer2 * em;
	long xPer3 = xPer2 * xPer;
	long xPer3Ep = xPer3 * ep;
	long xPer3Em = xPer3 * em;
  #endif

	#pragma ivdep // Go icc, you can vectorize it (no dependancy)
	for (long i = 0 ; i < n_sites ; ++i) {
		long k = (X + i) % n_sites;
		long er = state->occupations[k];

		observables[0][i] += er;
		observables[1][i] += ep * er;
		observables[2][i] += em * er;
		observables[3][i] += xPer * er;
		observables[4][i] += xPerEp * er;
		observables[5][i] += xPerEm * er;
		observables[6][i] += xPer;
  #ifdef LARGE_NOBS
		observables[7][i] += xPer2 * er;
		observables[8][i] += xPer2Ep * er;
		observables[9][i] += xPer2Em * er;
		observables[10][i] += xPer2;
		observables[11][i] += xPer3 * er;
		observables[12][i] += xPer3Ep * er;
		observables[13][i] += xPer3Em * er;
		observables[14][i] += xPer3;
  #endif
	}
#endif
}
