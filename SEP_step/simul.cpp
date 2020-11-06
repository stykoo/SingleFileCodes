#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <thread>
#include <iomanip>
#include <random>
#include "state.h"
#include "simul.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	std::vector<std::thread> threads;

	std::vector< ObservablesVec>
		allSumsObs(p.nbThreads, ObservablesVec(p.nbSteps));

	long nbSimulsPerThread = p.nbSimuls / p.nbThreads;
	
	if (p.nbSimuls % p.nbThreads != 0) {
		std::cerr << "Warning: nbSimuls is not a multiple of nbThreads."
			<< nbSimulsPerThread * p.nbThreads << " simulations will be done."
			<< std::endl;
	}

	// Threads
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
	        nbSimulsPerThread, std::ref(allSumsObs[i]), rd())); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	// Initialize the total sum
	ObservablesVec sumObs(p.nbSteps);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		sumObs.add(allSumsObs[k]);
	}

	int status = exportObservables(sumObs, p);

	return status;
}

// Run several simulation and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   ObservablesVec&sumObs,
						   const unsigned int seed) {
    // Random generator
	VSLStreamStatePtr stream;
	vslNewStream(&stream, CUSTOM_RNG, seed);

	ObservablesVec obs(p.nbSteps);

	for (long s = 0 ; s < nbSimuls ; ++s) {
		if (p.verbose || p.visu) {
			std::cout << "Simulation " << s + 1 
				<< " on thread " << std::this_thread::get_id() << std::endl;
		}

		// Run a single simulation
		runOneSimulation(p, obs, stream);

		// Add the observables to the sum
		sumObs.add(obs);

		if (p.visu) {
			std::cout << std::endl;
		}
	}

	vslDeleteStream(&stream);
}

// Run a single simulation and compute the moments.
void runOneSimulation(const Parameters &p, ObservablesVec &obs,
					  VSLStreamStatePtr stream) {
	// To generate the time
	double t = 0., tLast = 0.;
	long tDiscrete = 0;
	//double scalet = 1.0 / p.nbParticles;
	//double times[BATCH_SIZE];
	//int indices[BATCH_SIZE];
	double dt;
	double us[BATCH_SIZE];

	// Positions of the tracers
	State state;
	state.init(p, stream);
	const State initialState = state;

	// Compute the initial observables
	obs[tDiscrete].fromState(state, initialState);
	++tDiscrete;

	if (p.visu) {
		state.visualize(p, 0.);
	}

	// Loop over time
	// Generation of the random numbers in batches
	while (tDiscrete < p.nbSteps) {
		//vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream, BATCH_SIZE,
		//		         times, 0.0, scalet);
		//viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, indices,
		//			 0, p.nbParticles);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, us,
				     0.0, 1.0);

		for (int i = 0 ; i < BATCH_SIZE ; ++i) {
			// Evolve
			dt = state.update(p, us[i], stream);
			t += dt;

			while (t - tLast > p.dt && tDiscrete < p.nbSteps) {
				obs[tDiscrete].fromState(state, initialState);
				++tDiscrete;
				tLast += p.dt;
			}
			if (p.visu) {
				state.visualize(p, t);
			}
		}
	}
}

// Export the observables to a file.
int exportObservables(const ObservablesVec &sumObs, const Parameters &p) {
	std::ofstream file(p.output);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# SF_ContTime (" << __DATE__ <<  ", " << __TIME__ << "): ";
	p.print(file);
	file << "\n# t";
	for (size_t i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
			file << " x^" << i+1;
	}
	file << "\n";

	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

	// Data (we write the average and not the sum)
	for (long i = 0 ; i < p.nbSteps ; ++i) {
		file << i * p.dt << " ";
		sumObs[i].print(p.nbSimuls, file);
		file << "\n";
	}

	file.close();
	return 0;
}
