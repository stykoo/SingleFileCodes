#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <thread>
#include <iomanip>
#include <random>
#include "utils.h"
#include "state.h"
#include "simul.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	std::vector<std::thread> threads;

	CoeffsCustom coeffs;
	int status = loadCoeffsCustom(p, coeffs);
	if (status) {
		return status;
	}

	std::vector< ObservablesVecCustom >
		allSumsObs(p.nbThreads, ObservablesVecCustom(p.nbSteps, coeffs));

	long nbSimulsPerThread = p.nbSimuls / p.nbThreads;
	
	if (p.nbSimuls % p.nbThreads != 0) {
		std::cerr << "Warning: nbSimuls is not a multiple of nbThreads."
			<< nbSimulsPerThread * p.nbThreads << " simulations will be done."
			<< std::endl;
	}

	// Threads
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
	        nbSimulsPerThread, coeffs, std::ref(allSumsObs[i]), rd())); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	// Initialize the total sum
	ObservablesVecCustom sumObs(p.nbSteps, coeffs);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		sumObs.add(allSumsObs[k]);
	}

	status = exportObservables(sumObs, coeffs, p);

	return status;
}

// Run several simulation and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
						   const CoeffsCustom &coeffs,
		                   ObservablesVecCustom &sumObs,
						   const unsigned int seed) {
    // Random generator
	//std::mt19937 rndGen(seed);
	VSLStreamStatePtr stream;
	vslNewStream(&stream, CUSTOM_RNG, seed);

	ObservablesVecCustom obs(p.nbSteps, coeffs);

	for (long s = 0 ; s < nbSimuls ; ++s) {
		if (p.verbose || p.visu) {
			std::cout << "Simulation " << s + 1 
				<< " on thread " << std::this_thread::get_id() << std::endl;
		}

		// Run a single simulation
		//runOneSimulation(p, obs, rndGen);
		runOneSimulation(p, obs, stream);

		// Add the observables to the sum
		sumObs.add(obs);

		if (p.visu) {
			std::cout << std::endl;
		}
	}

	vslDeleteStream(&stream);
}

/*
// Run a single simulation and compute the moments.
void runOneSimulation(const Parameters &p, ObservablesVecCustom &obs,
					  std::mt19937 &rndGen) {
	// To generate the time
	double t = 0., tLast = 0.;
	long tDiscrete = 0;
	std::exponential_distribution<double> rndTime((double) p.nbParticles);
	std::uniform_int_distribution<long> distInt(0, p.nbParticles - 1);
	std::uniform_real_distribution<double> distReal(0.0, 1.0);

	// Positions of the tracers
	State state;
	state.init(p, rndGen);
	const State initialState = state;

	// Compute the initial observables
	obs[tDiscrete].fromState(state, initialState, p.nbSites);
	++tDiscrete;

	if (p.visu) {
		state.visualize(p, 0.);
	}

	// Loop over time
	while (tDiscrete < p.nbSteps) {
		// Evolve
		state.update(p, distInt(rndGen), distReal(rndGen));
		t += rndTime(rndGen);

		while (t - tLast > p.dt && tDiscrete < p.nbSteps) {
			obs[tDiscrete].fromState(state, initialState, p.nbSites);
			++tDiscrete;
			tLast += p.dt;
		}
		if (p.visu) {
			state.visualize(p, t);
		}
	}
}
*/

// Run a single simulation and compute the moments.
void runOneSimulation(const Parameters &p, ObservablesVecCustom &obs,
					  VSLStreamStatePtr stream) {
	// To generate the time
	double t = 0., tLast = 0.;
	long tDiscrete = 0;
	double scalet = 1.0 / p.nbParticles;
	double times[BATCH_SIZE];
	int indices[BATCH_SIZE];
	double us[BATCH_SIZE];

	// Positions of the tracers
	State state;
	state.init(p, stream);
	const State initialState = state;

	// Compute the initial observables
	obs[tDiscrete].fromState(state, initialState, p.nbSites);
	++tDiscrete;

	if (p.visu) {
		state.visualize(p, 0.);
	}

	// Loop over time
	// Generation of the random numbers in batches
	while (tDiscrete < p.nbSteps) {
		vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream, BATCH_SIZE,
				         times, 0.0, scalet);
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, indices,
					 0, p.nbParticles);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, us,
				     0.0, 1.0);

		for (int i = 0 ; i < BATCH_SIZE ; ++i) {
			// Evolve
			state.update(p, indices[i], us[i]);
			t += times[i];

			while (t - tLast > p.dt && tDiscrete < p.nbSteps) {
				obs[tDiscrete].fromState(state, initialState, p.nbSites);
				++tDiscrete;
				tLast += p.dt;
			}
			if (p.visu) {
				state.visualize(p, t);
			}
		}
	}
}

// Load the custom coefficients for the observables
int loadCoeffsCustom(const Parameters &p, CoeffsCustom &coeffs) {
	coeffs.clear();

	std::ifstream file(p.custom);
	if (!file.is_open()) {
		std::cerr << "Could not open " << p.custom << std::endl;
		return 1;
	}

	for (std::string line; std::getline(file, line);) {
    	std::istringstream str(line);
    	coeffs.emplace_back(std::istream_iterator<int>(str),
			                std::istream_iterator<int>{});
	}

	for(const auto &c : coeffs) {
		if (c.size() != (size_t) p.nbTracers) {
			std::cerr << "Wrong size of coefficients in " << p.custom
				      << std::endl;
			return 1;
		}
	}
	return 0;
}

// Export the observables to a file.
int exportObservables(const ObservablesVecCustom &sumObs,
		              const CoeffsCustom &coeffs,
		              const Parameters &p) {
	std::ofstream file(p.output);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# SF_ContTime (" << __DATE__ <<  ", " << __TIME__ << "): ";
	p.print(file);
	file << "\n# t";
	for (size_t i = 0 ; i < coeffs.size() ; ++i) {
		file << " ";
		for (long j = 0 ; j < p.nbTracers ; ++j) {
			file << "x" << j+1 << "^" << coeffs[i][j];
		}
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
