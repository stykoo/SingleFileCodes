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
		allSumsObs(p.nbThreads,
				   ObservablesVec(p.nbSteps, p.computeProfs, p.nbSitesProf,
					              p.rev));

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
	ObservablesVec sumObs(p.nbSteps, p.computeProfs, p.nbSitesProf, p.rev);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		sumObs.add(allSumsObs[k]);
	}

	int status = exportObservables(sumObs, p);
	if (status) return status;
	if (p.computeProfs) {
		status = exportProfs(sumObs, p);
	}
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

	ObservablesVec obs(p.nbSteps, p.computeProfs, p.nbSitesProf, p.rev);

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
	double us[BATCH_SIZE], dts[BATCH_SIZE];
	int parts[BATCH_SIZE];
	long nbParts;

	// Initialization
	State state;
	if (p.determ) {
		state.init_determ(p);
	} else {
		state.init(p, stream);
	}
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
		nbParts = state.getNbParts();
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, us,
				     0.0, 1.0);
		// bath particles + TP + left / right reservoirs -> nbParts + 3
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, parts,
					 0, nbParts+3);
		vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream, BATCH_SIZE,
						 dts, 0.0, 1./(nbParts+3));

		for (int i = 0 ; i < BATCH_SIZE ; ++i) {
			// Evolve
			t += dts[i];
			state.update(p, us[i], parts[i]);

			while (t - tLast > p.dt && tDiscrete < p.nbSteps) {
				obs[tDiscrete].fromState(state, initialState);
				++tDiscrete;
				tLast += p.dt;
				if (tDiscrete >= p.nbSteps)
					goto endloop;
			}
			if (p.visu) {
				state.visualize(p, t);
			}
			// Regenerate the random numbers if number of particles changed
			if (state.getNbParts() != nbParts) {
				break;
			}
		}
	}
endloop:
	std::cout << "";
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

int exportProfs(const ObservablesVec &sumObs, const Parameters &p) {
	for (long i = 0 ; i < p.nbSteps ; ++i) {
		double t = i * p.dt;
		std::ofstream file(p.output + "_profP_" + std::to_string(t) + ".dat");
		if (!file.is_open()) {
			return 1;
		}
		file << "# SF_ContTime (" << __DATE__ <<  ", " << __TIME__ << "): ";
		p.print(file);
		if (p.rev) {
			file << "\n# r er (1-e+)*er (1-e-)*er";
			file << " X*er X*(1-e+)*er X*(1-e-)*er X";
			file << " X^2*er X^2*(1-e+)*er X^2*(1-e-)*er X^2";
			file << " X^3*er X^3*(1-e+)*er X^3*(1-e-)*er X^3\n";
		} else {
			file << "\n# r er e+*er e-*er";
			file << " X*er X*e+*er X*e-*er X";
			file << " X^2*er X^2*e+*er X^2*e-*er X^2";
			file << " X^3*er X^3*e+*er X^3*e-*er X^3\n";
		}
		file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);
		// Data (we write the average and not the sum)
		sumObs[i].printProfsP(p.nbSimuls, file);
		file << "\n";
		file.close();
	}
	for (long i = 0 ; i < p.nbSteps ; ++i) {
		double t = i * p.dt;
		std::ofstream file(p.output + "_profM_" + std::to_string(t) + ".dat");
		if (!file.is_open()) {
			return 1;
		}
		file << "# SF_ContTime (" << __DATE__ <<  ", " << __TIME__ << "): ";
		p.print(file);
		if (p.rev) {
			file << "\n# r er (1-e+)*er (1-e-)*er";
			file << " X*er X*(1-e+)*er X*(1-e-)*er X";
			file << " X^2*er X^2*(1-e+)*er X^2*(1-e-)*er X^2";
			file << " X^3*er X^3*(1-e+)*er X^3*(1-e-)*er X^3\n";
		} else {
			file << "\n# r er e+*er e-*er";
			file << " X*er X*e+*er X*e-*er X";
			file << " X^2*er X^2*e+*er X^2*e-*er X^2";
			file << " X^3*er X^3*e+*er X^3*e-*er X^3\n";
		}
		file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);
		// Data (we write the average and not the sum)
		sumObs[i].printProfsM(p.nbSimuls, file);
		file << "\n";
		file.close();
	}
	return 0;
}
