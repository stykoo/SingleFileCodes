/*
Copyright Université Pierre et Marie Curie (2017)
Contributors: Alexis Poncet

alexis.poncet@ens.fr

This software is a computer program whose purpose is to simulate the
motion of interacting Brownian particles in a quasi-1d environment,
some of them being driven by an external force.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
/*
 * BrownianQuasi1d
 * simul.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Contain the implementation of the class Simulation and of the
 * generic functions related to the simulation.
 */

#include <algorithm>
#include <fstream>
#include <map>
#include <chrono>
#include <thread>
#include <iomanip>
#include "run.h"
#include "Simul.h"
#include "Simuls1d.h"
#include "SimulCanal.h"
#include "SimulPipe.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	std::vector<std::thread> threads;
	std::vector< std::vector<Observables> > allSumsObs(p.nbThreads);

	long nbSimulsPerThread = p.nbSimuls / p.nbThreads;
	
	if (p.nbSimuls % p.nbThreads != 0) {
		std::cerr << "Warning: nbSimuls is not a multiple of nbThreads."
			<< nbSimulsPerThread * p.nbThreads << " simulations will be done."
			<< std::endl;
	}

	// Threads
	bool failed = false;
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
			nbSimulsPerThread, std::ref(allSumsObs[i]), rd(), &failed)); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	if (failed) {
		return 2;
	}

	// Initialize the total sum
	std::vector<Observables> sumObs;
	initObservables(sumObs, p);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		addObservables(sumObs, allSumsObs[k], p);
	}

	int status;
	if (p.computeProfs) {
		status = exportMoments1(sumObs, p);
		if (status)
			return status;
		status = exportProfiles(sumObs, p);
	} else {
		status = exportObservables(sumObs, p);
	}
	return status;
}

// Run several simulations and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
						   std::vector<Observables> &sumObs,
						   const unsigned int seed, bool *failed) {
    // Random generator
#ifdef USE_MKL
	VSLStreamStatePtr rndGen;
	vslNewStream(&rndGen, VSL_BRNG_SFMT19937, seed);
#else
	std::mt19937 rndGen(seed);
#endif

	// Initialize sumObs
	initObservables(sumObs, p);

	for (long s = 0 ; s < nbSimuls ; ++s) {
		if (p.verbose) {
			std::cout << "Simulation " << s + 1 
				<< " on thread " << std::this_thread::get_id() << std::endl;
		}

		// Run a single simulation (of the right type!)
		Simulation *simul;
		if (p.simulName == "pipe") {
			simul = new SimulPipe(p);
		} else if (p.simulName == "tonks") {
			simul = new SimulTonks(p);
		} else if (p.simulName == "canal") {
			simul = new SimulCanal(p);
		} else if (p.simulName == "coulomb") {
			simul = new SimulCoulomb(p);
		} else if (p.simulName == "dipole") {
			simul = new SimulDipole(p);
		} else if (p.simulName == "coulCircle") {
			simul = new SimulCoulombCircle(p);
		} else if (p.simulName == "dipCircle") {
			simul = new SimulDipoleCircle(p);
		} else {
			return;
		}
		std::vector<Observables> obs;
		int status = simul->run(obs, rndGen);

		while (status == 1 || status == 3) {
			if (status == 1) {
				std::cerr << "Warning: Wrong order of the particles ";
			} else {
				std::cerr << "Warning: A NaN appeared in the observables ";
			}
			std::cerr << "(thread: " << std::this_thread::get_id()
				<< ", simulation: " << s+1 << "). "
				<< "Running simulation again." << std::endl;
			status = simul->run(obs, rndGen);
		}
		delete simul;

		if (status == 2) {
			std::cerr << "Could not initialize the system properly "
				<< "(thread: " << std::this_thread::get_id()
				<< ", simulation: " << s+1 << "). "
				<< "You should reduce the timestep." << std::endl;
			*failed = true;
			return;
		}

		// Add the observables to the sum
		addObservables(sumObs, obs, p);
	}
}
