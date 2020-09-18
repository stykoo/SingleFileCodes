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
#include "Simul.h"

// Load the parameters and initialize the distributions
Simulation::Simulation(const Parameters &params) : p(params),
	noise(std::sqrt(2. * p.temperature * p.timestep)),
	distribUnif(0., 1.), distribNormal(0., 1.) {
		initXTracers.assign(p.nbTracers, 0.);
}

// Run the simulation.
// Return 1 if particles cross.
// Return 2 if initialization failed.
// Return 3 if a NaN appears in the observables
int Simulation::run(std::vector<Observables> &obs, std::mt19937 &rndGen) {
	// Initialization of the positions
	int status = init(rndGen);
	if (status) {
		return 2;
	}

	// Thermalization
	for (long j=0 ; j<p.nbItersTh ; ++j) {
		update(rndGen, true);
		if (p.checkOrder && !isOrdered()) {
			return 1;
		}
	}

	// Initial observables
	setInitXTracers();
	initObservables(obs, p);
	computeObservables(obs[0]);

	// Main loop
	for (long j=0 ; j<p.nbIters-1 ; ++j) {
		update(rndGen, false);
		if ((j + 1) % p.skip == 0) {
			int status = computeObservables(obs[(j+1)/p.skip]);
			if (status != 0) {
				return 3;
			}

			if (p.checkOrder && !isOrdered()) {
				return 1;
			}
		}
	}
	return 0;
}

// Run the simulation.
// Return 1 if particles cross.
// Return 2 if initialization failed.
// Return 3 if a NaN appears in the observables
#ifdef USE_MKL
int Simulation::run(std::vector<Observables> &obs, VSLStreamStatePtr stream) {
	// Initialization of the positions
	int status = init(stream);
	if (status) {
		return 2;
	}

	// Thermalization
	for (long j=0 ; j<p.nbItersTh ; ++j) {
		update(stream, true);
		if (p.checkOrder && !isOrdered()) {
			return 1;
		}
	}

	// Initial observables
	setInitXTracers();
	initObservables(obs, p);
	computeObservables(obs[0]);

	// Main loop
	for (long j=0 ; j<p.nbIters-1 ; ++j) {
		update(stream, false);
		if ((j + 1) % p.skip == 0) {
			int status = computeObservables(obs[(j+1)/p.skip]);
			if (status != 0) {
				return 3;
			}

			if (p.checkOrder && !isOrdered()) {
				return 1;
			}
		}
	}
	return 0;
}
#endif

// Set the initial positions of the tracers
void Simulation::setInitXTracers() {
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		initXTracers[i] = getPosX(p.idTracers[i]);
	}
}

// Compute the observables.
// Return 1 if a NaN appears, 0 otherwise.
int Simulation::computeObservables(Observables &o) {
	if (p.computeProfs) {
		double dx = periodicBC(getPosX(p.idTracers[0]) - initXTracers[0],
				               p.length);
		o.moments1[0] = dx;
		getPosRel(o.profiles[0], p.nbPtsProfs);
		for (long i = 1 ; i < DEFAULT_NB_MOMENTS ; ++i) {
			o.moments1[i] = dx * o.moments1[i-1];
			for (long j = 0 ; j < p.nbPtsProfs ; ++j) {
				o.profiles[i][j] = dx * o.profiles[i-1][j];
			}
		}
	} else {
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			o.pos[i] = getPosX(p.idTracers[i]);
			o.displ[i] = periodicBC(getPosX(p.idTracers[i]) - initXTracers[i],
									p.length);
			if (std::isnan(o.pos[i]) || std::isnan(o.displ[i])) {
				return 1;
			}
		}
	}
	return 0;
}

// Check if the positions of the particles are ordered.
bool Simulation::isOrdered() {
	long c = 0;
	for (long i = 0 ; i < p.nbParticles-1 ; ++i) {
		if (getPosX(i) > getPosX(i+1)) {
			++c;
		}
	}
	if (getPosX(p.nbParticles-1) > getPosX(0)) {
		++c;
	}
	return (c < 2);
}

