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
#include "observables.h"

// Initialize a vector of observables
void initObservables(std::vector<Observables> &obs, const Parameters &p) {
	const long n = p.nbIters / p.skip;
	obs.resize(n);
	for (long t = 0 ; t < n ; ++t) {
		if (p.computeProfs) {
			obs[t].moments1.assign(DEFAULT_NB_MOMENTS, 0.);
			obs[t].profiles.assign(DEFAULT_NB_MOMENTS,
					               std::vector<double>(p.nbPtsProfs, 0.));
		} else {
			obs[t].pos.assign(p.nbTracers, 0.);
			obs[t].displ.assign(p.nbTracers, 0.);
		}
	}
}

// Add observables o2 to observables o1.
void addObservables(std::vector<Observables> &obs1,
					const std::vector<Observables> &obs2, const Parameters &p)
{
	const long n = p.nbIters / p.skip;
	for (long t = 0 ; t < n ; ++t) {
		if (p.computeProfs) {
			for (long i = 0 ; i < DEFAULT_NB_MOMENTS ; ++i) {
				obs1[t].moments1[i] += obs2[t].moments1[i];
				for (long j = 0 ; j < p.nbPtsProfs ; ++j) {
					obs1[t].profiles[i][j] += obs2[t].profiles[i][j];
				}
			}
		} else {
			for (long i = 0 ; i < p.nbTracers ; ++i) {
				obs1[t].pos[i] += obs2[t].pos[i];
				obs1[t].displ[i] += obs2[t].displ[i];
			}
		}
	}
}

// Export the observables to a file.
int exportObservables(const std::vector<Observables> &sumObs,
					  const Parameters &p) {
	std::ofstream file(p.output);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# BrownianQuasi1d  (" << __DATE__ <<  ", " << __TIME__
		<< "): ";
	printParameters(p, file);
	file << "\n# t";
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		file << " ";
		file << "x" << i+1;
	}
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		file << " ";
		file << "dx" << i+1;
	}
	file << "\n";

	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

	const long n = p.nbIters / p.skip;
	// Data (we write the average and not the sum)
	for (long k = 0 ; k < n ; ++k) {
		file << k * p.skip * p.timestep;
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			file << " " << sumObs[k].pos[i] / p.nbSimuls;
		}
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			file << " " << sumObs[k].displ[i] / p.nbSimuls;
		}
		file << "\n";
	}

	file.close();

	return 0;
}

// Export the moments of TP 1
int exportMoments1(const std::vector<Observables> &sumObs,
				   const Parameters &p) {
	std::string fname = p.output + "_moments.dat";

	std::ofstream file(fname);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# BrownianQuasi1d  (" << __DATE__ <<  ", " << __TIME__
		<< "): ";
	printParameters(p, file);
	file << "\n# t";
	for (size_t i = 0 ; i < DEFAULT_NB_MOMENTS ; ++i) {
		file << " X^" << i+1;
	}
	file << "\n";

	const long n = p.nbIters / p.skip;

	for (long i = 0 ; i < n ; ++i) {
		double t = i * p.skip * p.timestep;
		file << t;
		for (size_t k = 0 ; k < DEFAULT_NB_MOMENTS  ; ++k) {
			file << " " << sumObs[i].moments1[k] / p.nbSimuls;
		}
		file << "\n";
	}
	file.close();
	return 0;
}

// Export the profiles to files.
int exportProfiles(const std::vector<Observables> &sumObs,
				   const Parameters &p) {
	const long n = p.nbIters / p.skip;

	for (long i = 0 ; i < n ; ++i) {
		double t = i * p.skip * p.timestep;
		std::string fname = p.output + "_prof" + std::to_string(t) + ".dat";

		std::ofstream file(fname);
		if (!file.is_open()) {
			return 1;
		}

		// Header
		file << "# BrownianQuasi1d  (" << __DATE__ <<  ", " << __TIME__
			<< "): ";
		printParameters(p, file);
		file << "\n# x";
		for (size_t i = 0 ; i < DEFAULT_NB_MOMENTS ; ++i) {
			file << " rho(x)*X^" << i;
		}
		file << "\n";
		for (long j = 0 ; j < p.nbPtsProfs ; ++j) {
			file << j * p.precProfs;
			for (size_t k = 0 ; k < DEFAULT_NB_MOMENTS  ; ++k) {
				file << " "
					 << (sumObs[i].profiles[k][j] * p.nbPtsProfs
						 / p.length / p.nbSimuls);
			}
		  	file << "\n";
		}
		file.close();
	}
	return 0;
}
