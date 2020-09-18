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
 * Simuls1d.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Implementation of the simulation in a 1d system.
 * Different interactions are possible.
 */

#include <algorithm>
#include "Simuls1d.h"

// Generate an initial state.
int Simul1d::init(std::mt19937 &rndGen) {
	positions.resize(p.nbParticles);
	forces.resize(p.nbParticles);
	initXTracers.resize(p.nbTracers);

	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i] = p.length * (distribUnif(rndGen) - 0.5);
		forces[i] = 0;  // Arbitrary
	}

	// Sort by x value
	std::sort(positions.begin(), positions.end());

	// If the potential is strong, the order of the particles may not be
	// conserved in the first iterations: we do some thermalization.
	update(rndGen, true);
	for (int i = 0 ; i < MAX_ITERS_INIT ; ++i) {
		if (isOrdered()) {
			return 0;
		}
		std::sort(positions.begin(), positions.end());
		update(rndGen, true);
	}

	return 1;
}

#ifdef USE_MKL
// Generate an initial state.
int Simul1d::init(VSLStreamStatePtr stream) {
	positions.resize(p.nbParticles);
	forces.resize(p.nbParticles);
	aux.resize(p.nbParticles);
	initXTracers.resize(p.nbTracers);

	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, p.nbParticles,
			     positions.data(), -p.length/2, p.length/2);
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;  // Arbitrary
	}

	// Sort by x value
	std::sort(positions.begin(), positions.end());

	// If the potential is strong, the order of the particles may not be
	// conserved in the first iterations: we do some thermalization.
	update(stream, true);
	for (int i = 0 ; i < MAX_ITERS_INIT ; ++i) {
		if (isOrdered()) {
			return 0;
		}
		std::sort(positions.begin(), positions.end());
		update(stream, true);
	}

	return 1;
}
#endif

// Implement one step of the time evolution of the system.
void Simul1d::update(std::mt19937 &rndGen, const bool thermalization) {
	calcForcesBetweenParticles();
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		// Noise
		positions[i] += noise * distribNormal(rndGen);

		// Forces between particles
		positions[i] += p.timestep * forces[i];
	}

	// Forces on the tracers
	if (!thermalization) {
		for (long i=0 ; i<p.nbTracers ; ++i) {
			positions[p.idTracers[i]] += p.timestep * p.forces[i];
		}
	}
	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i] = periodicBC(positions[i], p.length);
	}
}

#ifdef USE_MKL
// Implement one step of the time evolution of the system.
void Simul1d::update(VSLStreamStatePtr stream, const bool thermalization) {
	calcForcesBetweenParticles();
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, p.nbParticles,
				  aux.data(), 0.0, noise);
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		// Noise + Forces between particles
		positions[i] += aux[i] + p.timestep * forces[i];
	}

	// Forces on the tracers
	if (!thermalization) {
		for (long i=0 ; i<p.nbTracers ; ++i) {
			positions[p.idTracers[i]] += p.timestep * p.forces[i];
		}
	}
	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i] = periodicBC(positions[i], p.length);
	}
}
#endif

// Get position in X of particle i
double Simul1d::getPosX(const long i) {
	return positions[i];
}

// Get profile
void Simul1d::getPosRel(std::vector<double> &posr, const long nbPts) {
	double x;
	long k;
	for (long j = 0 ; j < nbPts ; ++j) {
		posr[j] = 0.0;
	}
	for (long i = 1 ; i < p.nbParticles ; ++i) {
		x = periodicBCpos(positions[i] - positions[0], p.length);
		k = (long) (x / p.length * nbPts);
		posr[k] += 1.0;
	}
}

// Compute the forces between the particles.
// WE ASSUME THAT THE PARTICLES ARE ORDERED AND NEVER CROSS.
void SimulTonks::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		long iPrev = (i + p.nbParticles - 1) % p.nbParticles;

		double dx = periodicBC(positions[i] - positions[iPrev], p.length);
		if (dx < 1. && dx > 0.) {
			double f = p.eps * (1. - dx);
			forces[i] += f;
			forces[iPrev] -= f;
		}
	}
}

// Compute the forces between the particles.
void SimulCoulomb::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (long j=i+1 ; j<p.nbParticles ; ++j) {
			double dx = periodicBC(positions[i] - positions[j], p.length);
			if (dx != 0.) {
				double f = sign(dx) * p.eps / (dx * dx);
				forces[i] += f;
				forces[j] -= f;
			}
		}
	}
}

SimulCoulombCircle::SimulCoulombCircle(const Parameters &p) : Simul1d(p) {
	pref = std::sqrt(2) * M_PI * M_PI * p.eps / (p.length * p.length);
}

// Compute the forces between the particles.
void SimulCoulombCircle::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (long j=i+1 ; j<p.nbParticles ; ++j) {
			double t = periodicBC(positions[i] - positions[j], p.length)
				       / p.length;
			if (t != 0.) {
				double den = std::sqrt(1. - std::cos(t));
				double f = pref * std::sin(t) / (den * den);
				forces[i] += f;
				forces[j] -= f;
			}
		}
	}
}

// Compute the forces between the particles.
void SimulDipole::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (long j=i+1 ; j<p.nbParticles ; ++j) {
			double dx = periodicBC(positions[i] - positions[j], p.length);
			if (dx != 0.) {
				double f = sign(dx) * 3. * p.eps / mypow(dx, 4);
				forces[i] += f;
				forces[j] -= f;
			}
		}
	}
}

SimulDipoleCircle::SimulDipoleCircle(const Parameters &p) : Simul1d(p) {
	pref = 6. * M_SQRT2 * p.eps * mypow(M_PI, 4) / mypow(p.length, 4);
}

// Compute the forces between the particles.
void SimulDipoleCircle::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (long j=i+1 ; j<p.nbParticles ; ++j) {
			double t = periodicBC(positions[i] - positions[j], p.length);
			t *= 2. * M_PI / p.length;
			if (t != 0.) {
				double den = std::sqrt(1. - std::cos(t));
				double f = pref * std::sin(t) / mypow(den, 5);
				forces[i] += f;
				forces[j] -= f;
			}
		}
	}
}
