#include <cmath>
#include <map>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <thread>
#include "state.h"

State::State(const Parameters &_p, VSLStreamStatePtr stream) : p(_p) {
	positions.resize(p.nbParticles);
	displacements.assign(p.nbParticles, 0);
	occupations.assign(p.nbSites, 0);

	if (p.determ) {
		init_determ();
	} else {
		init(stream);
	}
}

void State::init(VSLStreamStatePtr stream) {
	// Assign the sites of the tracers
	positions[0] = 0;
	occupations[0] = 1;

	// Generate a random permutation of the remaining sites
	// and assign the first sites to the particles.
	int nbSitesNT = p.nbSites - 1;
	int nbPartsNT = p.nbParticles - 1;
	std::vector<long> seq(nbSitesNT);
	for (long i = 1 ; i < p.nbSites ; ++i) {
		seq[i-1] = i;
	}

	// Random numbers in [0,1)
	std::vector<double> ran(nbPartsNT);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, nbPartsNT, ran.data(),
			     0, 1);

	// Partial shuffling algorithm
	for (long i=0 ; i < nbPartsNT ; ++i) {
		// k is uniform in [i, nbSitesNT)
		long k = i + (long) (ran[i] * (double) (nbSitesNT - i));
		positions[i + 1] = seq[k]; // 1 tracer
		occupations[seq[k]] = 1;
		// No need to do the full swap
		seq[k] = seq[i];
	}
}

void State::init_determ() {
	// Density smaller or larger than 0.5
	if (2 * p.nbParticles <= p.nbSites) {
		const long step = p.nbSites / p.nbParticles;
		long u = 0;
		for (long i = 0 ; i < p.nbSites ; ++i) {
			if (i % step == 0) {
				positions[u++] = i;
				occupations[i] = 1;
			} else {
				occupations[i] = 0;
			}
		}
	} else {
		const long step = p.nbSites / (p.nbSites - p.nbParticles);
		long u = 0;
		for (long i = 0 ; i < p.nbSites ; ++i) {
			if (i % step == step / 2) {
				occupations[i] = 0;
			} else {
				positions[u++] = i;
				occupations[i] = 1;
			}
		}
	}
}

void State::update(long part, double u) {
	long pos = positions[part]; 

	long d = 2 * (u < p.proba) - 1;
	long pos_next = pos + d;
	if (pos_next == p.nbSites)
		pos_next = 0;
	else if (pos_next == -1)
		pos_next = p.nbSites - 1;

	// Exclusion rules
	int o = occupations[pos_next];
	positions[part] += (1 - o) * (pos_next - pos); // Move if not occupied
	displacements[part] += (1 - o) * d; // Add the displacement
	occupations[pos] = o; // Occupied only if doesn't move
	occupations[pos_next] = 1; // Next site is occupied anyways
}

double State::computeDispl() const {
	long w = 0;
	for (long i = 0 ; i < p.nbParticles ; ++i) {
		w += displacements[i];
	}
	return w / ((double) p.nbParticles);
}

double State::computeDisplSq() const {
	long w2 = 0;
	for (long i = 0 ; i < p.nbParticles ; ++i) {
		w2 += displacements[i] * displacements[i];
	}
	return w2 / ((double) p.nbParticles);
}

// Print the system into terminal.
void State::visualize(const double t) {
	std::map<long,char> mymap;
	std::map<long,char>::iterator it;
	for (long i = 0 ; i < p.nbParticles ; ++i) {
		mymap[positions[i]] = 'O';
	}

	std::cout << "\r";
	for (long i = 0 ; i < p.nbSites ; ++i) {
		it = mymap.find(i);
		if (it != mymap.end()) {
			std::cout << it->second;
		} else {
			std::cout << '.';
		}
	}
	std::cout << std::fixed << std::setprecision(5) <<  "  (t = " << t << ")"
		<< std::flush;
	
	// Sleep
	std::this_thread::sleep_for(std::chrono::milliseconds(p.sleep));
}
