#include <map>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <thread>
#include "state.h"

void State::init(const Parameters &p, VSLStreamStatePtr stream) {
	x_tp = 0;
	positions.resize(p.nbParticles);
	occupations.assign(p.nbSites, 0);

	// Assign the sites of the tracer
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

void State::init_determ(const Parameters &p) {
	positions.resize(p.nbParticles);
	occupations.assign(p.nbSites, 0);

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

void State::update(const Parameters &p, long part, double u) {
	long pos = positions[part]; 

	// Probability to jump to the right
	double pr = (part == 0) ? p.prob : DEFAULT_PROBA_RIGHT;

	// New version
	long dx = 2 * (u < pr) - 1;
	long pos_next = pos + dx;
	if (pos_next == p.nbSites)
		pos_next = 0;
	else if (pos_next == -1)
		pos_next = p.nbSites - 1;

	int o = occupations[pos_next];
	long dpos = (1 - o) * (pos_next - pos);
	positions[part] += dpos; // Move if next site not occupied
	occupations[pos] = o; // Occupied only if doesn't move
	occupations[pos_next] = 1; // Next site is occupied anyways
	
	if (part == 0) {
		x_tp += (1 - o) * dx;
	}
}

// Print the system into terminal.
void State::visualize(const Parameters &p, const double t) {
	std::map<long,char> mymap;
	std::map<long,char>::iterator it;
	mymap[positions[0]] = 'X';
	for (long i = 1 ; i < p.nbParticles ; ++i) {
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
