#include <map>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <thread>
#include "utils.h"
#include "state.h"

/*
void State::init(const Parameters &p, std::mt19937 &rndGen) {
	positions.resize(p.nbParticles);
	occupations.assign(p.nbSites, -1);

	// Assign the sites of the tracers
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		positions[i] = p.initPos[i];
		occupations[p.initPos[i]] = i;
	}

	// Generate a random permutation of the remaining sites
	// and assign the first sites to the particles.
	std::vector<long> seq(p.nbSites - p.nbTracers);
	long i = 0, j = 0;
	for (long k = 0 ; k < p.nbTracers ; ++k) {
		while (i < p.initPos[k]) {
			seq[j++] = i++;
		}
		++i;
	}
	while (i < p.nbSites) {
		seq[j++] = i++;
	}

	std::shuffle(seq.begin(), seq.end(), rndGen);
	for (long i=p.nbTracers ; i < p.nbParticles ; ++i) {
		positions[i] = seq[i - p.nbTracers];
		occupations[seq[i - p.nbTracers]] = i;
	}
}
*/

void State::init(const Parameters &p, VSLStreamStatePtr stream) {
	positions.resize(p.nbParticles);
	occupations.assign(p.nbSites, 0);

	// Assign the sites of the tracers
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		positions[i] = p.initPos[i];
		occupations[p.initPos[i]] = 1;
	}

	// Generate a random permutation of the remaining sites
	// and assign the first sites to the particles.
	int nbSitesNT = p.nbSites - p.nbTracers;
	int nbPartsNT = p.nbParticles - p.nbTracers;
	std::vector<long> seq(nbSitesNT);
	long i = 0, j = 0;
	for (long k = 0 ; k < p.nbTracers ; ++k) {
		while (i < p.initPos[k]) {
			seq[j++] = i++;
		}
		++i;
	}
	while (i < p.nbSites) {
		seq[j++] = i++;
	}

	// Random numbers in [0,1)
	std::vector<double> ran(nbPartsNT);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, nbPartsNT, ran.data(),
			     0, 1);

	// Partial shuffling algorithm
	for (long i=0 ; i < nbPartsNT ; ++i) {
		// k is uniform in [i, nbSitesNT)
		long k = i + (long) (ran[i] * (double) (nbSitesNT - i));
		positions[i + p.nbTracers] = seq[k];
		occupations[seq[k]] = 1;
		// No need to do the full swap
		seq[k] = seq[i];
	}
}

void State::update(const Parameters &p, long part, double u) {
	long pos = positions[part]; 

	// Probability to jump to the right
	double pr = (part < p.nbTracers) ? p.probas[part] : DEFAULT_PROBA_RIGHT;

	// New version
	long pos_next = pos + 2 * (u < pr) - 1;
	if (pos_next == p.nbSites)
		pos_next = 0;
	else if (pos_next == -1)
		pos_next = p.nbSites - 1;

	int o = occupations[pos_next];
	long dpos = (1 - o) * (pos_next - pos);
	positions[part] += dpos; // Move if next site not occupied
	occupations[pos] = o; // Occupied only if doesn't move
	occupations[pos_next] = 1; // Next site is occupied anyways

	/*if (!occupations[pos_next]) {
		positions[part] = pos_next;
		occupations[pos] = 0;
		occupations[pos_next] = 1;
	}*/

	/*
	// Old version
	if (u < pr) {
		long posR = periodicAdd1(pos, p.nbSites);
		if(occupations[posR] == -1){
			occupations[pos] = -1;
			occupations[posR] = part;
			positions[part] = posR;
		}
	} else {
		long posL = periodicSubs1(pos, p.nbSites);
		if(occupations[posL] == -1){
			occupations[pos] = -1;
			occupations[posL] = part;
			positions[part] = posL;
		}
	}
	*/
}

// Displacements of the n first particles with respect to an initial state
/* void State::displ(std::vector<long> &dxs, const State &s_ini, const long n,
		          const long len) const {
	for (long i = 0 ; i < n ; ++i) {
		dxs[i] = periodicBCsym(positions[i] - s_ini.positions[i], len);
	}
} */

// Print the system into terminal.
void State::visualize(const Parameters &p, const long t) {
	std::map<long,char> mymap;
	std::map<long,char>::iterator it;
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		mymap[positions[i]] = 'X';
	}
	for (long i = p.nbTracers ; i < p.nbParticles ; ++i) {
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
