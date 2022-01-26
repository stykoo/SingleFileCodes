#include <map>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <thread>
#include <cassert>
#include "utils.h"
#include "state.h"

State::State(const Parameters &p) :
		nbSitesX(p.nbSitesX), nbSitesY(p.nbSitesY),
		nbSitesTot(p.nbSitesX * p.nbSitesY), nbParticles(p.nbParticles) {
}

void State::init(VSLStreamStatePtr stream) {
	positions.resize(nbParticles);
	occupations.assign(nbSitesTot, 0);

	// Assign the sites of the tracer
	positions[0] = 0;
	occupations[0] = 1;

	// Generate a random permutation of the remaining sites
	// and assign the first sites to the particles.
	int nbSitesNT = nbSitesTot - 1;
	int nbPartsNT = nbParticles - 1;
	std::vector<long> seq(nbSitesNT);
	for (long i = 1 ; i < nbSitesTot ; ++i) {
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
	reset();
}

void State::init_determ() {
	positions.resize(nbParticles);
	occupations.resize(nbSitesTot);

	if (nbSitesTot >= 2 * nbParticles) {
		assert(nbSitesTot % nbParticles == 0);
		// Place the particles equidistantly
		const long step = nbSitesTot / nbParticles;
		long u = 0;
		for (long j = 0 ; j < nbSitesY ; ++j) {
			for (long i = 0 ; i < nbSitesX ; ++i) {
				long k = j * nbSitesX + i;
				if ((i + j) % step == 0 && u < nbParticles) {
					positions[u++] = k;
					occupations[k] = 1;
				} else {
					occupations[k] = 0;
				}
			}
		}
	} else {
		long n_vacs = nbSitesTot - nbParticles;
		assert(nbSitesTot % n_vacs == 0);
		long step = nbSitesTot / n_vacs; // At least 2
		long hstep = step / 2;
		long u = 0;
		for (long j = 0 ; j < nbSitesY ; ++j) {
			for (long i = 0 ; i < nbSitesX ; ++i) {
				long k = j * nbSitesX + i;
				if ((i + j) % step != hstep && u < nbParticles) {
					positions[u++] = k;
					occupations[k] = 1;
				} else {
					occupations[k] = 0;
				}
			}
		}
	}
	reset();
}

void State::reset() {
	// Initial position and initial number of turns
	initialX = positions[0] % nbSitesX;
	winding = 0;
}

void State::update(long part, double u) {
	long pos = positions[part]; 

	// Compute next position
	long pos_next, x, y;
	if (u < 0.5) { // Vertical motion
		pos_next = pos + nbSitesX * (2 * (u < 0.25) - 1);
		if (pos_next >= nbSitesTot) {
			pos_next -= nbSitesTot;
		} else if (pos_next < 0) {
			pos_next += nbSitesTot;
		}
	} else { // Horizontal motion
		x = (pos % nbSitesX) + 2 * (u < 0.75) - 1;
		y = pos / nbSitesX;
		if (x == nbSitesX) {
			x = 0;
			if (part==0 && !occupations[y * nbSitesX]) {
				winding++;
			}
		} else if (x == -1) {
			x = nbSitesX - 1;
			if (part==0 && !occupations[y * nbSitesX + x])
				winding--;
		}
		pos_next = y * nbSitesX + x;
	}

	// Exclusion rules
	int o = occupations[pos_next];
	long dpos = (1 - o) * (pos_next - pos);
	positions[part] += dpos; // Move if next site not occupied
	occupations[pos] = o; // Occupied only if doesn't move
	occupations[pos_next] = 1; // Next site is occupied anyways
}

// Print the system into terminal.
void State::visualize(const Parameters &p, const double t, const long hl) {
	std::cout << "\033[2J\033[1;1H";  // Clear screen

	std::map<long,std::string> mymap;
	std::map<long,std::string>::iterator it;
	mymap[positions[0]] = "\e[0;36mX\e[m";
	for (long i = 1 ; i < p.nbParticles ; ++i) {
		if (i == hl)
			mymap[positions[i]] = "\e[0;33mO\e[m";
		else
			mymap[positions[i]] = "O";
	}

	std::cout << "\r";
	for (long i = 0 ; i < p.nbSitesY ; ++i) {
		for (long j = 0 ; j < p.nbSitesX ; ++j) {
			it = mymap.find(i * nbSitesX + j);
			if (it != mymap.end()) {
				std::cout << it->second;
			} else {
				std::cout << '.';
			}
		}
		std::cout << "\n";
	}
	std::cout << "X=" << getXtp() << "\n";
	std::cout << std::fixed << std::setprecision(5) <<  "(t = " << t << ")"
		<< std::endl;
	
	// Sleep
	std::this_thread::sleep_for(std::chrono::milliseconds(p.sleep));
}
