#include <cmath>
#include <map>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <thread>
#include "state.h"

void State::init(const Parameters &p, VSLStreamStatePtr stream) {
	positions.clear();
	occupations.assign(p.nbSites, 0);

	long mid = p.nbSites / 2;
	pos_tp = mid;
	occupations[mid] = 1;
	viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, mid,
				   occupations.data(), p.rhoM);
	viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, p.nbSites-mid-1,
				   occupations.data()+mid+1, p.rhoP);

	nbParts = 0;
	for (long i = 0 ; i < mid ; ++i) {
		if (occupations[i]) {
			positions.push_back(i);
			nbParts++;
		}
	}
	for (long i = mid+1 ; i < p.nbSites ; ++i) {
		if (occupations[i]) {
			positions.push_back(i);
			nbParts++;
		}
	}

}

void State::init_determ(const Parameters &p) {
	positions.clear();
	occupations.assign(p.nbSites, 0);

	long mid = p.nbSites / 2;
	pos_tp = mid;
	occupations[mid] = 1;

	long l;

	nbParts = 0;
	if (p.rhoM <= 0.5) {
		l = (long) (1. / p.rhoM);
		for (long i = mid-l ; i>=0 ; i-=l) {
			occupations[i] = 1;
			positions.push_front(i);
			nbParts++;
		}
	} else {
		l = (long) std::round(1. / (1. - p.rhoM));
		for (long i = mid-1 ; i>=0 ; --i) {
			if ((mid-i) % l != 0) {
				occupations[i] = 1;
				positions.push_front(i);
				nbParts++;
			}
		}
	}

	if (p.rhoP <= 0.5) {
		l = (long) (1. / p.rhoP);
		for (long i = mid+l ; i<p.nbSites ; i+=l) {
			occupations[i] = 1;
			positions.push_back(i);
			nbParts++;
		}
	} else {
		l = (long) std::round(1. / (1. - p.rhoP));
		for (long i = mid+1 ; i<p.nbSites ; ++i) {
			if ((i-mid) % l != 0) {
				occupations[i] = 1;
				positions.push_back(i);
				nbParts++;
			}
		}
	}
}

double State::update(const Parameters &p, double u, VSLStreamStatePtr stream) {
	int part;
	double dt;
	// bath particles + TP + left / right reservoirs -> nbParts + 3
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &part,
				 0, nbParts+3);
	vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream, 1,
			         &dt, 0.0, 1./(nbParts+3));

	if (part < nbParts) {
		// Bath particles
		long pos = positions[part]; 
		long pos_next = pos + 2 * (u < DEFAULT_PROBA_RIGHT) - 1;

		if (pos_next == p.nbSites) { // Jump to right reservoir
			if (u < 0.5 * (1. - p.rhoP)) {
				occupations[pos] = 0;
				positions.pop_back();
				nbParts--;
			}
		} else if (pos_next == -1) { // Jump to left reservoir
			if (u-0.5 < 0.5 * (1. - p.rhoM)) {
				occupations[pos] = 0;
				positions.pop_front();
				nbParts--;
			}
		} else {
			int o = occupations[pos_next];
			long dpos = (1 - o) * (pos_next - pos);
			positions[part] += dpos; // Move if next site not occupied
			occupations[pos] = o; // Occupied only if doesn't move
			occupations[pos_next] = 1; // Next site is occupied anyways
		}
	} else if (part == nbParts) {
		// TP (assuming it doesn't jump to the reservoirs...)
		long pos_next = pos_tp + 2 * (u < DEFAULT_PROBA_RIGHT) - 1;
		if (pos_next == p.nbSites || pos_next == -1) {
			std::cerr << "The TP jumps to a reservoir!" << std::endl;
		}
		int o = occupations[pos_next];
		long dpos = (1 - o) * (pos_next - pos_tp);
		pos_tp += dpos; // Move if next site not occupied
		occupations[pos_tp] = o; // Occupied only if doesn't move
		occupations[pos_next] = 1; // Next site is occupied anyways
	} else if (part == nbParts+1 && u < 0.5 * p.rhoP &&
			   !occupations[p.nbSites-1]) {
		// Jump from right reservoir
		occupations[p.nbSites-1] = 1;
		positions.push_back(p.nbSites-1);
		nbParts++;
	} else if (u < 0.5 * p.rhoM && !occupations[0]) {
		// Jump from left reservoir
		occupations[0] = 1;
		positions.push_front(0);
		nbParts++;
	}

	return dt;
}

// Print the system into terminal.
void State::visualize(const Parameters &p, const double t) {
	std::map<long,char> mymap;
	std::map<long,char>::iterator it;
	mymap[pos_tp] = 'X';
	for (long i = 0 ; i < nbParts ; ++i) {
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
