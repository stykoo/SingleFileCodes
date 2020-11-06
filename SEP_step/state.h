#ifndef STATE_H
#define STATE_H

#include <deque>
#include <mkl.h>
#include <mkl_vsl.h>
#include "parameters.h"

struct State {
	public:
		void init(const Parameters &p, VSLStreamStatePtr stream);
		void init_determ(const Parameters &p);
		double update(const Parameters &p, double u, VSLStreamStatePtr stream);
		void visualize(const Parameters &p, const double t);
		long getX() const { return pos_tp; }
		int getOcc(const long i) const { return occupations[i]; }

	protected:
		long nbParts; // Number of particles (except the TP)
		long pos_tp; // Position of the TP
		std::deque<long> positions;  // Positions of the particles
		std::vector<int> occupations;  // Occupations of the sites
};

#endif
