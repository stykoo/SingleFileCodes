#ifndef STATE_H
#define STATE_H

#include <vector>
#include <mkl.h>
#include <mkl_vsl.h>
#include "parameters.h"

struct State {
	public:
		void init(const Parameters &p, VSLStreamStatePtr stream);
		void init_determ(const Parameters &p);
		void update(const Parameters &p, long part, double u);
		void visualize(const Parameters &p, const double t);
		long getXtp() const { return x_tp; }
		//long getX(const long i) const { return positions[i]; }
		//int getOcc(const long i) const { return occupations[i]; }

	protected:
		long x_tp;
		std::vector<long> positions;  // Positions of the particles
		std::vector<int> occupations;  // Occupations of the sites
};

#endif
