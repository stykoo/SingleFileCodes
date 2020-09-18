#ifndef STATE_H
#define STATE_H

#include <vector>
// #include <random>
#include <mkl.h>
#include <mkl_vsl.h>
#include "parameters.h"

struct State {
	public:
		// void init(const Parameters &p, std::mt19937 &rndGen);
		void init(const Parameters &p, VSLStreamStatePtr stream);
		void update(const Parameters &p, long part, double u);
		// void displ(std::vector<long> &dxs, const State &s_ini, const long n,
		//            const long len) const;
		void visualize(const Parameters &p, const long t);
		long getX(const long i) const { return positions[i]; }
		int getOcc(const long i) const { return occupations[i]; }

	protected:
		std::vector<long> positions;  // Positions of the particles
		std::vector<int> occupations;  // Occupations of the sites
};

#endif
