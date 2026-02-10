#ifndef STATE_H
#define STATE_H

#include <mkl.h>
#include <mkl_vsl.h>
#include "parameters.h"

struct State {
	public:
		State(const Parameters &_p, VSLStreamStatePtr stream);
		void init(VSLStreamStatePtr stream);
		void init_determ();
		void update(long part, double u);
		double computeDispl() const;
		double computeDisplSq() const;
		void visualize(const double t);
		int getOcc(const long i) const { return occupations[i]; }
		long getDispl(const long i) const { return displacements[i]; }
		long getCurrent() const { return current; }

	protected:
		const Parameters p;
		std::vector<long> positions;  // Positions of the particles
		std::vector<long> displacements;  // Displacements of the particles
		std::vector<int> occupations;  // Occupations of the sites
        long current;
};

#endif
