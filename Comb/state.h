#ifndef STATE_H
#define STATE_H

#include <vector>
#include "mkl.h"
#include "mkl_vsl.h"
#include "parameters.h"

struct State {
	public:
		State(const Parameters &p);
		State(const State &s) = default;
		void init(VSLStreamStatePtr stream);
		void init_determ();
		void reset();
		void update(long part, double u);
		void visualize(const Parameters &p, const double t, const long hl=-1);
		//long getX(const long i) const { return positions[i] % nbSitesX; }
		//long getY(const long i) const { return positions[i] / nbSitesX; }
		// int getOcc(const long i) const { return occupations[i]; }
		long getXtp() const {
			return positions[0] - initialX + (winding * nbSitesX);
		}		

	protected:
		const long nbSitesX;
		const long nbSitesY;
		const long nbSitesTot;
		const long nbParticles;
		double probTP; // Jump probability of the TP to the right
		double initialX; // Initial position of the TP
		long winding; // Number of turns of the TP
		std::vector<long> positions;  // Positions of the particles
		std::vector<int> occupations;  // Occupations of the sites
};

#endif
