#ifndef STATE_H
#define STATE_H

#include <cmath>
#include <vector>
#include "mkl.h"
#include "mkl_vsl.h"
#include "parameters.h"


template<typename T>
void pbc(T &x, const T L){
	x -= L * std::floor(x / L);
}

template<typename T>
inline void pbcSym(T &x, const T L) {
	x -= L * std::round(x / L);
}


struct State {
	public:
		State(const Parameters &p);
		State(const State &s) = default;
		void init(VSLStreamStatePtr stream);
		void init_determ();
		void reset();
		void update(long part, double u, int dir);
		void updateNoBias(long part, double u);
		double getXtp() const {
			return positions[0] - initialX + (winding * length);
		}
		void getPosRel(std::vector<double> &posr, const long nbPts) const;

	protected:
		const long nbParticles; // Number of particles
		const double length;  // Length
		double initialX; // Initial position of the TP
		long winding; // Number of turns of the TP
		std::vector<double> positions;  // Positions of the particles
};

#endif
