#include <algorithm>
#include "state.h"

State::State(const Parameters &p) :
		nbParticles(p.nbParticles), length((double) p.nbParticles),
		initialX(0.0), winding(0) {
}

void State::init(VSLStreamStatePtr stream) {
	positions.resize(nbParticles);

	// One sets the position of the 1st particle and assign the other at random
	// Be careful that assigning all positions at random gives a non uniform
	// distribution on the circle.
	positions[0] = 0.;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, nbParticles-1,
			     positions.data()+1, 0, length);
	std::sort(positions.begin()+1, positions.end());

	reset();
}

void State::init_stat(VSLStreamStatePtr stream) {
	positions.resize(nbParticles);
	// aux is actually not necessary (everything can be done in positions)
	aux.resize(nbParticles); 
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nbParticles,
			      aux.data(), 0., 1.);

	for (long i = 0 ; i < nbParticles ; ++i)
		aux[i] *= aux[i];

	double s = 0.;
	for (long i = 0 ; i < nbParticles ; ++i)
		s += aux[i];

	positions[0] = 0.;
	for (long i = 1 ; i < nbParticles ; ++i) {
		positions[i] = positions[i-1] + length * aux[i-1] / s;
	}

	reset();
}

void State::init_determ() {
	positions.resize(nbParticles);

	for (long i = 0 ; i < nbParticles ; ++i) {
		positions[i] = (i * length) / nbParticles;
	}
	reset();
}

void State::reset() {
	// Initial position and initial number of turns
	initialX = positions[0];
	winding = 0;
}

void State::update(long part, double u, int dir) {
	// This version is a little bit faster than the previous one
	// but is harder to read.
	int eps = 2 * dir - 1; // 1 or -1

	// Neighbor particle in given direction
	long part_nbr = part + eps;
	if (part_nbr == -1)
		part_nbr = nbParticles - 1;
	else if (part_nbr == nbParticles)
		part_nbr = 0;

	// Difference of positions (should be in the correct direction)
	double dpos = positions[part_nbr] - positions[part];
	dpos += (dpos * eps < 0) * eps * length;

	// Motion (in the direction of eps)
	positions[part] += u * dpos;

	// Periodic boundary conditions
	if (eps * positions[part] > dir * length) {
		positions[part] -= (eps * length);
		winding += eps * (part == 0); // The tracer performed one turn
	} 

	/*long part_nbr;
	double dpos;

	// Move right
	if (dir) {
		part_nbr = part + 1;
		if (part_nbr == nbParticles)
			part_nbr = 0;
		dpos = positions[part_nbr] - positions[part];
		if (dpos < 0) // The displacement should be positive
			dpos += length;
		positions[part] += u * dpos;
		// Periodic boundary conditions
		if (positions[part] >= length) {
			positions[part] -= length;
			if (part == 0) // The tracer performed one turn
				winding++;
		}
	// Move left
	} else {
		part_nbr = part - 1;
		if (part_nbr == -1)
			part_nbr = nbParticles - 1;
		dpos = positions[part_nbr] - positions[part];
		if (dpos > 0) // The displacement should be negative
			dpos -= length;
		positions[part] += u * dpos;
		// Periodic boundary conditions
		if (positions[part] < 0) {
			positions[part] += length;
			if (part == 0) // The tracer performed one turn
				winding--;
		}
	} */
}

void State::getPosRel(std::vector<double> &posr, const long nbPts) const {
	double x;
	for (long j = 0 ; j < nbPts ; ++j) {
		posr[j] = 0.0;
	}
	for (long i = 1 ; i < nbParticles ; ++i) {
		x = (positions[i] - positions[0]) / length;
		//x -= std::floor(x);
		//if (x < 0)  x += 1.0;
		x += (x < 0);
		posr[(size_t) (x * nbPts)] += 1.0;
	}
}
