#include <cassert>
#include "observables.h"

/* Observables */
Observables::Observables(const long nbMoms, const double len, const long nbPts,
		                 const bool computeProf) :
		nbMoments(nbMoms), lenProf(len), nbPointsProf(nbPts),
		computeProf(computeProf) {
	moments.assign(nbMoms, 0);
	if (computeProf)
		profiles.assign(nbMoms, std::vector<double>(nbPointsProf, 0.));
}

void Observables::compute(const double dx, const std::vector<double> &dpos,
		                  const long n_parts, const long mid) {
	moments[0] = dx;
	for (long i = 1 ; i < nbMoments ; ++i)
		moments[i] = moments[i-1] * dx;

	if (computeProf) {
		for (long i = 0 ; i < nbMoments ; ++i) {
			for (long j = 0 ; j < nbPointsProf ; ++j) {
				profiles[i][j] = 0.0;
			}
		}

		for (long i = 0 ; i < mid ; ++i) {
			//assert(dpos[i] <= 0.);
			if (-dpos[i] < lenProf) {
				int k = (int) (-dpos[i] * nbPointsProf / lenProf);
				double u = 1;
				for (long j = 0 ; j < nbMoments ; ++j) {
					profiles[j][k] += u;
					u *= -dx;
				}
			}
		}
		for (long i = mid + 1 ; i < n_parts ; ++i) {
			//assert(dpos[i] >= 0.);
			if (dpos[i] < lenProf) {
				int k = (int) (dpos[i] * nbPointsProf / lenProf);
				double u = 1;
				for (long j = 0 ; j < nbMoments ; ++j) {
					profiles[j][k] += u;
					u *= dx;
				}
			}
		}
	}
}

void Observables::add(const Observables &obs) {
	for (long i = 0 ; i < nbMoments ; ++i) {
		moments[i] += obs.moments[i];
		if (computeProf) {
			for (long j = 0 ; j < nbPointsProf ; ++j) {
				profiles[i][j] += obs.profiles[i][j];
			}
		}
	}
}

void Observables::print(const long N, std::ostream &stream) const {
	for (long i = 0 ; i < nbMoments ; ++i) {
		stream << moments[i] / N << " ";
	}
}
