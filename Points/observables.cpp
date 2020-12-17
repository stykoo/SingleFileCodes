#include <cassert>
#include "observables.h"

/* Observables */
Observables::Observables(const long nbMoms, const double len, const long nbPts,
		                 const bool computeProf, const bool computeCorrel) :
		nbMoments(nbMoms), lenProf(len), nbPointsProf(nbPts),
		computeProf(computeProf), computeCorrel(computeCorrel) {
	moments.assign(nbMoms, 0);
	if (computeProf)
		profiles.assign(nbMoms, std::vector<double>(nbPointsProf, 0.));
	if (computeCorrel) {
		correlsP.assign(nbMoms, std::vector<double>(nbPointsProf, 0.));
		correlsM.assign(nbMoms, std::vector<double>(nbPointsProf, 0.));
	}
}

void Observables::compute(const double dx, const std::vector<double> &dpos,
		                  const long n_parts, const long mid) {
	int nM = 0, nP = 0; // Number of particles in bins closest to the TP

	// Moments
	moments[0] = dx;
	for (long i = 1 ; i < nbMoments ; ++i)
		moments[i] = moments[i-1] * dx;

	// Profiles
	if (computeProf) {
		for (long i = 0 ; i < nbMoments ; ++i) {
			for (long j = 0 ; j < nbPointsProf ; ++j) {
				profiles[i][j] = 0.0;
				if (computeCorrel) {
					correlsP[i][j] = 0.0;
					correlsM[i][j] = 0.0;
				}
			}
		}

		if (computeCorrel) {
			for (long i = 0 ; i < mid ; ++i)
				if (-dpos[i] * nbPointsProf < lenProf)
					nM++;
			for (long i = mid + 1 ; i < n_parts ; ++i)
				if (dpos[i] * nbPointsProf < lenProf)
					nP++;
		}

		for (long i = 0 ; i < mid ; ++i) {
			if (-dpos[i] < lenProf) {
				int k = (int) (-dpos[i] * nbPointsProf / lenProf);
				double u = 1;
				for (long j = 0 ; j < nbMoments ; ++j) {
					profiles[j][k] += u;
					if (computeCorrel) {
						correlsP[j][k] += nM * u;
						correlsM[j][k] += nP * u;
					}
					u *= -dx;
				}
			}
		}
		for (long i = mid + 1 ; i < n_parts ; ++i) {
			if (dpos[i] < lenProf) {
				int k = (int) (dpos[i] * nbPointsProf / lenProf);
				double u = 1;
				for (long j = 0 ; j < nbMoments ; ++j) {
					profiles[j][k] += u;
					if (computeCorrel) {
						correlsP[j][k] += nP * u;
						correlsM[j][k] += nM * u;
					}
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
		if (computeCorrel) {
			for (long j = 0 ; j < nbPointsProf ; ++j) {
				correlsP[i][j] += obs.correlsP[i][j];
				correlsM[i][j] += obs.correlsM[i][j];
			}
		}
	}
}

void Observables::print(const long N, std::ostream &stream) const {
	for (long i = 0 ; i < nbMoments ; ++i) {
		stream << moments[i] / N << " ";
	}
}
