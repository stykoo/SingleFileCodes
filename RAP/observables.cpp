#include <cassert>
#include "observables.h"

/* Observables */
Observables::Observables(const long nbMoms, const long nbPts,
		                 const bool computeProf) :
		nbMoments(nbMoms), nbPointsProf(nbPts), computeProf(computeProf) {
	moments.assign(nbMoms, 0);
	if (computeProf)
		profiles.assign(nbMoms, std::vector<double>(nbPointsProf, 0.));
}

void Observables::fromState(const State &state) {
	double dx = state.getXtp();
	moments[0] = dx;

	for (long i = 1 ; i < nbMoments ; ++i)
		moments[i] = moments[i-1] * dx;

	if (computeProf) {
		state.getPosRel(profiles[0], nbPointsProf);
		for (long i = 1 ; i < nbMoments ; ++i) {
			for (long j = 0 ; j < nbPointsProf ; ++j) {
				profiles[i][j] = profiles[i-1][j] * dx;
			}
		}
	}

}

void Observables::add(const Observables &obs) {
	//assert(nbObs == obs.nbObs);
	//assert(nbPointsProf == obs.nbPointsProf);
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

/* ObservablesVec */
ObservablesVec::ObservablesVec(const long nbIt, const long nbMoms,
		                       const long nbPointsProf,
							   const bool computeProf) :
		nbIters(nbIt),
		obsVec(nbIt, Observables(nbMoms, nbPointsProf, computeProf)) {
}

void ObservablesVec::add(const ObservablesVec &ov) {
	//assert(nbIters == ov.nbIters);
	for (long i = 0 ; i < nbIters ; ++i) {
		obsVec[i].add(ov.obsVec[i]);
	}
}

void ObservablesVec::print(const long N, std::ostream &stream) const {
	for (long t = 0 ; t < nbIters ; ++t) {
		stream << t << " ";
		obsVec[t].print(N, stream);
		stream << "\n";
	}
}
