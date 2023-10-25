#include <cassert>
#include "utils.h"
#include "observables.h"

/* Observables */
Observables::Observables(const long nbMoms) : nbMoments(nbMoms) {
	moments.assign(nbMoms, 0);
}

void Observables::fromState(const State &state) {
	long dx = state.getXtp();
	moments[0] = dx;
	for (long i = 1 ; i < nbMoments ; ++i) {
		moments[i] = moments[i-1] * dx;
	}
}

void Observables::add(const Observables &obs) {
	//assert(nbObs == obs.nbObs);
	for (long i = 0 ; i < nbMoments ; ++i) {
		moments[i] += obs.moments[i];
	}
}

void Observables::print(const long N, std::ostream &stream) const {
	for (long i = 0 ; i < nbMoments ; ++i) {
		stream << ((double) moments[i]) / N << " ";
	}
}

/* ObservablesVec */
ObservablesVec::ObservablesVec(const long nbIt, const long nbMoms) :
		nbIters(nbIt), obsVec(nbIt, Observables(nbMoms)) {
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
