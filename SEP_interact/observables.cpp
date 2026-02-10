#include <cassert>
#include "observables.h"

Observables::Observables() {
	moments.assign(DEFAULT_N_MOMS, 0);
}

void Observables::fromState(const State &state) {
	double m = state.computeDispl();
	double m2 = state.computeDisplSq();
	double j = (double) state.getCurrent();
	moments[0] = m;
	moments[1] = m2;
	moments[2] = m2 - m * m;
	moments[3] = j;
	moments[4] = j * j;
}

void Observables::add(const Observables &obs) {
	for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
		moments[i] += obs.moments[i];
	}
}

void Observables::print(const long N, std::ostream &stream) const {
	for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
		stream << ((double) moments[i]) / N << " ";
	}
}

ObservablesVec::ObservablesVec(const long nbIt) :
		nbIters(nbIt), obsVec(nbIt, Observables()) {
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
