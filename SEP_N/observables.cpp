#include <cassert>
#include "utils.h"
#include "observables.h"

ObservablesCustom::ObservablesCustom(const CoeffsCustom &c) :
		nbObs(c.size()), nbTracers(c[0].size()), coeffs(c) {
	moments.assign(nbObs, 0);
}

void ObservablesCustom::fromState(const State &state,
		                          const State &initialState,
		                          const long nbSites) {
	//std::vector<long> xsPer(nbTracers);
	//state.displ(xsPer, initialState, nbTracers, nbSites);
	/*
	for (long i = 0 ; i < nbTracers ; ++i) {
		xsPer[i] = periodicBCsym(state.getX(i) - initialState.getX(i),
				                 nbSites);
	}
	*/
	for (long i = 0 ; i < nbObs ; ++i) {
		moments[i] = 1;
	}
	long dx;
	for (long j = 0 ; j < nbTracers ; ++j) {
		dx = periodicBCsym(state.getX(j) - initialState.getX(j), nbSites);
		for (long i = 0 ; i < nbObs ; ++i) {
			moments[i] *= mypow(dx, coeffs[i][j]);
		}
	}
}

void ObservablesCustom::add(const ObservablesCustom &obs) {
	//assert(nbObs == obs.nbObs);
	for (long i = 0 ; i < nbObs ; ++i) {
		moments[i] += obs.moments[i];
	}
}

void ObservablesCustom::print(const long N, std::ostream &stream) const {
	for (long i = 0 ; i < nbObs ; ++i) {
		stream << ((double) moments[i]) / N << " ";
	}
}

ObservablesVecCustom::ObservablesVecCustom(const long nbIt,
		                                   const CoeffsCustom &c) :
		nbIters(nbIt), obsVec(nbIt, ObservablesCustom(c)) {
}

void ObservablesVecCustom::add(const ObservablesVecCustom &ov) {
	//assert(nbIters == ov.nbIters);
	for (long i = 0 ; i < nbIters ; ++i) {
		obsVec[i].add(ov.obsVec[i]);
	}
}

void ObservablesVecCustom::print(const long N, std::ostream &stream) const {
	for (long t = 0 ; t < nbIters ; ++t) {
		stream << t << " ";
		obsVec[t].print(N, stream);
		stream << "\n";
	}
}
