#include <cassert>
#include "observables.h"

Observables::Observables(bool computeProfs, long nbSites) :
	computeProfs(computeProfs), nbSites(nbSites) {
	moments.assign(DEFAULT_N_MOMS, 0);
	if (computeProfs) {
		profilesP.assign(DEFAULT_N_MOMS, std::vector<long long>(nbSites, 0));
		profilesM.assign(DEFAULT_N_MOMS, std::vector<long long>(nbSites, 0));
	}
}

void Observables::fromState(const State &state, const State &initialState) {
	double X = state.getX();
	double dx = X - initialState.getX();
	moments[0] = dx;
	for (long i = 1 ; i < DEFAULT_N_MOMS ; ++i) {
		moments[i] = moments[i-1] * dx;
	}

	// TODO profiles
	if (computeProfs) {
		//long ep = state.occupations[X+1];
		//long em = state.occupations[X-1];
		for (long k = 0 ; k < nbSites ; ++k) {
			profilesP[0][k] = state.getOcc(X+k);
			profilesM[0][k] = state.getOcc(X-k);
			for (long i = 1 ; i < DEFAULT_N_MOMS ; ++i) {
				profilesP[i][k] = profilesP[i-1][k] * dx;
				profilesM[i][k] = profilesM[i-1][k] * dx;
			}
		}
	}
}

void Observables::add(const Observables &obs) {
	for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
		moments[i] += obs.moments[i];
	}
	if (computeProfs) {
		for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
			for (long j = 0 ; j < nbSites ; ++j) {
				profilesP[i][j] += obs.profilesP[i][j];
				profilesM[i][j] += obs.profilesM[i][j];
			}
		}
	}
}

void Observables::print(const long N, std::ostream &stream) const {
	for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
		stream << ((double) moments[i]) / N << " ";
	}
}

void Observables::printProfsP(const long N, std::ostream &stream) const {
	for (long j = 0 ; j < nbSites ; ++j) {
		stream << j << " ";
		for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
			stream << ((double) profilesP[i][j]) / N << " ";
		}
		stream << "\n";
	}
}

void Observables::printProfsM(const long N, std::ostream &stream) const {
	for (long j = 0 ; j < nbSites ; ++j) {
		stream << -j << " ";
		for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
			stream << ((double) profilesM[i][j]) / N << " ";
		}
		stream << "\n";
	}
}

ObservablesVec::ObservablesVec(const long nbIt, bool computeProfs,
		                       long nbSites) :
		nbIters(nbIt), obsVec(nbIt, Observables(computeProfs, nbSites)) {
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
