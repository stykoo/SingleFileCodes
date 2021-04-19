#include <cassert>
#include "observables.h"

Observables::Observables(bool computeProfs, long nbSites) :
	computeProfs(computeProfs), nbSites(nbSites) {
	moments.assign(DEFAULT_N_MOMS, 0);
	if (computeProfs) {
		profilesP.assign(DEFAULT_N_OBS, std::vector<long long>(nbSites, 0));
		profilesM.assign(DEFAULT_N_OBS, std::vector<long long>(nbSites, 0));
	}
}

void Observables::fromState(const State &state, const State &initialState) {
	double X = state.getX();
	double xPer = X - initialState.getX();
	moments[0] = xPer;
	for (long i = 1 ; i < DEFAULT_N_MOMS ; ++i) {
		moments[i] = moments[i-1] * xPer;
	}

	if (computeProfs) {
		int ep = state.getOcc(X+1);
		int em = state.getOcc(X-1);
		long xPerEp = xPer * ep;
		long xPerEm = xPer * em;
		long xPer2 = xPer * xPer;
		long xPer2Ep = xPer2 * ep;
		long xPer2Em = xPer2 * em;
		long xPer3 = xPer2 * xPer;
		long xPer3Ep = xPer3 * ep;
		long xPer3Em = xPer3 * em;

		for (long k = 0 ; k < nbSites ; ++k) {
			int er = state.getOcc(X+k);
			profilesP[0][k] = er;
			profilesP[1][k] = ep * er;
			profilesP[2][k] = em * er;
			profilesP[3][k] = xPer * er;
			profilesP[4][k] = xPerEp * er;
			profilesP[5][k] = xPerEm * er;
			profilesP[6][k] = xPer;
			profilesP[7][k] = xPer2 * er;
			profilesP[8][k] = xPer2Ep * er;
			profilesP[9][k] = xPer2Em * er;
			profilesP[10][k] = xPer2;
			profilesP[11][k] = xPer3 * er;
			profilesP[12][k] = xPer3Ep * er;
			profilesP[13][k] = xPer3Em * er;
			profilesP[14][k] = xPer3;
		}
		for (long k = 0 ; k < nbSites ; ++k) {
			int er = state.getOcc(X-k);
			profilesM[0][k] = er;
			profilesM[1][k] = ep * er;
			profilesM[2][k] = em * er;
			profilesM[3][k] = xPer * er;
			profilesM[4][k] = xPerEp * er;
			profilesM[5][k] = xPerEm * er;
			profilesM[6][k] = xPer;
			profilesM[7][k] = xPer2 * er;
			profilesM[8][k] = xPer2Ep * er;
			profilesM[9][k] = xPer2Em * er;
			profilesM[10][k] = xPer2;
			profilesM[11][k] = xPer3 * er;
			profilesM[12][k] = xPer3Ep * er;
			profilesM[13][k] = xPer3Em * er;
			profilesM[14][k] = xPer3;
		}
	}
}

void Observables::add(const Observables &obs) {
	for (long i = 0 ; i < DEFAULT_N_MOMS ; ++i) {
		moments[i] += obs.moments[i];
	}
	if (computeProfs) {
		for (long i = 0 ; i < DEFAULT_N_OBS ; ++i) {
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
		for (long i = 0 ; i < DEFAULT_N_OBS ; ++i) {
			stream << ((double) profilesP[i][j]) / N << " ";
		}
		stream << "\n";
	}
}

void Observables::printProfsM(const long N, std::ostream &stream) const {
	for (long j = 0 ; j < nbSites ; ++j) {
		stream << -j << " ";
		for (long i = 0 ; i < DEFAULT_N_OBS ; ++i) {
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
