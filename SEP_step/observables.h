#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <iostream>
#include <vector>
#include <string>
#include "state.h"

#define DEFAULT_N_MOMS 6

class Observables {
	public:
		Observables(bool computeProfs, long nbSites);
		Observables(const Observables &o) = default;
		void fromState(const State &state, const State &initialState);
		void add(const Observables &obs);
		void print(const long N, std::ostream &stream=std::cout) const;
		void printProfsP(const long N, std::ostream &stream=std::cout) const;
		void printProfsM(const long N, std::ostream &stream=std::cout) const;

	protected:
		const bool computeProfs; // Compute the profiles or not
		const long nbSites; // Number of sites
		std::vector<long long> moments;  // All the moments
		std::vector<std::vector<long long> > profilesP; // Generalized profiles
		std::vector<std::vector<long long> > profilesM; // Generalized profiles
};

class ObservablesVec {
	public:
		ObservablesVec(const long nbIt, bool computeProfs, long nbSites);
		void add(const ObservablesVec &ov);
		void print(const long N, std::ostream &stream=std::cout) const;
		Observables& operator[](std::size_t idx) {
			return obsVec[idx];
		}
    	const Observables& operator[](std::size_t idx) const {
			return obsVec[idx];
		}

	protected:
		const long nbIters;
		std::vector<Observables> obsVec;
};

#endif
