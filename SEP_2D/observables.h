#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "state.h"


class Observables {
	public:
		Observables(const long nbMoms);
		Observables(const Observables &o) = default;
		void fromState(const State &state);
		void add(const Observables &obs);
		void print(const long N, std::ostream &stream=std::cout) const;

	protected:
		const long nbMoments;  // Number of moments
		std::vector<long long> moments;  // All the moments
};

class ObservablesVec {
	public:
		ObservablesVec(const long nbIt, const long nbMoms);
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
