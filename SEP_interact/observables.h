#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <iostream>
#include <vector>
#include <string>
#include "state.h"

#define DEFAULT_N_MOMS 5

class Observables {
	public:
		Observables();
		Observables(const Observables &o) = default;
		void fromState(const State &state);
		void add(const Observables &obs);
		void print(const long N, std::ostream &stream=std::cout) const;

	protected:
		std::vector<double> moments;  // The moments of the displacement
};

class ObservablesVec {
	public:
		ObservablesVec(const long nbIt);
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
