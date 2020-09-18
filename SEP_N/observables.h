#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <iostream>
#include <vector>
#include <string>
#include "state.h"

typedef std::vector< std::vector<int> > CoeffsCustom;

class ObservablesCustom {
	public:
		ObservablesCustom(const CoeffsCustom &coeffs);
		ObservablesCustom(const ObservablesCustom &o) = default;
		void fromState(const State &state, const State &initialState,
				       const long nbSites);
		void add(const ObservablesCustom &obs);
		void print(const long N, std::ostream &stream=std::cout) const;

	protected:
		const long nbObs;  // Number of tracers
		const long nbTracers;  // Number of tracers
		std::vector< std::vector<int> > coeffs;  // Coefficients of the moments
		std::vector<long long> moments;  // All the moments
};

class ObservablesVecCustom {
	public:
		ObservablesVecCustom(const long nbIt, const CoeffsCustom &coeffs);
		void add(const ObservablesVecCustom &ov);
		void print(const long N, std::ostream &stream=std::cout) const;
		ObservablesCustom& operator[](std::size_t idx) {
			return obsVec[idx];
		}
    	const ObservablesCustom& operator[](std::size_t idx) const {
			return obsVec[idx];
		}

	protected:
		const long nbIters;
		std::vector<ObservablesCustom> obsVec;
};

#endif
