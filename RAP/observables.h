#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "state.h"

class Observables {
	public:
		Observables(const long nbMoms, const long nbPts,
				    const bool computeProf);
		Observables(const Observables &o) = default;
		void fromState(const State &state);
		void add(const Observables &obs);
		void print(const long N, std::ostream &stream=std::cout) const;
		double getProfiles(const long i, const long j) const {
			return profiles[i][j];
		}

	protected:
		const long nbMoments;  // Number of moments
		const long nbPointsProf;  // Number of points in profiles
		const double computeProf; // Compute the profiles
		std::vector<double> moments;  // All the moments
		std::vector<std::vector<double> > profiles;  // Generalized profiles
};

class ObservablesVec {
	public:
		ObservablesVec(const long nbIt, const long nbMoms,
				       const long nbPointsProf, const bool computeProf);
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
