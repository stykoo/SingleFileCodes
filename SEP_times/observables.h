#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "parameters.h"

class Observables {
	public:
		Observables(const long n);
		Observables(const Observables &o) = default;
		void add(const Observables &obs);
		void recordPos(const long i, const long x);
		void computeMoments();
		int exportMoments(const Parameters &p) const;

	protected:
		const long n1, n2, n3;  // Number of points
		std::vector<long long> mom1, mom2, mom3, mom4;  // Moments
};

#endif
