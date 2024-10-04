#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define ALIGNMENT 32

#include <boost/align/aligned_allocator.hpp>
#include "parameters.h"

template <typename T>
using aligned_vector = std::vector<T,
	  boost::alignment::aligned_allocator<T, ALIGNMENT>>;

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
		aligned_vector<long> mom1, mom2, mom3, mom4;  // Moments
};

void add_arrays(long *a1, const long *a2, long n);

#endif
