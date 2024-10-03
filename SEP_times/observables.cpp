#include <string>
#include <fstream>
#include <cassert>
#include "observables.h"

Observables::Observables(const long n) :
		n1(n), n2(n*(n+1)/2), n3(n*(n+1)*(n+2)/2) {
	mom1.assign(n1, 0);
	mom2.assign(n2, 0);
	mom3.assign(n3, 0);
	mom4.assign(n3, 0);
}

void Observables::recordPos(const long i, const long x) {
	mom1[i] = x;
}

void Observables::computeMoments() {
	// Two-point moment
	long q = 0;
	for (long i = 0 ; i < n1 ; ++i) {
		for (long j = 0 ; j <= i ; ++j) {
			mom2[q] = mom1[i] * mom1[j];
			q++;
		}
	}
	// Three-point and four-point moments (with final time)
	q = 0;
	long long xf = mom1[n1-1];
	for (long i = 0 ; i < n1 ; ++i) {
		for (long j = 0 ; j <= i ; ++j) {
			for (long k = 0 ; k <= j ; ++k) {
				mom3[q] = mom1[i] * mom1[j] * mom1[k];
				mom4[q] = mom3[q] * xf;
				q++;
			}
		}
	}
}

void Observables::add(const Observables &obs) {
	for (long i = 0 ; i < n1 ; ++i) {
	mom1[i] += obs.mom1[i];
	}
	for (long i = 0 ; i < n2 ; ++i) {
		mom2[i] += obs.mom2[i];
	}
	for (long i = 0 ; i < n3 ; ++i) {
		mom3[i] += obs.mom3[i];
		mom4[i] += obs.mom4[i];
	}
}

int Observables::exportMoments(const Parameters &p) const {
	std::ofstream file;
	
	// First moment
	file.open(p.output + "_mom1.dat");
	if (!file.is_open()) {
		return 1;
	}

	file << "# SEP_times (" << __DATE__ <<  ", " << __TIME__ << "): ";
	p.print(file);
	file << "\n# t x\n";
	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);
	for (long i = 0 ; i < n1 ; ++i) {
		file << i * p.dt << " " << mom1[i] / ((double) p.nbSimuls) << "\n";
	}
	file.close();

	// Second moment
	file.open(p.output + "_mom2.dat");
	if (!file.is_open()) {
		return 1;
	}
	file << "# SEP_times (" << __DATE__ <<  ", " << __TIME__ << "): ";
	p.print(file);
	file << "\n# t1 t2 x1*x2\n";
	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);
	long q = 0;
	for (long i = 0 ; i < n1 ; ++i) {
		for (long j = 0 ; j <= i ; ++j) {
			file << j * p.dt << " " << i * p.dt << " "
			   << mom2[q] / ((double) p.nbSimuls) << "\n";
			q++;
		}
	}
	file.close();

	// Second moment
	file.open(p.output + "_moms34.dat");
	if (!file.is_open()) {
		return 1;
	}
	file << "# SEP_times (" << __DATE__ <<  ", " << __TIME__ << "): ";
	p.print(file);
	file << "\n# t1 t2 t3 x1*x2*x3 x1*x2*x3*xf\n";
	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);
	q = 0;
	for (long i = 0 ; i < n1 ; ++i) {
		for (long j = 0 ; j <= i ; ++j) {
			for (long k = 0 ; k <= j ; ++k) {
				file << k * p.dt << " " << j * p.dt << " " << i * p.dt << " "
				   << mom3[q] / ((double) p.nbSimuls) << " "
				   << mom4[q] / ((double) p.nbSimuls) << "\n";
				q++;
			}
		}
	}
	file.close();

	return 0;
}
