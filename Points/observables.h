#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <iostream>
#include <vector>
#include <array>
#include <string>

class Observables {
	public:
		Observables(const long nbMoms, const double len, const long nbPts,
				    const bool computeProf, const bool computeCorrel);
		Observables(const Observables &o) = default;
		void compute(const double dx, const std::vector<double> &dpos,
	                 const long n_parts, const long mid);
		void add(const Observables &obs);
		void print(const long N, std::ostream &stream=std::cout) const;
		double getProfiles(const long i, const long j) const {
			return profiles[i][j];
		}
		double getCorrelsP(const long i, const long j) const {
			return correlsP[i][j];
		}
		double getCorrelsM(const long i, const long j) const {
			return correlsM[i][j];
		}

	protected:
		const long nbMoments;  // Number of moments
		const double lenProf; // Length for profiles
		const long nbPointsProf;  // Number of points in profiles
		const bool computeProf; // Compute the profiles
		const bool computeCorrel; // Compute the correlations
		std::vector<double> moments;  // All the moments
		std::vector<std::vector<double> > profiles;  // Generalized profiles
		std::vector<std::vector<double> > correlsP;  // Generalized correlations
		std::vector<std::vector<double> > correlsM;  // Generalized correlations
};

#endif
