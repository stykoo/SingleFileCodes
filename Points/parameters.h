#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_N_MOMENTS 4
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_OUTPUT_FILE "data"
#define DEFAULT_THREADS 1
#define DEFAULT_VISU_SLEEP 200
#define DEFAULT_LEN_TONKS 0

#include <iostream>
#include <string>
#include <vector>

struct Parameters {
	// Methods
	int check() const;
	void print(std::ostream &stream = std::cout) const;
	int fromCommandLine(int argc, char **argv);

	long nbParticles;  // Number of particles
	double duration;  // Duration of the simulation

	long nbSimuls;  // Number of simulations
	double lenProf;  // Length over which the profiles are computed
	long nbPtsProf;  // Number of points for profiles
	int nbThreads;  // Number of threads
	std::string output;  // Name of the output file
	double lenTonks; // Length of rods in Tonks gas

	bool determ; // Deterministic initial conditions
	bool export_prof;  // Export the profiles
	bool verbose;  // Verbose mode
};

#endif
