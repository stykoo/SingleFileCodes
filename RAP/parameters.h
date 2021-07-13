#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_PROBA_RIGHT 0.5
#define DEFAULT_N_MOMENTS 4
#define DEFAULT_N_PTS_PROF 10
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_OUTPUT_FILE "data"
#define DEFAULT_THREADS 1
#define DEFAULT_VISU_SLEEP 200

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
	double duration_therm;  // Duration of thermalization
	double dt;  // Timestep for export of observables
	long nbSteps;  // This is duration / dt

	long nbSimuls;  // Number of simulations
	long nbPtsProfRel;  // Number of points for profiles per particle
	long nbPtsProf;  // Number of points for profiles
	int nbThreads;  // Number of threads
	std::string output;  // Name of the output file

	double probTP; // Jump probability of the TP to the right
	bool determ; // Deterministic initial conditions
	bool stat; // Stationary-state initial conditions
	bool export_prof;  // Export the profiles
	bool verbose;  // Verbose mode
};

#endif
