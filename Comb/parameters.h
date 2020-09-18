#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_PROBA_RIGHT 0.5
#define DEFAULT_N_MOMENTS 4
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_OUTPUT_FILE "observables.dat"
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

	long nbSitesX;  // Number of sites in the x direction
	long nbSitesY;  // Number of sites in the y direction
	long nbParticles;  // Number of particles

	double duration;  // Duration of the simulation
	double dt;  // Timestep for export of observables
	long nbSteps;  // This is duration / dt

	long nbSimuls;  // Number of simulations
	int nbThreads;  // Number of threads
	std::string output;  // Name of the output file

	double probTP; // Jump probability of the TP to the right

	bool determ; // Deterministic initial conditions
	
	bool visu;  // Visualization
	int sleep;  // Number of milliseconds between two visualizations
	bool verbose;  // Verbose mode
};

#endif
