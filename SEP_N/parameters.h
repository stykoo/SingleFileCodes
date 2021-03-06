#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_PROBA_RIGHT 0.5
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_OUTPUT_FILE "observables.dat"
#define DEFAULT_CUSTOM_FILE "custom.dat"
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

	long nbSites;  // Number of sites
	long nbParticles;  // Number of vacancies

	double duration;  // Duration of the simulation
	double dt;  // Timestep for export of observables
	long nbSteps;  // This is duration / dt

	long nbSimuls;  // Number of simulations
	int nbThreads;  // Number of threads
	std::string output;  // Name of the output file
	std::string custom;  // Name of the file for custom observables

	long nbTracers;  // Number of tracers
	std::vector<long> initPos;  // Initial positions of the tracers
	std::vector<double> probas;  // Probabilities to jump to the right
	
	bool visu;  // Visualization
	int sleep;  // Number of milliseconds between two visualizations
	bool verbose;  // Verbose mode
};

#endif
