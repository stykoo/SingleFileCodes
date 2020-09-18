#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <thread>
#include <iomanip>
#include <random>
#include "simul.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	std::vector<std::thread> threads;

	std::vector< Observables >
		allSumsObs(p.nbThreads,
				   Observables(DEFAULT_N_MOMENTS, p.lenProf, p.nbPtsProf,
					           p.export_prof));

	long nbSimulsPerThread = p.nbSimuls / p.nbThreads;
	
	if (p.nbSimuls % p.nbThreads != 0) {
		std::cerr << "Warning: nbSimuls is not a multiple of nbThreads."
			<< nbSimulsPerThread * p.nbThreads << " simulations will be done."
			<< std::endl;
	}

	// Threads
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
	        nbSimulsPerThread, std::ref(allSumsObs[i]), rd())); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	// Initialize the total sum
	Observables sumObs(DEFAULT_N_MOMENTS, p.lenProf, p.nbPtsProf, p.export_prof);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		sumObs.add(allSumsObs[k]);
	}

	int status;
	status = exportObservables(sumObs, p);
	if (!status && p.export_prof)
		status = exportProfiles(sumObs, p);	
	return status;
}

// Run several simulation and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                    Observables &sumObs, const unsigned int seed) {
    // Random generator
	VSLStreamStatePtr stream;
	vslNewStream(&stream, CUSTOM_RNG, seed);

	Observables obs(DEFAULT_N_MOMENTS, p.lenProf, p.nbPtsProf, p.export_prof);

	for (long s = 0 ; s < nbSimuls ; ++s) {
		if (p.verbose) {
			std::cout << "Simulation " << s + 1 
				<< " on thread " << std::this_thread::get_id() << std::endl;
		}

		// Run a single simulation
		runOneSimulation(p, obs, stream);
		// Add the observables to the sum
		sumObs.add(obs);
	}

	vslDeleteStream(&stream);
}

// Run a single simulation and compute the moments.
void runOneSimulation(const Parameters &p, Observables &obs,
					  VSLStreamStatePtr stream) {
	//double length = (double) p.nbParticles;
	double length = p.nbParticles * (1 - p.lenTonks);
	double sigma = sqrt(p.duration); 
	long mid = p.nbParticles / 2;

	std::vector<double> pos_0(p.nbParticles), pos_t(p.nbParticles);

	// Initial configuration
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, p.nbParticles,
			     pos_0.data(), 0, length);
	// Final configuration
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, p.nbParticles,
				  pos_t.data(), 0.0, sigma);
	for (long i = 0 ; i < p.nbParticles ; ++i)
		pos_t[i] += pos_0[i];

	if (p.lenTonks) { // Tonks gas
		std::sort(pos_0.begin(), pos_0.end());
		std::sort(pos_t.begin(), pos_t.end());
		for (long i = 1 ; i < p.nbParticles ; ++i) {
			pos_0[i] += i * p.lenTonks;
			pos_t[i] += i * p.lenTonks;
		}
	} else { // Pointlike particles
		// Select median (no need to do the full sort)
		std::nth_element(pos_0.begin(), pos_0.begin() + mid, pos_0.end());
		std::nth_element(pos_t.begin(), pos_t.begin() + mid, pos_t.end());
	}

	double xfin = pos_t[mid];
	double dx = xfin - pos_0[mid];
	
	if (p.export_prof) {
		for (long i = 0 ; i < p.nbParticles ; ++i) {
			pos_t[i] -= xfin;
		}
	}
	
	// Compute observables
	obs.compute(dx, pos_t, p.nbParticles, mid);
}

// Export the observables to a file.
int exportObservables(const Observables &sumObs, const Parameters &p) {
	std::ofstream file(p.output + ".dat");
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# SFPoint (" << __DATE__ <<  ", " << __TIME__ << "): ";
	p.print(file);
	file << "\n# t";
	for (size_t i = 0 ; i < DEFAULT_N_MOMENTS ; ++i) {
		file << " x^" << i+1;
	}
	file << "\n";

	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

	// Data (we write the average and not the sum)
	file << p.duration << " ";
	sumObs.print(p.nbSimuls, file);
    file << "\n";
	file.close();
	return 0;
}

int exportProfiles(const Observables &sumObs, const Parameters &p) {
	std::string fname = p.output + "_prof.dat";
	double fac = p.lenProf / p.nbPtsProf;

	std::ofstream file(fname);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# RAP_Profiles (" << __DATE__ <<  ", " << __TIME__ << "): ";
	p.print(file);
	file << "\n# x";
	for (size_t i = 0 ; i < DEFAULT_N_MOMENTS ; ++i) {
		file << " rho(x)*X^" << i;
	}
	file << "\n";
	for (long j = 0 ; j < p.nbPtsProf ; ++j) {
		file << j * fac;
		for (size_t k = 0 ; k < DEFAULT_N_MOMENTS ; ++k) {
			// Factor 2 because of symmetrization
			file << " " << sumObs.getProfiles(k, j) / fac / p.nbSimuls / 2;
		}
		file << "\n";
	}

	file.close();

	return 0;
}
