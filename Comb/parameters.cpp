#include <exception>
#include <boost/program_options.hpp>
#include "parameters.h"

namespace po = boost::program_options;
#include "parameters.h"

// Check if the parameters are valid. Return 0 if they are, 1 otherwise.
int Parameters::check() const {
	if (nbSitesX <= 0) {
		std::cerr << "Error: nbSitesX should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (nbSitesY <= 0) {
		std::cerr << "Error: nbSitesY should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (nbParticles <= 0) {
		std::cerr << "Error: nbParticles should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (duration <= 0.) {
		std::cerr << "Error: duration should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (dt <= 0.) {
		std::cerr << "Error: dt should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (nbSimuls <= 0) {
		std::cerr << "Error: nbSimuls should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (probTP < 0. || probTP > 1.) {
		std::cerr << "Error: the probabilities should be between 0"
			<< " and 1." << std::endl;
		return 1;
	}
	return 0;
}

// Print the parameters to stream.
void Parameters::print(std::ostream &stream) const {
	stream << "sitesX=" << nbSitesX << ", sitesY=" << nbSitesY
		<< ", particles=" << nbParticles
		<< ", duration=" << duration << ", simuls=" << nbSimuls
		<< ", prob=" << probTP << ", determ=" << determ;
}

int Parameters::fromCommandLine(int argc, char **argv) {
	po::options_description opts("Options");
	opts.add_options()
		("sitesX,n", po::value<long>(&nbSitesX)->required(),
		 "Number of sites in x direction")
		("sitesY,m", po::value<long>(&nbSitesY)->required(),
		 "Number of sites in y direction")
		("particles,k", po::value<long>(&nbParticles)->required(),
		 "Number of particles")
		("duration,T", po::value<double>(&duration)->required(),
		 "Duration of the simulation")
		("dt,t", po::value<double>(&dt)->required(),
		 "Timestep for export")
		("simuls,s", po::value<long>(&nbSimuls)->required(),
		 "Number of repetitions of the simulation")
		("prob,p", po::value<double>(&probTP)->default_value(
		  DEFAULT_PROBA_RIGHT),
		 "Probability for the TP to jump to the right.")
        ("determ,D", po::bool_switch(&determ),
		 "Deterministic initial conditions")
		("threads,c", po::value<int>(&nbThreads)->default_value(
			DEFAULT_THREADS), "Number of threads")
		("output,o", po::value<std::string>(&output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")
        ("visu,V", po::bool_switch(&visu), "Visualisation")
        ("sleep", po::value<int>(&sleep)->default_value(DEFAULT_VISU_SLEEP),
		 "Number of ms between visualizations")
        ("verbose,v", po::bool_switch(&verbose), "Verbose mode")
        ("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

        // Display help and exit
        if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
            return 1;
        }

        po::notify(vars);

		nbSteps = (long) (duration / dt); // Number of timesteps
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	return 0;
}
