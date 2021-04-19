#include <exception>
#include <boost/program_options.hpp>
#include "parameters.h"

namespace po = boost::program_options;
#include "parameters.h"

// Check if the parameters are valid. Return 0 if they are, 1 otherwise.
int Parameters::check() const {
	if (nbSites <= 0) {
		std::cerr << "Error: nbSites should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (rhoP < 0 || rhoP > 1) {
		std::cerr << "Error: rhoP should be between 0 and 1."
			<< std::endl;
		return 1;
	}
	if (rhoM < 0 || rhoM > 1) {
		std::cerr << "Error: rhoM should be between 0 and 1."
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
	if (determ) {
		if ((rhoM <= 0.5 && std::floor(1/rhoM) != 1/rhoM)
			|| (rhoM > 0.5 &&
				std::abs(std::round(1/(1-rhoM))-1/(1-rhoM)) > 1e-6)) {
			std::cerr << "Error: 1/rhoM or 1/(1-rhoM) should be an integer"
				<< std::endl;
			return 1;
		}
		if ((rhoP <= 0.5 && std::floor(1/rhoP) != 1/rhoP)
			|| (rhoP > 0.5 &&
				std::abs(std::round(1/(1-rhoP))-1/(1-rhoP)) > 1e-6)) {
			std::cerr << "Error: 1/rhoP or 1/(1-rhoP) should be an integer"
				<< std::endl;
			return 1;
		}
	}
	return 0;
}

// Print the parameters to stream.
void Parameters::print(std::ostream &stream) const {
	stream << "sites=" << nbSites << ", rhoP=" << rhoP << ", rhoM=" << rhoM
		<< ", duration=" << duration << ", simuls=" << nbSimuls;
	if (rev)
		stream << ", rev";
	if (determ)
		stream << ", determ";
}

int Parameters::fromCommandLine(int argc, char **argv) {
	po::options_description opts("Options");
	opts.add_options()
		("sites,n", po::value<long>(&nbSites)->required(), "Number of sites")
		("rhoP", po::value<double>(&rhoP)->required(), "Density in front")
		("rhoM", po::value<double>(&rhoM)->required(), "Density behind")
		("duration,T", po::value<double>(&duration)->required(),
		 "Duration of the simulation")
		("dt,t", po::value<double>(&dt)->required(),
		 "Timestep for export")
        ("determ", po::bool_switch(&determ),
		 "Deterministic initial conditions")
        ("prof", po::bool_switch(&computeProfs),
		 "Compute generalized profiles")
		("sitesProf", po::value<long>(&nbSitesProf)->default_value(0),
		 "Number of sites for profiles")
        ("rev", po::bool_switch(&rev), "Compute (1-eta_1) eta_r, etc.")
		("simuls,s", po::value<long>(&nbSimuls)->required(),
		 "Number of repetitions of the simulation")
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

		nbSteps = (long) (duration / dt);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	return 0;
}
