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
	if (nbTracers > nbParticles) {
		std::cerr << "Error: The number of tracers should be smaller than"
			<< " the number of particles." << std::endl;
		return 1;
	}
	if ((long) initPos.size() != nbTracers) {
		std::cerr << "Critical error: the size of the vector of positions"
			<< " is wrong." << std::endl;
		return 1;
	}
	for (auto po : initPos) {
		if (po < 0 || po >= nbSites) {
			std::cerr << "Error: the positions should in the right range."
				<< std::endl;
			return 1;
		}
	}
	for (long i = 0 ; i < nbTracers - 1 ; ++i) {
		if (initPos[i] >= initPos[i+1]) {
			std::cerr << "Error: please sort the positions."
				<< std::endl;
			return 1;
		}
	}
	if ((long) probas.size() != nbTracers) {
		std::cerr << "Error: the number of probabilities and the number"
			<< " of tracers should be equal." << std::endl;
		return 1;
	}
	for (auto pr : probas) {
		if (pr < 0. || pr > 1.) {
			std::cerr << "Error: the probabilities should be between 0"
				<< " and 1." << std::endl;
			return 1;
		}
	}
	return 0;
}

// Print the parameters to stream.
void Parameters::print(std::ostream &stream) const {
	stream << "sites=" << nbSites << ", particles=" << nbParticles
		<< ", tracers=" << nbTracers << ", duration=" << duration
		<< ", simuls=" << nbSimuls;
	stream << ", initPos=";
	for (auto po : initPos) {
		stream << po << ":";
	}
	stream << ", probas=";
	for (auto pr : probas) {
		stream << pr << ":";
	}
}

int Parameters::fromCommandLine(int argc, char **argv) {
	po::options_description opts("Options");
	opts.add_options()
		("sites,n", po::value<long>(&nbSites)->required(), "Number of sites")
		("particles,m", po::value<long>(&nbParticles)->required(),
		 "Number of particles")
		("duration,T", po::value<double>(&duration)->required(),
		 "Duration of the simulation")
		("dt,t", po::value<double>(&dt)->required(),
		 "Timestep for export")
		("simuls,s", po::value<long>(&nbSimuls)->required(),
		 "Number of repetitions of the simulation")
		("initPos,d",
		 po::value< std::vector<long> >(&initPos)->multitoken()->required(),
		 "Initial positions of the tracers. It fixes the number of tracers.")
		("probas,p",
		 po::value< std::vector<double> >(&probas)->multitoken()->required(),
		 "Probabilities to jump to the right.")
		("threads,c", po::value<int>(&nbThreads)->default_value(
			DEFAULT_THREADS), "Number of threads")
		("output,o", po::value<std::string>(&output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")
		("custom", po::value<std::string>(&custom)->default_value(
			DEFAULT_CUSTOM_FILE), "File for custom observables")
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

		// The number of tracers is given by the number of positions
		nbTracers = (long) initPos.size();
		// Number of timesteps
		nbSteps = (long) (duration / dt);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	return 0;
}
