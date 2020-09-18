#include <iostream>
#include "parameters.h"
#include "simul.h"

int main(int argc, char **argv) {
	// Load parameters
    Parameters params;
    int status = params.fromCommandLine(argc, argv);
	if (status) {
		return status - 1;
	}
	status = params.check();
	if (status) {
		return status;
	}

    if (params.verbose) {
		std::cout << "--- This is " << argv[0] << ", compiled on " << __DATE__
			<< " at " << __TIME__ << " ---" << std::endl;
        params.print();
		std::cout << std::endl;
    }

	status = runSimulations(params);

	return status;
}
