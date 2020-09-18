/*
Copyright Université Pierre et Marie Curie (2017)
Contributors: Alexis Poncet

alexis.poncet@ens.fr

This software is a computer program whose purpose is to simulate the
motion of interacting Brownian particles in a quasi-1d environment,
some of them being driven by an external force.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
/*
 * BrownianQuasi1d
 * SimulCanal.h
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Definitions related to the simulation in a canal.
 */

#ifndef SIMUL_CANAL_H
#define SIMUL_CANAL_H

#include "parameters.h"
#include "Simul.h"

#define DIM_CANAL 2

class SimulCanal : public Simulation {
	public:
		SimulCanal(const Parameters &p);
		~SimulCanal();

	protected:
		// Positions of the particles
		std::vector< std::array<double, DIM_CANAL> > positions;
		// Forces between the particles
		std::vector< std::array<double, DIM_CANAL> > forces;

		// Methods to implement from Simulation
		int init(std::mt19937 &rndGen) override;
		void update(std::mt19937 &rndGen, const bool thermalization) override;
#ifdef USE_MKL
		int init(VSLStreamStatePtr stream) override;
		void update(VSLStreamStatePtr stream,
				    const bool thermalization = false) override;
#endif
		double getPosX(const long i) override;
		void getPosRel(std::vector<double> &posr, const long nbPts) override;

		// New methods
		void calcForcesBetweenParticles();
		void keepInCanal();
};

#endif
