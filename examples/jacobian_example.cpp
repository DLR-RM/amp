/**
 * \brief Test inverse kinematics planner.
 */

#include <Eigen/Dense>
#include <random>

#include "amp/config.hpp"
#include "amp/simulation.hpp"
#include "amp/nodes.hpp"

#include <chrono>
#include <iostream>
#include <string>
#include <memory>

#ifdef DEBUG
#define _DEBUG_
#endif

int main (int argc, char* argv[])
{
    typedef amp::parameters<double> config_t;
    typedef amp::node<double> node_t;
    typedef amp::simulation<double> sim_t;

    /** Optionally Load all of the heuristic parameters from config file */
    const char* config_file;
    if (argc>1) config_file = argv[1];
    config_t configuration(config_file);
    configuration.set_from_yaml();

    /** Create some local objects for the desired initial joint states and tcp poses */
    std::shared_ptr<sim_t> sim = configuration.create_simulation();
    const int ndofj = sim->get_model_dof();
    sim_t::vector_t state = sim_t::vector_t::Zero(ndofj);
    sim_t::jacobian_t jacob;

    /** Setup rng for joint values between -2 and 2 radians */
    typedef std::mt19937_64 gen_t;
    typedef std::uniform_real_distribution<double> uni_t;
	gen_t generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	uni_t dist(-2.,2.);

    /** Setup rng for joint values between -2 and 2 radians */
	auto start = std::chrono::steady_clock::now();
    for (int idx = 0; idx<1e8; ++idx) {
		for (int idxj=0; idxj<ndofj; ++idxj) {
			state(idxj) = dist(generator);
		}	
        jacob = sim->get_J_with_update(state);
	}
    auto end = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cerr << "Computation time: " << diff.count() << " milliseconds.\n";

    return(0);
}

