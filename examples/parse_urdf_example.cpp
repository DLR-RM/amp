/**
 * \brief Example loading a urdf model.
 */

#include <Eigen/Dense>

#include "amp/config.hpp"
#include "amp/simulation.hpp"

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
    typedef amp::simulation<double> sim_t;

    auto start = std::chrono::steady_clock::now();

    /** Optionally Load all of the heuristic parameters from config file */
    const char* config_file;
    if (argc>1) config_file = argv[1];
    config_t configuration(config_file);
    configuration.set_from_yaml();

    /** Create some local objects for the desired initial joint states and tcp poses */
    std::shared_ptr<sim_t> sim = configuration.create_simulation();

#ifdef _DEBUG_
    std::cerr << "model num dof: " << sim->get_model_dof() << std::endl;
#endif // _DEBUG_

    auto end = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cerr << "Computation time: " << diff.count() << " milliseconds.\n";

    return 0;
}
