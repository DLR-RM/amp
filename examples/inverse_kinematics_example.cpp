/**
 * \brief Test inverse kinematics planner.
 */

#include <Eigen/Dense>

#include "amp/config.hpp"
#include "amp/simulation.hpp"
#include "amp/nodes.hpp"
#include "amp/interpolate.hpp"
#include "amp/planner.hpp"

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
    typedef amp::planner<double> planner_t;
    typedef amp::simulation<double> sim_t;
    typedef amp::interpolator<double> interp_t;

    auto start = std::chrono::steady_clock::now();

    /** Optionally Load all of the heuristic parameters from config file */
    const char* config_file;
    if (argc>1) config_file = argv[1];
    config_t configuration(config_file);
    configuration.set_from_yaml();

    /** Create some local objects for the desired initial joint states and tcp poses */
    std::shared_ptr<sim_t> sim = configuration.create_simulation();

    const int nder = 3;
    const int ndofj = sim->get_model_dof();
    sim_t::matrix_t tcp_0 = sim_t::matrix_t::Zero(7,1);
    sim_t::matrix_t tcp_f = sim_t::matrix_t::Zero(7,1);
    sim_t::vector_t state_f = sim_t::vector_t::Zero(ndofj);
    sim_t::vector_t state_0 = sim_t::vector_t::Zero(ndofj);
    sim_t::vector_t dstate_0 = sim_t::vector_t::Zero(ndofj);

    /** Initial Configuration */
    state_0 <<  1.5, 1.3, 1.5,  1.6, 0., -0.8, 0.;

    /** Setup the initial tcp pose */
    sim->fwd_kin(state_0,dstate_0);
    tcp_0 = sim->tcp_pose(state_0);
    std::shared_ptr<node_t> node_initial = std::make_shared<node_t>(ndofj,nder);
    node_initial->time_curr = 0.;
    node_initial->state.block(0,0,7,1) = tcp_0;
    node_initial->joints.block(0,0,ndofj,1) = state_0;

    /** Final tcp pose is 20 cm above initial */
    std::shared_ptr<node_t> node_final = std::make_shared<node_t>(ndofj,nder);
    node_final->time_curr = configuration.local_trajectory_duration();
    node_final->state.block(0,0,7,1) = tcp_0;
    node_final->state(2,0) = 0.2;


    /** Setup all the required amp modules */
    std::shared_ptr<interp_t> interpolator =
            configuration.create_interpolator(ndofj, nder);

    std::shared_ptr<planner_t> local_planner =
            configuration.create_planner(sim, interpolator);

    /** Use the local planner to find a trajectory using the Inverse Kinematics */
    local_planner->set_desired_tcp_trajectory(node_final->time_curr, node_initial, node_final);
    local_planner->run_adaptive(state_0, sim_t::vector_t::Zero(ndofj));

    auto end = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cerr << "Computation time: " << diff.count() << " milliseconds.\n";

    return(0);
}
