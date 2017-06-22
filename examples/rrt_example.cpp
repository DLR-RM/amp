/**
 * \brief Example using the rrt planner with inverse kinematics planner as edge steering
 * method.
 */
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <iomanip>

#include "amp/config.hpp"
#include "amp/nodes.hpp"
#include "amp/trajectories.hpp"
#include "amp/interpolate.hpp"
#include "amp/simulation.hpp"
#include "amp/planner.hpp"
#include "amp/rrt.hpp"
#include "utils/read_write_matrices.hpp"


/** Writes some information about the solution to a log file */
void print_log_file(const int sim_idx, const int iter, const int edge, const bool success, 
		const std::chrono::milliseconds& times, std::ofstream& output, const int num_try);


int main (int argc, char* argv[])
{
    typedef amp::parameters<double> config_t;
    typedef amp::node<double> node_t;
    typedef amp::trajectory<double> traj_t;
    typedef amp::planner<double> planner_t;
    typedef amp::simulation<double> sim_t;
    typedef amp::interpolator<double> interp_t;
    typedef amp::rrt<double> rrt_t;

    /** Optionally Load all of the heuristic parameters from config file */
    const char* config_file;
    if (argc>1) config_file = argv[1];
    config_t configuration(config_file);
    configuration.set_from_yaml();

    /** Create some local objects for the desired initial joint states and tcp poses */
    std::shared_ptr<sim_t> sim = configuration.create_simulation();
    const int ndofj = sim->get_model_dof();
    const int nder = 3;
    sim_t::matrix_t tcp_0 = sim_t::matrix_t::Zero(7,1);
    sim_t::matrix_t tcp_f = sim_t::matrix_t::Zero(7,1);
    sim_t::vector_t state_f = sim_t::vector_t::Zero(ndofj);
    sim_t::vector_t state_0 = sim_t::vector_t::Zero(ndofj);
    sim_t::vector_t dstate_0 = sim_t::vector_t::Zero(ndofj);

    /** Initial Configuration */
    state_0 <<  1.5, 1.3, 1.5,  1.6, 0., -0.8, 0.;

	/** Final Configuration */
    state_f << -1.5, 1.3, 1.5, -1.6, 0.,  0.8, 0.;


    /** Run a fixed number of simulations until a solution is found */
    std::ofstream log;
    log.open("./output_log.dat", std::ofstream::out);
    for (int nsim=0; nsim<configuration.rrt_max_simulations(); ++nsim) {
		std::cerr << "Simulation: " << nsim << std::endl;
		auto start = std::chrono::steady_clock::now();
		
        for (int ntry=0; ntry<configuration.rrt_max_seeds(); ++ntry) {
			std::cerr << "Try: " << ntry << std::endl;

            /** Setup trajectory and planner */
            sim->fwd_kin(state_0,dstate_0);
            tcp_0 = sim->tcp_pose(state_0);
			std::shared_ptr<node_t> node_initial = std::make_shared<node_t>(ndofj,nder);
			node_initial->time_curr = 0.;
			node_initial->state.block(0,0,7,1) = tcp_0;
			node_initial->joints.block(0,0,ndofj,1) = state_0;

            sim->fwd_kin(state_f,dstate_0);
            tcp_f = sim->tcp_pose(state_f);
			std::shared_ptr<node_t> node_final = std::make_shared<node_t>(ndofj,nder);
			node_final->time_curr = 1.;
			node_final->state.block(0,0,7,1) = tcp_f;

            /** Setup all the required amp modules */
            std::shared_ptr<interp_t> interpolator =
                    configuration.create_interpolator(ndofj, nder);

            std::shared_ptr<planner_t> local_planner =
                    configuration.create_planner(sim, interpolator);

            std::shared_ptr<rrt_t> rrt_planner =
                    configuration.create_rrt(local_planner);

            rrt_planner->set_goal(node_final);

            /** Run the Global Planner */
            bool connect_to_goal_on_first_=false;
            if (ntry==0) connect_to_goal_on_first_=true;
            rrt_t::rrt_result res_ =
                    rrt_planner->run(node_initial, configuration.rrt_max_iterations(),
                                     connect_to_goal_on_first_);

            /** Parse out the solution from RRT tree */
            if (res_.goal_connected) {
                rrt_planner->parse_solutions();
                const std::shared_ptr<const traj_t> solution = rrt_planner->get_solution(0);
                rrt_planner->write_traj_to_file(solution, "solution", nsim);
				auto end = std::chrono::steady_clock::now();
				auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
				print_log_file(nsim, res_.iterations, res_.edges, res_.goal_connected, diff, log, ntry);
				break;
			}
		}
    }

    log.close();
    return(0);
}


void print_log_file(const int sim_idx, const int iter, const int edge, const bool success, 
		const std::chrono::milliseconds& times, std::ofstream& output, const int num_try)
{
    std::ostringstream line;
    line.setf(std::ios::fixed, std::ios::floatfield);
    line.precision(4);
    line << std::setw(12) << sim_idx << " ";
    line << std::setw(12) << iter << " ";
    line << std::setw(12) << edge << " ";
    line << std::setw(12) << success << " ";
    line << std::setw(12) << times.count() << " ";
    line << std::setw(12) << num_try << " ";
    std::string strs = line.str();
    output << strs << std::endl;
}
