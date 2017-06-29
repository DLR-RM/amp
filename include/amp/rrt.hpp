/*
FILE amp/rrt.hpp

LICENSE

Copyright 2017 Samantha Stoneman
Research Associate, German Aerospace Center (DLR)

This file is part of the Articulated-robot Motion Planner (Amp) library.

Amp is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Amp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Amp. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _AMP_RRT_HPP_
#define _AMP_RRT_HPP_

#include <iostream>
#include <memory>
#include <chrono>
#include <list>
#include <utility>
#include <random>
#include <limits>


#ifdef USE_EIGEN_PKGCONFIG_HEADERS
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#else
#include <Eigen/Dense>
#include <Eigen/StdVector>
#endif

#include "kdtree/kdtree.h"

#include "amp/nodes.hpp"
#include "amp/trajectories.hpp"
#include "amp/planner.hpp"
#include "utils/read_write_matrices.hpp"


#ifdef DEBUG
#define _DEBUG_
#endif

#ifdef VERBOSE
#define _VERBOSE_
#endif

#ifdef LOG_ALL_DATA
#define _LOG_ALL_DATA_
#endif

namespace amp {

/**
 * @class branches
 * Holds pointer to a node and its incoming and outgoing
 * trajectories. From any branch object the tree can be parsed forwards and
 * backwards.
 */
template<typename scalar>
class branches {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef node<scalar> node_t;
    typedef trajectory<scalar> traj_t;

    // is that necessary since they are pointers?
    typedef std::list<std::shared_ptr<const traj_t> > connections_t;
    typedef std::pair<const branches<scalar>*, typename connections_t::const_iterator> incoming_t;


    /**
     * @brief node_
     */
    const std::shared_ptr<const node_t> node_;

    /**
     * @brief outgoing_branches_
     */
    connections_t outgoing_branches_;

    /**
     * @brief incoming_branches
     */
    const incoming_t incoming_branches;

    /**
     * @brief branches Instantiates a new branch and sets its preceding branch
     * @param node_in
     * @param branch_in
     * @param traj_in
     */
    branches(const std::shared_ptr<const node_t>& node_in,
             const branches<scalar>* branch_in,
             typename connections_t::const_iterator traj_in) :
        node_(node_in),
        incoming_branches(std::make_pair(branch_in,traj_in)) { }

    /**
     * @brief branches Instantiates a new branch with no predecessor. This is a
     * convenience overload for the root node.
     * @param node_in
     */
    branches(const std::shared_ptr<const node_t>& node_in) :
        node_(node_in) { }

    /**
     * @brief insert_branch
     * @param traj_in
     */
    void insert_branch(const std::shared_ptr<const traj_t>& traj_in) {
        outgoing_branches_.push_back(traj_in);
    }

    /**
     * @brief remove_branch
     * @param traj_in
     */
    void remove_branch(const std::shared_ptr<const traj_t>& traj_in) {
        outgoing_branches_.remove(traj_in); // calls destructor of
    }
};





/**
 * @class rrt Rapidly exploring random trees class which owns the graph mechanism
 * branches and methods to interact with the graph.
 * rrt can act as a traditional rrt, rrt* or rrt*-gbo.
 *
 */
template<typename scalar>
class rrt {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef node<scalar> node_t;
    typedef trajectory<scalar> traj_t;
    typedef planner<scalar> planner_t;
    typedef branches<scalar> branches_t;
    typedef std::list<std::shared_ptr<branches_t> > tree_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,1> vector_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> matrix_t;
    typedef std::mt19937_64 gen_t;
    typedef std::uniform_real_distribution<scalar> uni_t;
    typedef rwmatrix<scalar,matrix_t> rwm;
    typedef std::list<std::shared_ptr<const node_t> > sample_list_t;
    typedef std::list<std::shared_ptr<const traj_t> > traj_list_t;
    typedef std::list<std::shared_ptr<const branches_t> > goal_list_t;

    struct rrt_result {
        bool goal_connected;
        unsigned int iterations;
        unsigned int edges;

        rrt_result() :
            goal_connected(false),
            iterations(0),
            edges(0) { }

        rrt_result(const rrt_result& res_) :
            goal_connected(res_.goal_connected),
            iterations(res_.iterations),
            edges(res_.edges) { }
    };

    /**
     * @brief The search_type_t enum
     */
    enum search_type_t {
        DEFAULT = -1,
        RRT = 0,
        RRTSTAR,
        RRTSTARGBO
    };

    /**
     * @brief rrt
     * @param planner_in
     */
    rrt(const std::shared_ptr<planner_t>& planner_in,
        scalar r_radius_max_in=std::numeric_limits<scalar>::max(),
        scalar q_radius_max_in=std::numeric_limits<scalar>::max());

    /**
     *
     */
    ~rrt();


    /**
     * @brief run Run the search algorithm and fill in the tree graph
     */
    rrt_result run(const std::shared_ptr<const node_t>& root_, const int stop_iter, const bool connect_to_goal_first=true);

    /**
     * @brief set_algorithm
     */
    void set_algorithm(const search_type_t search_type_in) {
        search_type = search_type_in;
    }

    /**
     * @brief set_algorithm
     */
    void set_goal(const std::shared_ptr<node_t>& goal_in) {
        goal_ = goal_in;
    }

    /**
     * @brief get_solution TODO: fix error case, current empty pointer return
     * isn't safe.
     * @param sol_idx
     * @return
     */
    const std::shared_ptr<const traj_t> get_solution(const int sol_idx) const {
        if (sol_idx==0)
            return(trajs_solution.front());
        else if (sol_idx<trajs_solution.size()) {
            auto traj=trajs_solution.begin();
            for(int i=1; i<sol_idx; ++i)
                ++traj;
            return(*traj);
        }
        else
            return(std::shared_ptr<traj_t>());
    }

    /**
     * @brief parse_tree
     */
    void parse_tree();

    /**
     * @brief parse_solutions
     */
    void parse_solutions();

    /**
     * @brief write_traj_to_file
     * @param traj_
     * @param filename
     */
    void write_traj_to_file(const std::shared_ptr<const traj_t>& traj_,
                            const char* filename, const int idx);

private:    
    /**
     * @brief planner_
     */
    const std::shared_ptr<planner_t> planner_;

    tree_t tree_;

    kdtree* kdtree_;

    search_type_t search_type;

    std::shared_ptr<node_t> goal_;

    sample_list_t samples_failed;
    traj_list_t trajs_failed;

    sample_list_t samples_succeeded;
    traj_list_t trajs_succeeded;

    traj_list_t trajs_solution;
    goal_list_t solutions;


    /**
     * @brief generator
     */
    gen_t generator;

    /**
     * @brief distributions
     */
    std::vector<uni_t> distributions;

    /**
     * @brief radius_max
     */
    Eigen::Matrix<scalar,2,1> radius_max;

    /**
     * @brief seed_rng
     */
    void seed_rng();

    /**
     * @brief shorten_sample_distance
     * @param sample_in
     */
    void shorten_sample_distance(const std::shared_ptr<node_t>& sample_in,
                                 const std::shared_ptr<const node_t>& nearest_in);

    /**
     * @brief create_sample
     * @return
     */
    std::shared_ptr<node_t> create_sample();

    /**
     * @brief create_branch
     * @return
     */
    std::shared_ptr<branches_t> create_branch(const std::shared_ptr<const node_t> &root_in);

    /**
     * @brief insert_node
     * @param branch_in
     */
    void insert_node(const std::shared_ptr<branches_t>& branch_in);

    /**
     * @brief nearest_node
     * @param node_in
     * @return
     */
    branches_t *nearest_node(const std::shared_ptr<node_t> &node_in);


    /**
     * @brief nearest_nodes_in_radius
     * @param node_in
     * @return
     */
    tree_t nearest_nodes_in_radius(std::shared_ptr<node_t> node_in, const scalar);

    /**
     * @brief connect
     * @param traj_result
     * @param node_to
     * @param node_from
     * @return
     */
    typename planner_t::return_type_t connect(std::shared_ptr<traj_t> &traj_result,
                                              const std::shared_ptr<const node_t> &node_from,
                                              const std::shared_ptr<const node_t> &node_to);

    /**
     * @brief log_failed
     * @param sample_failed
     * @param traj_failed
     */
    void log_failed(const std::shared_ptr<const node_t>& sample_failed,
                    const std::shared_ptr<const traj_t>& traj_failed);

    /**
     * @brief write_data_to_file
     */
    void write_data_to_file(sample_list_t& samples, traj_list_t& trajs, const char* filename);

};

} //namespace amp

#endif //_AMP_RRT_HPP_
