/*
FILE amp/config.hpp

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

#ifndef _AMP_CONFIG_HPP_
#define _AMP_CONFIG_HPP_

#include "yaml-cpp/yaml.h"

#include <string>
#include <iostream>
#include <memory>

#include "simulation.hpp"
#include "simulation_rbdl.hpp"
#include "amp/interpolate.hpp"
#include "amp/planner.hpp"
#include "amp/rrt.hpp"

#ifdef DEBUG
#define _DEBUG_
#endif

namespace amp {

/**
 *
 */
template<typename scalar>
class parameters {
public:
    typedef simulation<scalar> sim_t;
    typedef rbdl_interface<scalar> sim_rbdl_t;
    typedef interpolator<scalar> interp_t;
    typedef planner<scalar> planner_t;
    typedef rrt<scalar> rrt_t;

    /**
     * \brief Initializes all parameters to the default values.
     */
    parameters(const char* config_in);

    /**
     * \brief Sets the parameters defined in the yaml config_in file.
     */
    bool set_from_yaml();

    /**
     * \brief Returns the simulation pointer.
     */
    std::shared_ptr<sim_t> create_simulation();

    /**
     * @brief create_interpolator
     * @param n_dof
     * @param n_der
     * @return
     */
    std::shared_ptr<interp_t> create_interpolator(const int n_dof, const int n_der);

    /**
     * @brief create_planner
     * @param sim_in
     * @param interp_in
     * @return
     */
    std::shared_ptr<planner_t> create_planner(const std::shared_ptr<sim_t>& sim_in,
                                              const std::shared_ptr<interp_t>& interp_in);

    /**
     * @brief create_rrt
     * @param planner_in
     * @return
     */
    std::shared_ptr<rrt_t> create_rrt(const std::shared_ptr<planner_t>& planner_in);

    /**
     * @brief The local_trajectory_planner struct
     * Groups parameters pertaining to the local trajectory planner.
     */
    struct local_trajectory_planner {
        inline local_trajectory_planner() :
            m_initial_discretization(100),
            m_num_fixed_bcs(3),
            m_duration(10.),
            m_linear_solver(-1),
            m_int_absolute_tolerance(1e-8),
            m_int_relative_tolerance(1e-8),
            m_translation_error_gain(0.),
            m_rotation_error_gain(10.),
            m_translation_tolerance(1e-3),
            m_rotation_tolerance(1e-3),
            m_translation_vel_tolerance(1e-1),
            m_angular_vel_tolerance(1e-1) {}

        int m_initial_discretization;
        int m_num_fixed_bcs;
        scalar m_duration;
        int m_linear_solver;
        scalar m_int_absolute_tolerance;
        scalar m_int_relative_tolerance;
        scalar m_translation_error_gain;
        scalar m_rotation_error_gain;
        scalar m_translation_tolerance;
        scalar m_rotation_tolerance;
        scalar m_translation_vel_tolerance;
        scalar m_angular_vel_tolerance;
    };

    /**
     * @brief The global_planner struct
     * Groups the parameters pertaining to the global planner.
     */
    struct global_planner {
        global_planner():
            m_max_iterations(100),
            m_max_seeds(10),
            m_max_simulations(1),
            m_rrt_max_translation_radius(0.5),
            m_rrt_max_rotation_radius(0.5) {}

        int m_max_iterations;
        int m_max_seeds;
        int m_max_simulations;
        scalar m_rrt_max_translation_radius;
        scalar m_rrt_max_rotation_radius;
        std::vector<scalar> m_bounds_upper;
        std::vector<scalar> m_bounds_lower;
    };

    inline int local_trajectory_duration() const {return(traj_planner.m_duration);}

    inline int rrt_max_simulations() const {return(rrt_planner.m_max_simulations);}

    inline int rrt_max_iterations() const {return(rrt_planner.m_max_iterations);}

    inline int rrt_max_seeds() const {return(rrt_planner.m_max_seeds);}

private:
    YAML::Node params;
    bool m_isloaded;
    local_trajectory_planner traj_planner;
    global_planner rrt_planner;
    std::string robot_model_;
    std::string kd_lib_;
    int tcp_parent_joint_;
    std::vector<scalar> tcp_translation_;
    std::vector<scalar> tcp_rotation_;
};

}
#endif //_AMP_CONFIG_HPP_

















