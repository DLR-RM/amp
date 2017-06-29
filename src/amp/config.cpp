/*
FILE amp/config.cpp

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

#include "amp/config.hpp"

namespace amp {

template<typename scalar>
parameters<scalar>::parameters(const char* config_in) :
    params(YAML::LoadFile(config_in)),
    m_isloaded(false),
    tcp_parent_joint_(0),
    tcp_translation_(std::vector<scalar>(3,0.)),
    tcp_rotation_(std::vector<scalar>(4,0.))
{
    if (!params.IsNull()) m_isloaded = true;
}

template<typename scalar>
bool parameters<scalar>::set_from_yaml()
{
    if (!m_isloaded) return(false);

    YAML::Node model = params["robot_urdf_model"];
    if (!model.IsNull()) {
        robot_model_ = model.as<std::string>();
    }

    /* Get the local trajectory planner parameters */
    YAML::Node local_planner = params["local_trajectory_planner"];
    if (!local_planner.IsNull()) {
        traj_planner.m_initial_discretization = local_planner["initial_discretization"].as<int>();

        traj_planner.m_num_fixed_bcs = local_planner["num_fixed_boundary_conditions"].as<int>();

        traj_planner.m_duration = local_planner["duration"].as<scalar>();

        YAML::Node inv_kin = local_planner["inverse_kinematics"];
        if (!inv_kin.IsNull())
            traj_planner.m_linear_solver = inv_kin["solver"].as<int>();

        YAML::Node integration = local_planner["integration"];
        if (!integration.IsNull()) {
            traj_planner.m_int_absolute_tolerance = integration["absolute_tolerance"].as<scalar>();
            traj_planner.m_int_relative_tolerance = integration["relative_tolerance"].as<scalar>();
        }

        YAML::Node pos_error = local_planner["position_error"];
        if (!pos_error.IsNull()) {
            traj_planner.m_translation_error_gain = pos_error["translation_gain"].as<scalar>();
            traj_planner.m_rotation_error_gain = pos_error["rotation_gain"].as<scalar>();
            traj_planner.m_translation_tolerance = pos_error["translation_tolerance"].as<scalar>();
            traj_planner.m_rotation_tolerance = pos_error["rotation_tolerance"].as<scalar>();
        }

        YAML::Node vel_error = local_planner["velocity_error"];
        if (!vel_error.IsNull()) {
            traj_planner.m_translation_vel_tolerance = vel_error["translation_tolerance"].as<scalar>();
            traj_planner.m_angular_vel_tolerance = vel_error["angular_tolerance"].as<scalar>();
        }
    }

    /* Get the global planner parameters */
    YAML::Node global_planner = params["global_planner"];
    if (!global_planner.IsNull()) {
        rrt_planner.m_max_seeds = global_planner["max_seeds"].as<int>();
        rrt_planner.m_max_iterations = global_planner["max_iterations"].as<int>();

        YAML::Node bounds = global_planner["boundary"];
        if (!bounds.IsNull()) {
            unsigned int bounds_size = bounds["lower"].size();
            for (unsigned int idx=0; idx<bounds_size; ++idx) {
                rrt_planner.m_bounds_lower.push_back(bounds["lower"][idx].as<scalar>());
                rrt_planner.m_bounds_upper.push_back(bounds["upper"][idx].as<scalar>());
            }
        }

        YAML::Node radius = global_planner["radius"];
        if (!radius.IsNull()) {
            rrt_planner.m_rrt_max_translation_radius = radius["translation_max_radius"].as<scalar>();
            rrt_planner.m_rrt_max_rotation_radius = radius["rotation_max_radius"].as<scalar>();
        }
    }
    return(true);
}


template<typename scalar>
std::shared_ptr<typename parameters<scalar>::sim_t>
parameters<scalar>::create_simulation()
{
    std::shared_ptr<sim_rbdl_t> sim = std::make_shared<sim_rbdl_t>(robot_model_);

    /* Get the parent joint index, otherwise default to last joint in chain */
    YAML::Node tcp_parent = params["tcp_parent_joint"];
    if (!tcp_parent.IsNull()) {
        tcp_parent_joint_ = tcp_parent.as<int>();
    }
    else tcp_parent_joint_ = sim->get_model_dof();

    /* Get the tcp transformation from parent joint, otherwise default to frame
     * of the parent joint */
    YAML::Node tcp_trafo = params["tcp_transformation_from_parent"];
    if(!tcp_trafo.IsNull()) {
        for (unsigned int idx=0; idx<tcp_trafo["translation"].size(); ++idx)
            tcp_translation_[idx] = tcp_trafo["translation"][idx].as<scalar>();
        for (unsigned int idx=0; idx<tcp_trafo["rotation"].size(); ++idx)
            tcp_rotation_[idx] = tcp_trafo["rotation"][idx].as<scalar>();
    }

    sim->set_error_gains(traj_planner.m_translation_error_gain, traj_planner.m_rotation_error_gain);
    sim->set_tcp_id(tcp_parent_joint_);
    sim->set_tcp_offset(&(tcp_translation_[0]),&(tcp_rotation_[0]));
    return(sim);
}


template<typename scalar>
std::shared_ptr<typename parameters<scalar>::interp_t>
parameters<scalar>::create_interpolator(const int n_dof, const int n_der)
{
    std::shared_ptr<interp_t> interpolate =
            std::make_shared<interp_t>(traj_planner.m_initial_discretization, n_dof, n_der);
    interpolate->set_num_fixed_bcs(traj_planner.m_num_fixed_bcs);
    return(interpolate);
}


template<typename scalar>
std::shared_ptr<typename parameters<scalar>::planner_t>
parameters<scalar>::create_planner(const std::shared_ptr<sim_t>& sim_in,
                                   const std::shared_ptr<interp_t>& interp_in)
{
    std::shared_ptr<planner_t> local_planner =
            std::make_shared<planner_t>(static_cast<typename sim_t::lin_alg_t>(traj_planner.m_linear_solver),
                                        traj_planner.m_translation_tolerance,
                                        traj_planner.m_rotation_tolerance,
                                        traj_planner.m_translation_vel_tolerance);
    local_planner->set_simulator(sim_in);
    local_planner->set_interpolator(interp_in);
    local_planner->set_observers();
    local_planner->set_discretization(traj_planner.m_initial_discretization,traj_planner.m_duration);
    local_planner->set_integration_accuracies(traj_planner.m_int_absolute_tolerance,
                                              traj_planner.m_int_relative_tolerance);
    return(local_planner);
}


template<typename scalar>
std::shared_ptr<typename parameters<scalar>::rrt_t>
parameters<scalar>::create_rrt(const std::shared_ptr<planner_t>& planner_in)
{
    std::shared_ptr<rrt_t> global_planner = std::make_shared<rrt_t>(planner_in,
                                                                    rrt_planner.m_rrt_max_translation_radius,
                                                                    rrt_planner.m_rrt_max_rotation_radius);
    global_planner->set_algorithm(rrt_t::RRT);
    return(global_planner);
}


template class parameters<double>;
}

