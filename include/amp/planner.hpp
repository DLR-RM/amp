/*
FILE amp/planner.hpp

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

#ifndef _AMP_PLANNER_HPP_
#define _AMP_PLANNER_HPP_

#include <iostream>
#include <memory>
#include <cmath>

#ifdef USE_EIGEN_PKGCONFIG_HEADERS
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#else
#include <Eigen/Dense>
#include <Eigen/StdVector>
#endif

#include <boost/numeric/odeint.hpp>

#include "amp/nodes.hpp"
#include "amp/trajectories.hpp"
#include "amp/interpolate.hpp"
#include "amp/simulation.hpp"
#include "amp/integrator.hpp"
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

template <typename scalar>
class planner
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef simulation<scalar> sim_t;
    typedef rb_integrator<scalar> rbi_t;
    typedef rb_observer<scalar> rbo_t;
    typedef node<scalar> node_t;
    typedef trajectory<scalar> traj_t;
    typedef interpolator<scalar> interp_t;
    typedef typename sim_t::vector_t vector_t;
    typedef typename sim_t::matrix_t matrix_t;
    typedef observers<scalar> obs_t;
    typedef rwmatrix<scalar,matrix_t> rwm;

    /**
     * @brief The return_type_t enum
     */
    enum return_type_t {
        SUCCESS = 0,
        JOINT_LIMIT,
        JOINT_VEL_LIMIT,
        TCP_VEL_LIMIT,
        COLLISION,
        TCP_ERROR,
        TCP_VEL_ERROR,
        INT_MIN_STEPSIZE,
        INT_MAX_STEPS
    };

    /**
     * @brief planner
     * @param sim_in
     * @param interp_in
     * @param filename_
     */
    planner(typename sim_t::lin_alg_t lin_solver_in = sim_t::lin_alg_t::LU,
            const double pos_tol=1.,
            const double rot_tol=1.,
            const double vel_tol=1.);

    /**
     * @brief destructor
     */
    ~planner();

    /**
     * @brief set_simulator
     * @param sim_in
     */
    void set_simulator(const std::shared_ptr<sim_t>& sim_in);

    /**
     * @brief set_interpolator
     * @param interp_in
     */
    void set_interpolator(const std::shared_ptr<interp_t>& interp_in);

    /**
     * @brief set_observers
     */
    void set_observers();

    /**
     * @brief set_discretization
     * @param int_via_points_in
     * @param ik_via_points_in
     * @param tf_
     */
    void set_discretization(const int int_via_points_in,
                            const scalar tf_);


    /**
     * @brief set_desired_tcp_trajectory
     * @param np_
     * @param node_initial_
     * @param node_final_
     */
    void set_desired_tcp_trajectory(const int np_,
                                    const std::shared_ptr<const node_t> &node_initial_,
                                    const std::shared_ptr<const node_t> &node_final_);

    /**
     * @brief planner<scalar>::set_desired_tcp_trajectory
     * @param np_
     * @param traj_
     * @param node_initial_
     * @param node_final_
     */
    void set_desired_tcp_trajectory(const int np_,
                                    std::shared_ptr<traj_t>& traj_,
                                    const std::shared_ptr<const node_t>& node_initial_,
                                    const std::shared_ptr<const node_t>& node_final_);

    /**
     * @brief dof
     * @return
     */
    inline int dof() const { return(dof_); }

    /**
     * @brief eval_stop_conditions
     * @param m_measure
     */
    void eval_stop_conditions(const scalar m_measure = 1.);


    /**
     * @brief run constant step euler integrator
     * this is just a simple for loop!
     */
    void run(const vector_t& state_in,
             const vector_t& dstate_in,
             const bool output_states=false);



    /**
     * @brief run adaptive step integrator from boost odeint module
     */
    return_type_t run_adaptive(const vector_t& state_in,
                               const vector_t& dstate_in,
                               const bool write_data_flag=false,
                               const char *filename=NULL);


    /**
     * @brief tcp_traj
     * @return
     */
    std::shared_ptr<traj_t> tcp_traj() const {
        return(tcp_traj_);
    }


    /**
     * @brief set_integration_accuracies
     * @param eta_abs_in
     * @param eta_rel_in
     */
    void set_integration_accuracies(const scalar eta_abs_in,
                                    const scalar eta_rel_in) {
        int_eta_abs = eta_abs_in;
        int_eta_rel = eta_rel_in;
    }


    scalar state_bound;
    scalar dstate_bound;
    scalar tcp_bound;
    scalar dtcp_bound;

private:
    /**  @brief sim_ */
    std::shared_ptr<sim_t> sim_;

    /** @brief interpolate_ */
    std::shared_ptr<interp_t> interpolate_;

    /**  @brief tcp_traj_ */
    std::shared_ptr<traj_t> tcp_traj_;

    /**  @brief observer_ */
    obs_t observer_;

    /** @brief observer_states */
    obs_t observer_states;

    /** @brief rb_integrator_ */
    rbi_t* rb_integrator_;

    /**  @brief rb_observer_ */
    rbo_t* rb_observer_;

    /** @brief robot model dof */
    int dof_;

    /**  @brief use stopping condition for end-effector position error  */
    bool flag_epos = false;


    /**  @brief use stopping condition for end-effector velocity error */
    bool flag_evel = false;


    /**  @brief discretization for joint velocity integration calls */
    int int_via_points;


    /** @brief optional fixed joint in the kinematic chain */
    int fixed_joint_idx;


    /** @brief lin_solver to use for inverse kinematics computation */
    typename sim_t::lin_alg_t lin_solver;

    const double pos_tol;
    const double rot_tol;
    const double vel_tol;

    scalar tf;
    scalar dt_int;

    scalar int_eta_abs;
    scalar int_eta_rel;

    void write_data(const Eigen::Map<matrix_t> &times_in,
                    const int n_result_in,
                    const char* filename);
};

} //namespace amp

#endif // _AMP_PLANNER_HPP_
