/*
FILE amp/simulation_rbdl.hpp

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

#ifndef _AMP_RBDL_INTERFACE_HPP_
#define _AMP_RBDL_INTERFACE_HPP_

#include <iostream>
#include <string>
#include <cmath>

#ifdef USE_EIGEN_PKGCONFIG_HEADERS
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif

#include "amp/simulation.hpp"

#include <rbdl/rbdl.h>
#include <rbdl/addons/urdfreader/urdfreader.h>
#include <rbdl/Logging.h>

#ifdef DEBUG
#define _DEBUG_
#endif

namespace amp {


template<typename scalar>
class rbdl_interface: public simulation<scalar>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef simulation<scalar> sim_t;
    typedef RigidBodyDynamics::Model model_t;
    typedef Eigen::Matrix<scalar,3,1> vector3_t;
    typedef Eigen::Matrix<scalar,3,3> matrix3_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,1> vector_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> matrix_t;
    typedef matrix_t jacobian_t;

    typedef Eigen::LLT<matrix_t> llt_t;
    typedef Eigen::LDLT<matrix_t> ldlt_t;
    typedef Eigen::HouseholderQR<matrix_t> qr_t;
    typedef Eigen::FullPivLU<matrix_t> lu_t;
    typedef Eigen::JacobiSVD<matrix_t> svd_t;

    rbdl_interface(const std::string& filename);
    
    ~rbdl_interface();

    const int get_model_dof() const;

    /**
     * @brief Runs the forward kinematics, simplify accelerations to zero,
     * model accelerations, angular accelerations, corioli will not be
     * calculated.
     */
    void fwd_kin(const vector_t& state_, const vector_t& dstate_);

    /**
     * @brief Runs the forward kinematics
     */
    void fwd_dyn(vector_t& ddstate_, const vector_t& state_,
                 const vector_t& dstate_, const vector_t& tau_);

    /**
     * @brief Runs the forward kinematics and sets current
     * tcp pose
     */
    matrix_t tcp_pose(const vector_t& state_);

    matrix3_t rotation_to_tcp_frame(const vector_t& state_);

    matrix3_t rotation_to_base_frame(const vector_t& state_);

    /**
     * @brief Computes the manipulability
     * WARNING: performs the determinant on J
     */
    scalar manipulability(const jacobian_t& J_) const;

    /**
     * @brief Sets the current Jacobian
     */
    jacobian_t get_J(const vector_t& state_) const;

    /**
     * @brief Sets the current Jacobian
     */
    jacobian_t get_J_with_update(const vector_t& state_) const;

    /**
     * @brief solves dstate: dtcp = J(state)*dstate
     * using one of the available linear solver
     */
    vector_t inv_kin(const matrix_t& dtcp_, const vector_t& state_,
                     const typename sim_t::lin_alg_t lin_=sim_t::DEFAULT_ALG);

    /**
     * @brief Difference beteween desired tcp pose \a tcp_ and actual tcp
     * pose from forward kinematics, \a tcp_e
     */
    matrix_t tcp_error(const matrix_t& tcp_, const vector_t& state_);

    /**
     * @brief Composes the 6x1 gain * tcp-position-error term used to
     * correct for inverse kinematics error (see Siciliano).
     */
    matrix_t tcp_error_correction(const matrix_t& tcp_e);

    /**
     * @brief integrate states forward one time step
     */
    vector_t diff_kin(const vector_t& state_, const vector_t& dstate_,
                      scalar dt_);

    /**
     * @brief runs the forward kinematics using numerically computed dstate
     */
    matrix_t inv_kin_error_check(const vector_t& state_, const vector_t& dstate_);

private:
    model_t* model;

};

} // namespace amp
#endif // _AMP_RBDL_INTERFACE_HPP_
