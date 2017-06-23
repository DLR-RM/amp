/*
FILE amp/integrator.cpp

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

#include "amp/integrator.hpp"

namespace amp {


template<typename scalar>
rb_integrator<scalar>
::rb_integrator(sim_t* sim_in,
                const std::shared_ptr<interp_t> &interpolate_in,
                observers<scalar>& observers_in,
                lin_alg_t lin_solver_in,
                const int fixed_joint_in) :
    sim_(sim_in),
    interpolate_(interpolate_in),
    observers_(observers_in),
    lin_solver_(lin_solver_in),
    tcp_0(matrix_t::Zero(7,1)),
    tcp_i(observers_.tcp_i_),
    tcp_star(observers_.tcp_star_),
    dtcp_i(observers_.dtcp_i_),
    dtcp_star(observers_.dtcp_star_),
    dtcp_c(observers_.dtcp_c_),
    state_star(observers_.state_star_),
    dstate_star(observers_.dstate_star_),
    manipulability(observers_.manipulability_)
{

    tcp_i = matrix_t::Zero(7,1);
    tcp_star = matrix_t::Zero(7,1);
    dtcp_i = matrix_t::Zero(6,1);
    dtcp_star = matrix_t::Zero(6,1);
    dtcp_c = matrix_t::Zero(6,1);
    state_star = vector_t::Zero(sim_->get_model_dof());
    dstate_star = vector_t::Zero(sim_->get_model_dof());
    manipulability = 0.;
}



template<typename scalar>
void rb_integrator<scalar>
::initialize(const vector_t& state_, const vector_t& dstate_, const lin_alg_t lin_solver_in)
{
    lin_solver_ = lin_solver_in;
    tcp_i.block(0,0,3,1) = tcp_traj_->dimensions(0,0,1).block(0,0,1,3).transpose();
    tcp_i.block(3,0,4,1) = tcp_traj_->dimensions(0,0,1).block(0,3,1,4).transpose();
    dtcp_i.block(0,0,3,1) = tcp_traj_->dimensions(1,0,1).block(0,0,1,3).transpose();
    dtcp_i.block(3,0,3,1) = tcp_traj_->angular_rates(1,0,1).transpose();

    state_star = state_;
    dstate_star = dstate_;
    sim_->fwd_kin(state_star,dstate_star);
    tcp_star = sim_->tcp_pose(state_star);
    dtcp_star = sim_->inv_kin_error_check(state_star,dstate_star);

    tcp_0 = tcp_star;

#ifdef _DEBUG_
    /** compute manipulability */
    manipulability = sim_->manipulability(sim_->get_J(state_star));
#endif //_DEBUG
}



template<typename scalar>
void rb_integrator<scalar>
::set_tcp_traj(std::shared_ptr<traj_t> tcp_)
{
    tcp_traj_=tcp_;
}



template<typename scalar>
void rb_integrator<scalar>
::operator() (const state_type& y , state_type& dydt , const scalar t)
{
    MapStateConst y_map(y.data(), y.size(), 1);
    MapState dydt_map(dydt.data(), dydt.size(), 1);

    /** get curve values at current t, needed if time points not known apriori */
    interpolate_->sample_one_time(tcp_traj_,t+interpolate_->time_initial());
    tcp_traj_->normalize_quaternions();
    tcp_traj_->check_antipodal_quaternion(tcp_i.block(3,0,4,1).array());
    tcp_traj_->compute_angular_rates();

    /** set the reference point */
    tcp_i.block(0,0,3,1) = tcp_traj_->dimensions(0,0,1).block(0,0,1,3).transpose();
    tcp_i.block(3,0,4,1) = tcp_traj_->dimensions(0,0,1).block(0,3,1,4).transpose();
    dtcp_i.block(0,0,3,1) = tcp_traj_->dimensions(1,0,1).block(0,0,1,3).transpose();
    dtcp_i.block(3,0,3,1) = tcp_traj_->angular_rates(1,0,1).transpose();

    /** compute current tcp error, also runs fwd kinematics to get current Jacobian */
    tcp_star = sim_->tcp_pose(y_map);
    dtcp_c = sim_->tcp_error_correction(tcp_i-tcp_star);

    /** run inverse kinematics */
    matrix_t dtcp_i_e = dtcp_i;
    //    dtcp_i_e.block(0,0,3,1) = sim_->rotation_to_tcp_frame(y_map) * dtcp_i.block(0,0,3,1);
    dtcp_i_e = dtcp_i;
    dydt_map = sim_->inv_kin((dtcp_i_e+dtcp_c), y_map, lin_solver_);

    /** reset data, calculate velocity error */
    state_star = y_map;
    dstate_star = dydt_map;
    sim_->fwd_kin(state_star,dstate_star);
    dtcp_star = sim_->inv_kin_error_check(state_star,dstate_star);
    //    dtcp_star.block(0,0,3,1) = sim_->rotation_to_base_frame(state_star) * dtcp_star.block(0,0,3,1);


#ifdef _DEBUG_
    /** compute manipulability */
    manipulability = sim_->manipulability(sim_->get_J(state_star));
#endif //_DEBUG
}

// instantiate robot class
template class rb_integrator<double>;

}
