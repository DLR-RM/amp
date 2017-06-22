/*
FILE amp/planner.cpp

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

#include "amp/planner.hpp"

int global_count = 0;

namespace amp {

template<typename scalar>
planner<scalar>
::planner(typename sim_t::lin_alg_t lin_solver_in,
          const double pos_tol_in,
          const double rot_tol_in,
          const double vel_tol_in) :
    state_bound(120.*M_PI/(180.)),
    dstate_bound(120.*M_PI/(180.)),
    tcp_bound(1.),
    dtcp_bound(1.),
    dof_(0),
    flag_epos(false),
    flag_evel(false),
    int_via_points(0),
    fixed_joint_idx(-1),
    lin_solver(lin_solver_in),
    pos_tol(pos_tol_in),
    rot_tol(rot_tol_in),
    vel_tol(vel_tol_in),
    tf(0.),
    dt_int(0.),
    int_eta_abs(0.),
    int_eta_rel(0.)
{
}


template<typename scalar>
planner<scalar>::~planner()
{
    delete rb_integrator_;
    delete rb_observer_;
}


template<typename scalar>
void planner<scalar>::set_simulator(const std::shared_ptr<sim_t>& sim_in)
{
    sim_ = sim_in;
    dof_ = sim_->get_model_dof();
}



template<typename scalar>
void planner<scalar>::set_interpolator(const std::shared_ptr<interp_t>& interp_in)
{
    interpolate_ = interp_in;
}

/**
 * todo: make integrator take shared ptr to sim!
 */
template<typename scalar>
void planner<scalar>::set_observers()
{
    rb_integrator_ = new rbi_t(sim_.get(), interpolate_, observer_states);
    rb_observer_ = new rbo_t(observer_states);
}


template<typename scalar>
void planner<scalar>::set_discretization(const int int_via_points_in,
                                         const scalar tf_)
{
    int_via_points = int_via_points_in;
    tf = tf_;
    dt_int = tf/(static_cast<scalar>(int_via_points-1));
}





template<typename scalar>
void planner<scalar>::set_desired_tcp_trajectory(const int np_,
                                                 const std::shared_ptr<const node_t>& node_initial_,
                                                 const std::shared_ptr<const node_t>& node_final_)
{
    std::shared_ptr<typename interp_t::result_t> interp_result;
    interp_result = interpolate_->run(node_initial_, node_final_, np_);
    tcp_traj_ = interp_result->traj_result;
}


template<typename scalar>
void planner<scalar>::set_desired_tcp_trajectory(const int np_,
                                                 std::shared_ptr<traj_t>& traj_,
                                                 const std::shared_ptr<const node_t>& node_initial_,
                                                 const std::shared_ptr<const node_t>& node_final_)
{
    interpolate_->set_guess_flag(interp_t::DIRECT_TRAJECTORY);
    std::shared_ptr<typename interp_t::result_t> interp_result;
    interp_result = interpolate_->run(traj_, node_initial_, node_final_, np_);
    tcp_traj_ = interp_result->traj_result;
}




template<typename scalar>
void planner<scalar>::run(const vector_t& state_in,
                          const vector_t& dstate_in,
                          const bool output_states)
{
}



template<typename scalar>
typename planner<scalar>::return_type_t
planner<scalar>::run_adaptive(const vector_t& state_in,
                              const vector_t& dstate_in,
                              const bool write_data_flag,
                              const char* filename)
{
    //preallocate trajectory with one via point
    std::shared_ptr<traj_t> tcp_point = std::make_shared<traj_t>(dof_,3,1);
    interpolate_->sample_one_time(tcp_point,tcp_traj_->time_initial());
    rb_integrator_->set_tcp_traj(tcp_point);
    rb_integrator_->initialize(state_in,dstate_in,lin_solver);
    matrix_t tcp_avg_e = matrix_t::Zero(7,1);
    matrix_t dtcp_avg_e = matrix_t::Zero(6,1);

    // initialize st::vectors
    typename rbi_t::state_type state_init(state_in.rows());

    // copy data (need mutable object)
    for (int idx=0; idx<state_in.rows(); ++idx) {
        state_init[idx] = state_in(idx);
    }

    return_type_t ret_val = return_type_t::SUCCESS;
    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5<typename rbi_t::state_type> error_stepper_type;
    try{        
         integrate_adaptive(make_controlled<error_stepper_type>(int_eta_abs,int_eta_rel), *rb_integrator_,
                            state_init, 0., interpolate_->duration(), dt_int, *rb_observer_);
    }
    catch(typename rbo_t::error_msg_t errmsg) {
#ifdef _VERBOSE_
        std::cerr << "\nIntegration caught exception: ";
        std::cerr << rb_observer_->out_error_msg(errmsg) << "!\n" << std::endl;
#endif //_VERBOSE_
        switch(errmsg) {
        case(rbo_t::MIN_STEPSIZE_EXCEEDED): {
            ret_val = return_type_t::INT_MIN_STEPSIZE;
            break;
        }
        case(rbo_t::MAX_INT_STEPS_EXCEEDED): {
            ret_val = return_type_t::INT_MAX_STEPS;
            break;
        }
        default:
            break;
        }
    }

    // Save the result in trajectory
    const int n_result = observer_states.m_times.size();
    Eigen::Map<matrix_t> times_(observer_states.m_times.data(), n_result, 1);

    // re-adjust times to correct values
    times_.array() += interpolate_->time_initial();

    for(int idx=0; idx<n_result; ++idx) {
        tcp_avg_e += (observer_states.m_tcps_ref[idx]-observer_states.m_tcps[idx]).cwiseAbs();
        dtcp_avg_e += (observer_states.m_dtcps_ref[idx]-observer_states.m_dtcps[idx]).cwiseAbs();
    }
    tcp_avg_e *= 1./static_cast<scalar>(n_result);
    dtcp_avg_e *= 1./static_cast<scalar>(n_result);

    if ((tcp_avg_e.block(0,0,3,1).array() > pos_tol).any())
        ret_val = return_type_t::TCP_ERROR;

    if ((tcp_avg_e.block(3,0,4,1).array() > rot_tol).any())
        ret_val = return_type_t::TCP_ERROR;

    if ((dtcp_avg_e.array() > vel_tol).any())
        ret_val = return_type_t::TCP_VEL_ERROR;

    if (ret_val == return_type_t::SUCCESS) {
        if (n_result!=tcp_traj_->via_points()) {
            tcp_traj_->resize_trajectory(times_);
        }
        for(int idx=0; idx<n_result; ++idx) {
            tcp_traj_->set_dimensions(observer_states.m_tcps[idx].transpose(),0,0,dof_,idx,1);
            tcp_traj_->set_dimensions(observer_states.m_dtcps[idx].transpose(),1,0,3,idx,1);
            tcp_traj_->set_angular_rates(observer_states.m_dtcps[idx].block(3,0,3,1).transpose(),1,idx,1);
            tcp_traj_->set_joint_dimensions(observer_states.m_states[idx].transpose(),0,0,dof_,idx,1);
            tcp_traj_->set_joint_dimensions(observer_states.m_dstates[idx].transpose(),1,0,dof_,idx,1);
            tcp_traj_->set_manipulability(observer_states.m_manips[idx],idx);
        }
    }
    if (write_data_flag)
        write_data(times_, n_result, filename);
#ifdef _VERBOSE_
    std::cerr << "avg tcp error:  " << tcp_avg_e.transpose() << " m - q" << std::endl;
    std::cerr << "avg dtcp error: " << dtcp_avg_e.transpose() << " m/s - rad/s" << std::endl;
#endif //_VERBOSE_

    observer_states.clear_all();
    return(ret_val);
}




template<typename scalar>
void planner<scalar>::write_data(const Eigen::Map<matrix_t>& times_in,
                                 const int n_result_in,
                                 const char* filename)
{
    matrix_t manip_and_time = matrix_t::Zero(n_result_in, 2);
    manip_and_time.block(0,0,manip_and_time.rows(),1) = times_in;
    for(int idx=0; idx<n_result_in; ++idx) {
        manip_and_time(idx,1) = observer_states.m_manips[idx];
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename, "manip_measure", 0).c_str(),manip_and_time);

    matrix_t tcpa_and_time = matrix_t::Zero(n_result_in, 8);
    tcpa_and_time.block(0,0,tcpa_and_time.rows(),1) = times_in;
    for(int idx=0; idx<n_result_in; ++idx) {
        tcpa_and_time.block(idx,1,1,7) = observer_states.m_tcps_ref[idx].transpose();
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename, "tcp_ref",0).c_str(),tcpa_and_time);

    matrix_t dtcpa_and_time = matrix_t::Zero(n_result_in, 7);
    dtcpa_and_time.block(0,1,1,6) = observer_states.m_dtcps_ref[0].transpose();
    dtcpa_and_time.block(0,0,dtcpa_and_time.rows(),1) = times_in;
    for(int idx=0; idx<n_result_in; ++idx) {
        dtcpa_and_time.block(idx,1,1,6) = observer_states.m_dtcps_ref[idx].transpose();
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename, "dtcp_ref",0).c_str(),dtcpa_and_time);

    matrix_t tcp_and_time = matrix_t::Zero(n_result_in, 8);
    tcp_and_time.block(0,1,1,7) = observer_states.m_tcps[0].transpose();
    tcp_and_time.block(0,0,tcp_and_time.rows(),1)= times_in;
    for(int idx=0; idx<n_result_in; ++idx) {
        tcp_and_time.block(idx,1,1,7) = observer_states.m_tcps[idx].transpose();
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename, "tcp_star",0).c_str(),tcp_and_time);

    matrix_t dtcp_and_time = matrix_t::Zero(n_result_in, 7);
    dtcp_and_time.block(0,1,1,6) = observer_states.m_dtcps[0].transpose();
    dtcp_and_time.block(0,0,dtcp_and_time.rows(),1) = times_in;
    for(int idx=0; idx<n_result_in; ++idx) {
        dtcp_and_time.block(idx,1,1,6) = observer_states.m_dtcps[idx].transpose();
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename, "dtcp_star",0).c_str(),dtcp_and_time);

    matrix_t state_and_time = matrix_t::Zero(observer_states.m_times.size(), dof_+1);
    state_and_time.block(0,0,state_and_time.rows(),1) = times_in;
    for(int idx=0; idx<n_result_in; ++idx) {
        state_and_time.block(idx,1,1,dof_) = observer_states.m_states[idx].transpose();
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename, "state_star",0).c_str(),state_and_time);

    matrix_t dstate_and_time = matrix_t::Zero(n_result_in, dof_+1);
    dstate_and_time.block(0,0,dstate_and_time.rows(),1) = times_in;
    for(int idx=0; idx<n_result_in; ++idx) {
        dstate_and_time.block(idx,1,1,dof_) = observer_states.m_dstates[idx].transpose();
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename, "dstate_star",0).c_str(),dstate_and_time);
    ++global_count;
}


// instantiate robot class
template class planner<double>;

}//namespace amp
