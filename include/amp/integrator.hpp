/*
FILE amp/integrator.hpp

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

#ifndef _AMP_INTEGRATOR_HPP
#define _AMP_INTEGRATOR_HPP

#include <iostream>
#include <memory>

#ifdef USE_EIGEN_PKGCONFIG_HEADERS
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#else
#include <Eigen/Dense>
#include <Eigen/StdVector>
#endif

#include "amp/simulation.hpp"
#include "amp/trajectories.hpp"
#include "amp/interpolate.hpp"

#ifdef DEBUG
#define _DEBUG_
#endif

namespace amp {



/**
 *
 */
template <typename scalar>
struct observers
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,1> vector_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> matrix_t;
    typedef std::vector<scalar> manip_observer_t;
    typedef std::vector<int> rank_observer_t;
    typedef std::vector<matrix_t, Eigen::aligned_allocator<matrix_t> > error_observer_t;
    typedef std::vector<vector_t, Eigen::aligned_allocator<vector_t> > state_observer_t;

    scalar tcp_error_avg;
    scalar dtcp_error_avg;
    matrix_t tcp_i_;
    matrix_t tcp_star_;
    matrix_t dtcp_i_;
    matrix_t dtcp_star_;
    matrix_t dtcp_c_;
    vector_t state_star_;
    vector_t dstate_star_;
    scalar manipulability_;

    /**
     * @brief m_times
     */
    std::vector<scalar> m_times;


    /**
     * @brief m_states
     */
    state_observer_t m_states, m_dstates;

    /**
     * @brief m_tcps
     */
    error_observer_t m_tcps, m_dtcps, m_tcps_ref, m_dtcps_ref;

    /**
     * @brief m_manips
     */
    manip_observer_t m_manips;

    observers() :
        tcp_error_avg(0.),
        dtcp_error_avg(0.) { }

    void clear_all () {
        tcp_error_avg =0.;
        dtcp_error_avg =0.;
        m_times.clear();
        m_states.clear();
        m_dstates.clear();
        m_tcps.clear();
        m_dtcps.clear();
        m_tcps_ref.clear();
        m_dtcps_ref.clear();
        m_manips.clear();
    }
};



template<typename scalar>
struct rb_observer
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,1> vector_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> matrix_t;
    typedef observers<scalar> observers_t;
    typedef std::vector<scalar> state_type;
    typedef Eigen::Map<const vector_t> MapStateConst;

    enum error_msg_t {
        MAX_INT_STEPS_EXCEEDED = 0,
        MIN_STEPSIZE_EXCEEDED
    };

    const char* out_error_msg(const error_msg_t msg_) const{
        const char* out_char;
        switch(msg_)
        {
        case(MAX_INT_STEPS_EXCEEDED):
        {
            out_char = "MAX_INT_STEPS_EXCEEDED";
            break;
        }
        case(MIN_STEPSIZE_EXCEEDED):
        {
            out_char = "MIN_STEPSIZE_EXCEEDED";
            break;
        }
        default:
            out_char = "UNDEFINED";
        }
        return(out_char);
    }

    observers_t& observers_;
    const int max_int_steps_allowed;
    const scalar min_int_stepsize_allowed;

    rb_observer(observers_t& observers_in) :
        observers_(observers_in),
        max_int_steps_allowed(1000),
        min_int_stepsize_allowed(1e-6)
    { }

    void operator()( const state_type &x , double t )
    {

        if(observers_.m_times.size()>1 &&
                (t-observers_.m_times.back())<min_int_stepsize_allowed) {
            throw(MIN_STEPSIZE_EXCEEDED);
        }

        MapStateConst x_map(x.data(), x.size(), 1);
        observers_.m_times.push_back( t );
        observers_.m_states.push_back(x_map);
        observers_.m_dstates.push_back(observers_.dstate_star_);
        observers_.m_tcps.push_back(observers_.tcp_star_);
        observers_.m_dtcps.push_back(observers_.dtcp_star_);
        observers_.m_tcps_ref.push_back(observers_.tcp_i_);
        observers_.m_dtcps_ref.push_back(observers_.dtcp_i_);
        observers_.m_manips.push_back(observers_.manipulability_);

        // if(observers_.m_states.size() > max_int_steps_allowed) {
        //     throw(MAX_INT_STEPS_EXCEEDED);
        // }
    }
};



/**
 * @class rb_integrator
 * Function wrapper for integration using boost/odeint library.
 */
template<typename scalar>
class rb_integrator
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef simulation<scalar> sim_t;
    typedef trajectory<scalar> traj_t;
    typedef interpolator<scalar> interp_t;
    typedef typename sim_t::vector_t vector_t;
    typedef typename sim_t::matrix_t matrix_t;
    typedef std::vector<scalar> state_type;
    typedef Eigen::Map<vector_t> MapState;
    typedef Eigen::Map<const vector_t> MapStateConst;
    typedef std::vector<scalar> manip_observer_t;
    typedef std::vector<matrix_t, Eigen::aligned_allocator<matrix_t> > error_observer_t;
    typedef std::vector<vector_t, Eigen::aligned_allocator<vector_t> > state_observer_t;

    typedef typename sim_t::lin_alg_t lin_alg_t;

    sim_t* sim_;
    const std::shared_ptr<interp_t> interpolate_;
    std::shared_ptr<traj_t> tcp_traj_;

    observers<scalar>& observers_;

    lin_alg_t lin_solver_;
    matrix_t tcp_0;
    matrix_t& tcp_i;
    matrix_t& tcp_star;
    matrix_t& dtcp_i;
    matrix_t& dtcp_star;
    matrix_t& dtcp_c;
    vector_t& state_star;
    vector_t& dstate_star;
    scalar& manipulability;

    rb_integrator(sim_t* sim_in,
                  const std::shared_ptr<interp_t>& interpolate_in,
                  observers<scalar>& observers_in,
                  lin_alg_t lin_solver_in = lin_alg_t::LU,
                  const int fixed_joint_in=-1);

    ~rb_integrator(){ }

    void initialize(const vector_t& state, const vector_t& dstate, const lin_alg_t lin_solver_in);

    void set_tcp_traj(std::shared_ptr<traj_t> tcp_);

    void operator() (const state_type& y , state_type& dydt , const scalar t);
};
}
#endif // INTEGRATOR_HPP
