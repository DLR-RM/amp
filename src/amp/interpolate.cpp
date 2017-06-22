/*
FILE amp/interpolate.cpp

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

#include "amp/interpolate.hpp"


namespace amp {




template<typename scalar>
interpolator<scalar>
::interpolator(const int via_points_, const int dofj_, const int ders_, const guess_flag_t guess_type_in) :
    edge_flag(DEFAULT_EDGE),
    guess_flag(guess_type_in),
    n_fbcs(0),
    dofj(dofj_),
    ders(ders_),
    n_p(0),
    n_pd(node_t::dof),
    via_points(via_points_),
    duration_(0.),
    time_step(0.),
    time_final_(0.),
    time_initial_(0.)
{ }





template<typename scalar>
interpolator<scalar>
::~interpolator ()
{
}





template<typename scalar>
void interpolator<scalar>
::set_via_points(const int via_points_in)
{
    via_points = via_points_in;
}





template<typename scalar>
void interpolator<scalar>
::set_num_params(const int n_p_in) 
{
	n_p = n_p_in;
}





template<typename scalar>
std::shared_ptr<typename interpolator<scalar>::result_t>
interpolator<scalar>
::run(const std::shared_ptr<const node_t>& state_from_in,
      const std::shared_ptr<const node_t>& state_to_in,
      const int n_p_)
{
    n_p = n_p_;
    state_from = state_from_in;
    state_to = state_to_in;

    time_initial_ = state_from->time_curr;
    time_final_ = state_to->time_curr;
    duration_ = time_final_ - time_initial_;
    time_step = duration_ / static_cast<scalar>(via_points-1);

    // set initial guess
    splines.clear();
    Matrix parameters = Matrix::Zero(n_p, n_pd);

    // why is that necessary?

    std::shared_ptr<trajectory_t> traj_guess =
            std::make_shared<trajectory_t>(dofj,ders,via_points);
    traj_guess->set_time_sequence(time_initial_,time_final_);

    set_spline_parameterization(parameters, traj_guess);
    traj_guess->normalize_quaternions();
    traj_guess->check_antipodal_quaternions();
    traj_guess->compute_angular_rates();

    // set result
    std::shared_ptr<result_t> result_ = std::make_shared<result_t>(3,dofj,ders,n_p,n_pd);
    result_->traj_result = traj_guess;
    result_->node_start = state_to;
    result_->node_end = state_from;
    result_->parameters = parameters;

    return(result_);
}


template<typename scalar>
std::shared_ptr<typename interpolator<scalar>::result_t>
interpolator<scalar>
::run(std::shared_ptr<trajectory_t>& traj_guess,
      const std::shared_ptr<const node_t>& state_from_in,
      const std::shared_ptr<const node_t>& state_to_in,
      const int n_p_)
{
    n_p = n_p_;
    state_from = state_from_in;
    state_to = state_to_in;

    time_initial_ = state_from->time_curr;
    time_final_ = state_to->time_curr;
    duration_ = time_final_ - time_initial_;
    time_step = duration_ / static_cast<scalar>(via_points-1);

    std::cerr << "state from: t=" << state_from->time_curr << std::endl;
    std::cerr << state_from->state.transpose() << std::endl;
    std::cerr << "state to: t=" << state_to->time_curr << std::endl;
    std::cerr << state_to->state.transpose() << std::endl;

    // set initial guess
    splines.clear();
    Matrix parameters = Matrix::Zero(n_p, n_pd);

    // why is that necessary?
    set_spline_parameterization(parameters, traj_guess);
    traj_guess->normalize_quaternions();
    traj_guess->check_antipodal_quaternions();
    traj_guess->compute_angular_rates();

    // set result
    std::shared_ptr<result_t> result_ =
            std::make_shared<result_t>(3,dofj,ders,n_p,n_pd);
    result_->traj_result = traj_guess;
    result_->node_start = state_to;
    result_->node_end = state_from;
    result_->parameters = parameters;

    return(result_);
}


template<typename scalar>
void interpolator<scalar>
::set_spline_parameterization(Matrix& parameters,
                              std::shared_ptr<trajectory_t> &trajectory_)
{
    typename sgt::knot_vector knots;
    set_knots(knots);

    if (guess_flag==DIRECT_TRAJECTORY){
        set_guess_from_trajectory(parameters, trajectory_, knots);
    }
    else if (guess_flag==DIRECT_PARAMETERS) {
        set_guess_from_parameters(parameters, trajectory_, knots);
    }
    else {
        set_guess_with_qp(parameters, trajectory_, knots);
    }
}






template<typename scalar>
void interpolator<scalar>
::set_knots(typename sgt::knot_vector& knots)
{

    const int n_kn = n_p+(2*n_fbcs)+ order;
    knots = sgt::knot_vector::Zero(n_kn);
    knots.segment(n_kn-order,order) = rowVector::Constant(order,1.);

    knots.segment(order-1,n_kn-2*(order-1))
            = rowVector::LinSpaced(n_kn-2*(order-1),0.,1.);
}






template<typename scalar>
std::shared_ptr<typename interpolator<scalar>::trajectory_t> interpolator<scalar>
::resample(Matrix& parameters_, const int via_in)
{
    std::shared_ptr<trajectory_t> traj_ =
            std::make_shared<trajectory_t>(dofj,ders,via_in);

    traj_->set_time_sequence(time_initial_,time_final_);

    // make spline objects for each dimension
    int dim=0;
    for (auto spline=splines.begin(); spline!=splines.end(); ++spline, ++dim) {
        std::shared_ptr<sgt> myspline = *spline;
        myspline->set_via_points(via_in);
        myspline->compute_bases(order);

        // set new values
        traj_->set_derivatives(myspline->get_spline_values(order-1).block(0,1,via_in,order-1),
                                     dim,0,order-1);
    }
    traj_->normalize_quaternions();
    traj_->check_antipodal_quaternions();
    traj_->compute_angular_rates();

    return(traj_);
}





template<typename scalar>
void interpolator<scalar>
::sample_one_time(std::shared_ptr<trajectory_t>& traj_, const scalar via_in)
{
    int dim=0;
    for (auto spline=splines.begin(); spline!=splines.end(); ++spline, ++dim) {
        std::shared_ptr<sgt> spl = *spline;
        // adjust time sequence to 0
        typename sgt::curve point = spl->get_one_spline_value(via_in-time_initial_,order-1);
        traj_->set_derivatives(point.block(0,1,1,order-1),dim,0,order-1,0,1);
    }
}






template<typename scalar>
void interpolator<scalar>
::set_guess_with_qp(Matrix& parameters_,
                    std::shared_ptr<trajectory_t>& trajectory_,
                    const typename sgt::knot_vector& knots)
{
    // make spline objects for each dimension
    for (int dim=0; dim<n_pd; ++dim) {
        std::shared_ptr<sgt> myspline = std::make_shared<sgt>(n_fbcs, n_p+2*n_fbcs, knots, duration_);
        splines.push_back(myspline);
        myspline->set_via_points(via_points);
        myspline->compute_bases(order);
        myspline->init_generator(state_from->state.block(dim,0,1,n_fbcs).transpose(),
                                 state_to->state.block(dim,0,1,n_fbcs).transpose(),
                                 sgt::p_vector::Zero(n_p));

        // solve min jerk QP and set resultant parameter values as initial guess
        myspline->qp_initial_guess(2,0,Eigen::Matrix<scalar,0,2>::Zero());
        parameters_.block(0,dim,n_p,1) = myspline->parameters();
        trajectory_->set_derivatives(myspline->get_spline_values(order-1).block(0,1,via_points,order-1),
                                     dim,0,order-1);
    }
}






template<typename scalar>
void interpolator<scalar>
::set_guess_from_parameters(const Matrix& parameters_,
                            std::shared_ptr<trajectory_t>& trajectory_,
                            const typename sgt::knot_vector& knots)
{
    for (int dim=0; dim<n_pd; ++dim) {
        std::shared_ptr<sgt> myspline = std::make_shared<sgt>(n_fbcs, n_p+2*n_fbcs, knots, duration_);
        splines.push_back(myspline);
        myspline->set_via_points(via_points);
        myspline->compute_bases(order);

        myspline->init_generator(state_from->state.block(dim,0,1,n_fbcs).transpose(),
                    state_to->state.block(dim,0,1,n_fbcs).transpose(),
                    sgt::p_vector::Zero(n_p) );

        myspline->update_vertices(parameters_.block(0,dim,n_p,1));
    }
}






template<typename scalar>
void interpolator<scalar>
::set_guess_from_trajectory(Matrix& parameters_,
                            std::shared_ptr<trajectory_t>& trajectory_,
                            const typename sgt::knot_vector& knots)
{
    // adjusted time sequence to 0
    if (trajectory_->time_initial()>1e-6) {
        trajectory_->set_times(trajectory_->times().array()-time_initial_);
    }

    for (int dim=0; dim<n_pd; ++dim) {
        std::shared_ptr<sgt> myspline = std::make_shared<sgt>(n_fbcs, n_p+2*n_fbcs, knots, duration_);
        splines.push_back(myspline);
        myspline->set_via_points(trajectory_->via_points());
        myspline->compute_bases_nonuniform(trajectory_->times(), order);
        myspline->init_generator( state_from->state.block(dim,0,1,n_fbcs).transpose(),
                    state_to->state.block(dim,0,1,n_fbcs).transpose(),
                    sgt::p_vector::Zero(n_p));

        // fit spline to position curve of current DoF
        myspline->llsq_initial_guess((trajectory_->derivatives(dim)).block(0,0,trajectory_->via_points(),1), 0);
        parameters_.block(0,dim,n_p,1) = myspline->parameters();

    }

    int dim=0;
    trajectory_->resize_trajectory(Vector::LinSpaced(via_points,time_initial_,time_final_));
    for (auto spline=splines.begin(); spline!=splines.end(); ++spline, ++dim) {
        (*spline)->set_via_points(via_points);
        (*spline)->compute_bases(order);
        trajectory_->set_derivatives((*spline)->get_spline_values(order-1).block(0,1,via_points,order-1),
                                     dim,0,order-1);
    }
}



// instantiate robot class
template class interpolator<float>;
template class interpolator<double>;


} // end namespace amp
