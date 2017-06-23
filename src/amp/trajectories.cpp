/*
FILE amp/trajectories.cpp

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

#include "amp/trajectories.hpp"


namespace amp {

template<typename scalar>
trajectory<scalar>
::trajectory(const int dofj_,
             const int ders_,
             const int viapoints_) :
    dofj(dofj_),
    ders(ders_),
    n_state(dof*ders),
    n_angle(3*(ders-1)),
    n_joint(dofj*ders),
    n_cols(2+n_state+n_angle+n_joint+n_control),
    t0(0.),
    tf(0.),
    step(0.),
    viapoints(viapoints_),
    Wp(Matrix34::Zero())
{
    trajectory_ = Matrix::Zero(viapoints,n_cols);
}



template<typename scalar>
trajectory<scalar>
::trajectory(const std::shared_ptr<const trajectory> &traj_in) :
    dofj(traj_in->dofj),
    ders(traj_in->ders),
    n_state(dof*ders),
    n_angle(3*(ders-1)),
    n_joint(dofj*ders),
    n_cols(2+n_state+n_angle+n_joint+n_control),
    t0(traj_in->time_initial()),
    tf(traj_in->time_final()),
    step(0.),
    viapoints(traj_in->via_points()),
    Wp(Matrix34::Zero())
{
    set_data(traj_in->whole_trajectory());
}



template<typename scalar>
trajectory<scalar>
::trajectory(const trajectory&  traj_in):
    dofj(traj_in.dofj),
    ders(traj_in.ders),
    n_state(dof*ders),
    n_angle(3*(ders-1)),
    n_joint(dofj*ders),
    n_cols(2+n_state+n_angle+n_joint+n_control),
    t0(traj_in.time_initial()),
    tf(traj_in.time_final()),
    step(0.),
    viapoints(traj_in.via_points()),
    Wp(Matrix34::Zero())
{
    set_data(traj_in.whole_trajectory());
}



template<typename scalar>
void trajectory<scalar>
::set_time_sequence(const scalar t0_, const scalar tf_)
{
    t0 = t0_;
    tf = tf_;
    step = (tf-t0)/static_cast<scalar>(viapoints-1);
    trajectory_.block(0,0,viapoints,1) =
            Vector::LinSpaced(viapoints,t0_,tf_);
}



template<typename scalar>
void trajectory<scalar>
::set_times(const Matrix& values)
{
    trajectory_.block(0,0,viapoints,1) = values;
}




template<typename scalar>
void trajectory<scalar>
::compute_angular_rates()
{
    Vector3 angular_temp = Vector3::Zero();
    for(int idx=0; idx<ders-1; ++idx) {
        MapMatrixConst qs(trajectory_.data()+1+(doft*ders)+idx,
                          viapoints, dofr,
                          Eigen::Stride<dynamic,dynamic>(ders,n_cols));
        MapMatrixConst qds(trajectory_.data()+1+(doft*ders)+idx+1,
                           viapoints, dofr,
                           Eigen::Stride<dynamic,dynamic>(ders,n_cols));
        for (int t=0;t<viapoints; ++t) {
            angular_temp = get_angular_rate(qs.block(t,0,1,4).transpose(),
                                            qds.block(t,0,1,4).transpose());
            trajectory_(t, 1 + n_state + idx) = angular_temp(0);
            trajectory_(t, 1 + n_state + (ders-1) + idx) = angular_temp(1);
            trajectory_(t, 1 + n_state + 2*(ders-1) + idx) = angular_temp(2);
        }
    }
}





template<typename scalar>
typename trajectory<scalar>::Vector3 trajectory<scalar>
::get_angular_rate(const Vector4& q, const Vector4& q_rate)
{
    get_Wp(q);
    return( 2. * Wp * q_rate );
}





template<typename scalar>
void trajectory<scalar>
::get_Wp(const Vector4& q)
{
    Wp(0,0) =      q(3);
    Wp(0,1) =      q(2);
    Wp(0,2) =  -1.*q(1);
    Wp(0,3) =  -1.*q(0);
    Wp(1,0) =  -1.*q(2);
    Wp(1,1) =      q(3);
    Wp(1,2) =      q(0);
    Wp(1,3) =  -1.*q(1);
    Wp(2,0) =      q(1);
    Wp(2,1) =  -1.*q(0);
    Wp(2,2) =      q(3);
    Wp(2,3) =  -1.*q(2);
}




template<typename scalar>
void trajectory<scalar>
::normalize_quaternions()
{
    Quaternion quater = Quaternion::Identity();
    MapArray quater_array(trajectory_.data()+1+doft*ders,
                          4, 1,
                          Eigen::Stride<dynamic,dynamic>(1,ders));
    if((quater_array.abs() < 1e-6).all()) {
        quater_array = quater.coeffs();
    }
    else {
        quater.coeffs() = quater_array;
        quater.normalize();
        quater_array = quater.coeffs();
    }
    for(int idx=1; idx<viapoints; ++idx) {
        new (&quater_array) MapArray(trajectory_.data()+idx*n_cols+1+doft*ders,
                                     4, 1,
                                     Eigen::Stride<dynamic,dynamic>(1,ders));
        if((quater_array.abs() < 1e-6).all())
            quater_array = quater.coeffs();
        else {
            quater.coeffs() = quater_array;
            quater.normalize();
            quater_array = quater.coeffs();
        }
    }
}




template<typename scalar>
void trajectory<scalar>
::check_antipodal_quaternion(const Array4& array_q0)
{
    MapArray array_q1(trajectory_.data()+1+doft*ders,
                      4, 1,
                      Eigen::Stride<dynamic,dynamic>(1,ders));

    if ( ((array_q1+array_q0).abs()<0.1).all() )  {
        array_q1 *= -1.0;
    }
}




template<typename scalar>
void trajectory<scalar>
::check_antipodal_quaternions()
{
    MapArray array_q0(trajectory_.data()+1+doft*ders,
                      4, 1,
                      Eigen::Stride<dynamic,dynamic>(1,ders));
    MapArray array_q1(trajectory_.data()+n_cols+1+doft*ders,
                      4, 1,
                      Eigen::Stride<dynamic,dynamic>(1,ders));

    for (int idx=1; idx<viapoints-1; ++idx) {
        // check the previous pair and flip pole of idx if antipodel
        if ( (array_q1-array_q0).abs().maxCoeff()
             > (-1.0*array_q1-array_q0).abs().maxCoeff() ) {
            array_q1 *= -1.0;
        }
        // rest mapped arrays to next pair.
        new (&array_q0) MapArray(trajectory_.data()+idx*n_cols+1+doft*ders,
                                 4, 1,
                                 Eigen::Stride<dynamic,dynamic>(1,ders));
        new (&array_q1) MapArray(trajectory_.data()+(idx+1)*n_cols+1+doft*ders,
                                 4, 1,
                                 Eigen::Stride<dynamic,dynamic>(1,ders));
    }
}





template<typename scalar>
void trajectory<scalar>
::set_derivatives(const Matrix& values,
                  const int dim_,
                  const int der_,
                  const int nder_,
                  const int via_,
                  int total_)
{
    trajectory_.block(via_, 1+dim_*ders+der_,
                      (total_==0 ? viapoints : total_), nder_) = values;
}






template<typename scalar>
void trajectory<scalar>
::set_dimensions(const Matrix& values,
                 const int der_,
                 const int dim_,
                 const int ndof_,
                 const int via_,
                 int total_)
{
    for(int idx=0; idx<ndof_; ++idx) {
        trajectory_.block(via_,1+(dim_+idx)*ders+der_,
                          (total_==0 ? viapoints : total_),1) =
                values.block(0, idx, (total_==0 ? viapoints : total_), 1);
    }
}





template<typename scalar>
void trajectory<scalar>
::set_angular_rates(const Matrix& values,
                    const int der_,
                    const int via_,
                    int total_)
{
    for(int idx=0; idx<3; ++idx) {
        trajectory_.block(via_, 1+n_state+idx*(ders-1)+(der_-1),
                          (total_==0 ? viapoints : total_),1) =
                values.block(0, idx, (total_==0 ? viapoints : total_), 1);
    }
}




template<typename scalar>
void trajectory<scalar>
::set_controls(const Matrix& values,
               const int dim_,
               const int ndof_,
               const int via_,
               int total_)
{
    trajectory_.block(via_, 1+n_state+n_angle+dim_,
                      (total_==0 ? viapoints : total_), ndof_) =
            values.block(0, 0, (total_==0 ? viapoints : total_), ndof_);
}





template<typename scalar>
void trajectory<scalar>
::set_joint_derivatives(const Matrix& values,
                        const int dim_,
                        const int der_,
                        const int nder_,
                        const int via_,
                        int total_)
{
    trajectory_.block(via_, 1+n_state+n_angle+n_control+dim_*ders+der_,
                      (total_==0 ? viapoints : total_), nder_) = values;
}





template<typename scalar>
void trajectory<scalar>
::set_joint_dimensions(const Matrix& values,
                       const int der_,
                       const int dim_,
                       const int ndof_,
                       const int via_,
                       int total_)
{
    for(int idx=0; idx<ndof_; ++idx) {
        trajectory_.block(via_, 1+n_state+n_angle+n_control+(dim_+idx)*ders+der_,
                          (total_==0 ? viapoints : total_), 1) =
                values.block(0, idx, (total_==0 ? viapoints : total_), 1);
    }
}





template<typename scalar>
void trajectory<scalar>
::set_manipulability(const Matrix& values,
                     const int via_,
                     int total_)
{
    trajectory_.block(via_, 1+n_state+n_angle+n_control+n_joint,
                      (total_==0 ? viapoints : total_), 1) =
            values.block(0, 0, (total_==0 ? viapoints : total_), 1);
}




template<typename scalar>
void trajectory<scalar>
::set_manipulability(const scalar value, const int via_)
{
    trajectory_(via_,1+n_state+n_angle+n_control+n_joint) = value;
}




template<typename scalar>
void trajectory<scalar>
::node(std::shared_ptr<node_t> node_i, const int via_) const
{
    node_i->time_curr = times(via_,1)(0);

    for (int idx=0; idx<ders; ++idx) {
        node_i->state.block(0,idx,dof,1) = dimensions(idx,via_,1).transpose();
    }

    for (int idx=0; idx<ders-1; ++idx) {
        node_i->angular.block(0,idx,3,1) = angular_rates(idx,via_,1).transpose();
    }

    for (int idx=0; idx<ders-1; ++idx) {
        node_i->joints.block(0,idx,dofj,1) = joint_dimensions(idx,via_,1).transpose();
    }

    node_i->control = controls(via_,1).transpose();
}






template<typename scalar>
void trajectory<scalar>
::set_node(std::shared_ptr<const node_t> node_i, const int via_)
{
    trajectory_(via_,0) = node_i->time_curr;

    for (int idx=0; idx<ders; ++idx) {
        set_dimensions(node_i->state.block(0,idx,dof,1).transpose(),idx,0,dof,via_,1);
    }

    for (int idx=0; idx<ders-1; ++idx) {
        set_angular_rates(node_i->angular.block(0,idx,3,1).transpose(),idx,via_,1);
    }

    for (int idx=0; idx<ders-1; ++idx) {
        set_joint_dimensions(node_i->joints.block(0,idx,dofj,1).transpose(),idx,0,dofj,via_,1);
    }

    set_controls(node_i->control.transpose(),0,n_control,via_,1);
}




// instantiate trajectory
template class trajectory<float>;
template class trajectory<double>;

} // namespace amp

