/*
FILE amp/simulation.cpp

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

#include "amp/simulation_rbdl.hpp"


namespace amp {


template<typename scalar>
rbdl_interface<scalar>::rbdl_interface(const std::string& filename) :
    model(new model_t)
{
    if (!RigidBodyDynamics::Addons::
            URDFReadFromFile(filename.c_str(), model, false)) {
        std::cerr << "Error loading model: " << filename.c_str() << std::endl;
        abort();
    }
    model->gravity = vector3_t(0., 0., -9.81);
    this->tcp_id_ = model->dof_count;
}


template<typename scalar> 
rbdl_interface<scalar>::~rbdl_interface() 
{
    delete model;
}


template<typename scalar> 
const int rbdl_interface<scalar>::get_model_dof() const {
    return(model->dof_count);
}


template<typename scalar> 
void rbdl_interface<scalar>::fwd_kin(const vector_t& state_, 
                                     const vector_t& dstate_)
{
    vector_t ddstate_ = vector_t::Zero(model->dof_count);
    RigidBodyDynamics::UpdateKinematics(*model, state_, dstate_, ddstate_);
}


template<typename scalar> 
void rbdl_interface<scalar>::fwd_dyn(vector_t& ddstate_, 
                                     const vector_t& state_, const vector_t& dstate_, const vector_t& tau_)
{
}


template<typename scalar> 
typename rbdl_interface<scalar>::matrix_t
rbdl_interface<scalar>::tcp_pose(const vector_t& state_)
{
    matrix_t tcp = matrix_t::Zero(7,1);
    matrix3_t rot =
            RigidBodyDynamics::CalcBodyWorldOrientation(*model,state_,this->tcp_id_,false);
    Eigen::Quaternion<scalar> q_(rot);
    tcp.block(0,0,3,1) =
            RigidBodyDynamics::CalcBodyToBaseCoordinates(*model,state_,this->tcp_id_,
                                                         this->tcp_offset_.head(3),false);
    tcp.block(3,0,4,1) = q_.inverse().coeffs();
    return(tcp);
}



template<typename scalar> 
typename rbdl_interface<scalar>::matrix3_t
rbdl_interface<scalar>::rotation_to_tcp_frame(const vector_t& state_)
{
    return(RigidBodyDynamics::
           CalcBodyWorldOrientation(*model,state_,this->tcp_id_,false));
}


template<typename scalar> 
typename rbdl_interface<scalar>::matrix3_t
rbdl_interface<scalar>::rotation_to_base_frame(const vector_t& state_)
{
    return(RigidBodyDynamics::
           CalcBodyWorldOrientation(*model,state_,this->tcp_id_,false).transpose());
}


template<typename scalar>
scalar rbdl_interface<scalar>::manipulability(const jacobian_t& J_) const
{
    scalar det = 0.;
    scalar manipbty = 0.;
    det =(J_*J_.transpose()).determinant();
    if (det >= 1e-6)  {
        manipbty = sqrt(det);
    }
    return(manipbty);
}


template<typename scalar>
typename rbdl_interface<scalar>::jacobian_t
rbdl_interface<scalar>::get_J(const vector_t& state_) const
{
    jacobian_t J_ = jacobian_t::Zero(6,model->dof_count);
    RigidBodyDynamics::CalcPointJacobian6D(*model, state_, this->tcp_id_,
                                           this->tcp_offset_.head(3), J_, false);
    return (J_);
}


template<typename scalar> 
typename rbdl_interface<scalar>::vector_t
rbdl_interface<scalar>::rbdl_interface::inv_kin(const matrix_t& dtcp_,
                                                const vector_t& state_, const typename sim_t::lin_alg_t lin_)
{
    jacobian_t J_i = get_J(state_);
    matrix_t dtcp_rbdl_= matrix_t::Zero(6,1);
    dtcp_rbdl_.block(0,0,3,1) = dtcp_.block(3,0,3,1);
    dtcp_rbdl_.block(3,0,3,1) = dtcp_.block(0,0,3,1);
    vector_t dstate_= vector_t::Zero(model->dof_count);
    switch(lin_)
    {
    case(this->SVD): {
        svd_t svd_decomp(J_i,Eigen::ComputeThinU|Eigen::ComputeThinV);
        dstate_= svd_decomp.solve(dtcp_rbdl_);
        break;
    }
    case(this->QR): {
        qr_t qr_decomp(J_i);
        dstate_= qr_decomp.solve(dtcp_rbdl_);
        break;
    }
    case(this->LDLT): {
        ldlt_t ldlt_decomp(J_i);
        dstate_= ldlt_decomp.solve(dtcp_rbdl_);
        break;
    }
    case(this->LLT): {
        llt_t llt_decomp(J_i);
        dstate_= llt_decomp.solve(dtcp_rbdl_);
        break;
    }
    case(this->LU):
    default: {
        lu_t lu_decomp(J_i);
        lu_decomp.setThreshold(1e-5);
        dstate_= lu_decomp.solve(dtcp_rbdl_);
    }
    }
    return(dstate_);
}


template<typename scalar> 
typename rbdl_interface<scalar>::matrix_t
rbdl_interface<scalar>::tcp_error(const matrix_t& tcp_,
                                  const vector_t& state_)
{
    return(tcp_-tcp_pose(state_));
}


template<typename scalar> 
typename rbdl_interface<scalar>::matrix_t
rbdl_interface<scalar>::tcp_error_correction(const matrix_t& tcp_e)
{
    return(tcp_e.block(0,0,6,1).cwiseProduct(this->error_gains_));
}


template<typename scalar> 
typename rbdl_interface<scalar>::vector_t
rbdl_interface<scalar>::diff_kin(const vector_t& state_,
                                 const vector_t& dstate_, scalar dt_)
{
    return(state_+(dt_*dstate_));
}


template<typename scalar> 
typename rbdl_interface<scalar>::matrix_t rbdl_interface<scalar>
::inv_kin_error_check(const vector_t& state_, const vector_t& dstate_)
{
    matrix_t dtcp_ = matrix_t::Zero(6,1);
    matrix_t dtcp_rbdl_ = matrix_t::Zero(6,1);
    dtcp_rbdl_ = get_J(state_)*dstate_;
    dtcp_.block(0,0,3,1) = dtcp_rbdl_.block(3,0,3,1);
    dtcp_.block(3,0,3,1) = dtcp_rbdl_.block(0,0,3,1);
    return(dtcp_);
}


template class rbdl_interface<double>;
/* template class rbdl_interface<float>; */

} // namespace amp
