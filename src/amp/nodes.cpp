/*
FILE amp/nodes.cpp

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

#include "amp/nodes.hpp"

namespace amp {


template<typename scalar>
node<scalar>
::node(const int dofj_, const int ders_) :
    doa(6),
    dofj(dofj_),
    ders(ders_),
    time_curr(0.),
    state(state_t::Zero(dof,ders)),
    control(control_t::Zero(doa,1)),
    angular(angular_t::Zero(3,ders)),
    joints(joint_state_t::Zero(dofj,ders))
{
    state.block(doft,0,dofr,1) = Quaternion(Quaternion::Identity()).coeffs();
}

template<typename scalar>
node<scalar>
::node(const std::shared_ptr<const node<scalar> > &m_) :
    doa(6),
    dofj(m_->dofj),
    ders(m_->ders),
    time_curr(m_->time_curr),
    state(m_->state),
    control(m_->control),
    angular(m_->angular),
    joints(m_->joints)
{ }


template<typename scalar>
void node<scalar>
::copy_data(const std::shared_ptr<const node<scalar> > &m_)
{
    time_curr = m_->time_curr;
    state = m_->state;
    control = m_->control;
    angular = m_->angular;
    joints = m_->joints;
}


// instantiate node
template class node<float>;
template class node<double>;

} // end namespace amp
