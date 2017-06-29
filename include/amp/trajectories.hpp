/*
FILE amp/trajectories.hpp

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

#ifndef TRAJECTORIES_HPP
#define TRAJECTORIES_HPP

#ifdef USE_EIGEN_PKGCONFIG_HEADERS
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif

#include <memory>
#include "amp/nodes.hpp"


#ifdef VERBOSE
#define _VERBOSE_OUTPUT_
#endif //VERBOSE

#ifdef DEBUG
#define _DEBUG_
#endif //DEBUG

namespace amp {


/**
 * \brief Class which holds trajectory data as an eigen matrix and provides
 * selector methods to read and write logical blocks of the matrix.
 *
 * The Trajectory and node objects are dependent on eachother such that the
 * trajectory object stores the 'persistent' trajectory information and state
 * objects are used as local containers, to compute or be computed from a
 * trajectory instance.
 *
 * NOTE: selectors --> picking sliced/ reshaped data from trajectory class
 * using eigen Map objects
 */
template<typename scalar>
class trajectory {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    static const int doft = 3;
    static const int dofr = 4;
    static const int n_control = 6;
    static const int dof = doft+dofr;
    const int dofj;
    const int ders;
    const int n_state;   // = dof*ders;
    const int n_angle;   // = 3*(ders-1);
    const int n_joint;   // = dofj*ders
    const int n_cols;    // = 1+n_state+n_angle+n_control+n_joint;
    static const int dynamic = Eigen::Dynamic;

    enum time_derivatives {
        position = 0,
        velocity,
        acceleration,
        jerk,
        snap,
        crackle,
        pop
    };

    typedef amp::node<scalar> node_t;
    typedef Eigen::Matrix<scalar,3,4> Matrix34;
    typedef Eigen::Matrix<scalar,4,1> Vector4;
    typedef Eigen::Matrix<scalar,3,1> Vector3;
    typedef Eigen::Array<scalar,4,1> Array4;
    typedef Eigen::Quaternion<scalar> Quaternion;
    typedef Eigen::Transform<scalar,3,Eigen::AffineCompact> Transform;
    typedef Eigen::Matrix<scalar,dynamic,1> Vector;
    typedef Eigen::Matrix<scalar,dynamic,dynamic> Matrix;
    typedef Eigen::Matrix<scalar,dynamic,dynamic,Eigen::RowMajor> MatrixRowMajor;
    typedef Eigen::Map<Array4,0,Eigen::Stride<dynamic,dynamic> > MapArray;
    typedef Eigen::Map<const Matrix,0,Eigen::Stride<dynamic,dynamic> > MapMatrixConst;


    trajectory(const int dofj_,
               const int ders_,
               const int viapoints_);

    trajectory(const std::shared_ptr<const trajectory>&  traj_in);

    trajectory(const trajectory&  traj_in);

    inline scalar const time_initial() const {
        return(t0);
    }

    inline scalar const time_final() const {
        return(tf);
    }

    inline scalar const time_step() const {
        return(step);
    }

    int const via_points() const {
        return(viapoints);
    }

    int const n_columns() const {
        return(n_cols);
    }

    /**
     * @brief set_time_sequence sets a linear sequence
     * @param t0_
     * @param tf_
     */
    void set_time_sequence(const scalar t0_, const scalar tf_);

    /**
     * @brief set_times sets sequence defined in \a values
     * @param values
     */
    void set_times(const Matrix& values);

    /**
     * @brief sets all data for 1 or multiple dof with dim_ (starting) and
     * ndof_ (total)
     */
    void set_zero() {
        trajectory_.setZero(trajectory_.rows(), trajectory_.cols());
    }

    /**
     * @brief resize_trajectory resizes trajectory to new number of via points.
     *  This deletes all data currently stored in the matrix. Class member data
     *  is otherwise unaffected.
     * @param times_new
     */
    inline void resize_trajectory(const Vector& times_new) {
        viapoints = times_new.rows();
        trajectory_.setZero(trajectory_.rows(), trajectory_.cols());
        trajectory_.resize(viapoints,Eigen::NoChange);
        trajectory_.block(0, 0, viapoints, 1) = times_new;
    }


    /**
     * @brief sets all data for 1 or multiple dof with dim_ (starting) and
     * ndof_ (total)
     */
    void set_data(const Matrix& values,
                  const int via_=0,
                  int total_=0) {
        trajectory_.block(via_, 0, (total_==0 ? viapoints : total_), n_cols) = values;
    }


    /**
     * @brief this method is for easy debug output of trajectories
     */
    Matrix const whole_trajectory() const {
        return(trajectory_);
    }


    /**
     * @brief
     */
    Matrix const times(const int via_=0,
                       int total_=0) const {
        return(trajectory_.block(via_, 0, (total_==0 ? viapoints : total_), 1));
    }


    /**
     * @brief returns reference to row vector or matrix of row vectors
     * containing all dimensions for given time derivative, starting at via_ to
     * size of mapped is (total_via x dof) defaults to entire time sequence of
     * position.
     *  e.g. |x(t0)  y(t0)  z(t0)  qx(t0)  qy(t0)  qz(t0)  qw(t0)|
     *       | ...    ...    ...    ...     ...     ...     ...  |
     *       |x(tf)  y(tf)  z(tf)  qx(tf)  qy(tf)  qz(tf)  qw(tf)|
     */
    MapMatrixConst const dimensions(const int der_=position,
                                    const int via_=0,
                                    int total_=0) const {
        return( MapMatrixConst(trajectory_.data()+(via_*n_cols)+1+der_,
                               (total_==0 ? viapoints : total_), dof,
                               Eigen::Stride<dynamic,dynamic>(ders,n_cols)) );
    }


    /**
     * @brief returns reference to row vector or matrix of row vectors
     * containing all derivatives for given dimension, starting at via_ to size
     * of mapped
     *  e.g. |x(t0)  xd(t0) ...  xdddd(t0)|
     *       | ...    ...   ...     ...   |
     *       |x(tf)  xd(tf) ...  xdddd(tf)|
     */
    MapMatrixConst const derivatives(const int dim_=0,
                                     const int via_=0,
                                     int total_=0) const {
        return( MapMatrixConst(trajectory_.data()+(via_*n_cols)+1+(dim_*ders),
                               (total_==0 ? viapoints : total_), ders,
                               Eigen::Stride<dynamic,dynamic>(1,n_cols)) );
    }


    /**
     * @brief returns reference to column vector or matrix of column vectors
     * containing all angular rates of given order (e.g. 1==velocity), starting
     * at via_ for total_ time steps, mapped size (total_via x 3)
     *  e.g. |wx(t0) wy(t0) wz(t0)|
     *       | ...    ...    ...  |
     *       |wx(tf) wy(tf) wz(tf)|
     */
    MapMatrixConst const angular_rates(const int der_=velocity,
                                       const int via_=0,
                                       int total_=0) const {
        return(MapMatrixConst(trajectory_.data()+(via_*n_cols)+1+n_state+(der_-1),
                              (total_==0 ? viapoints : total_), 3,
                              Eigen::Stride<dynamic,dynamic>((ders-1),n_cols)) );
    }


    /**
     * @brief returns reference to column vector or matrix of column vectors
     * containing all angular rates of given order (e.g. 1==velocity), starting
     * at via_ for total_ time steps, mapped size (total_via x 3)
     *  e.g. |tx(t0) ty(t0) tz(t0)|
     *       | ...    ...    ...  |
     *       |tx(tf) ty(tf) tz(tf)|
     */
    MapMatrixConst const controls(const int via_=0,
                                  int total_=0) const {
        return(MapMatrixConst(trajectory_.data()+(via_*n_cols)+1+n_state+n_angle,
                              (total_==0 ? viapoints : total_), n_control,
                              Eigen::Stride<dynamic,dynamic>(1,n_cols)) );
    }


    /**
     * @brief returns reference to row vector or matrix of row vectors
     * containing all dimensions for given time derivative, starting at via_ to
     * size of mapped is (total_via x dof) defaults to entire time sequence of
     * position.
     *  e.g. |j0(t0)  j1(t0)  ...  jn(t0)|
     *       | ...     ...    ...   ...  |
     *       |j0(tf)  j1(tf)  ...  jn(tf)|
     */
    MapMatrixConst const joint_dimensions(const int der_=position,
                                          const int via_=0,
                                          int total_=0) const {
        return( MapMatrixConst(trajectory_.data()+(via_*n_cols)+1+n_state+n_angle+n_control+der_,
                               (total_==0 ? viapoints : total_), dof,
                               Eigen::Stride<dynamic,dynamic>(ders,n_cols)) );
    }


    /**
     * @brief returns reference to row vector or matrix of row vectors
     * containing all joint derivatives for given dimension, starting at via_
     * to size of mapped
     *  e.g. |x(t0)  xd(t0) ...  xdddd(t0)|
     *       | ...    ...   ...     ...   |
     *       |x(tf)  xd(tf) ...  xdddd(tf)|
     */
    MapMatrixConst const joint_derivatives(const int dim_=0,
                                           const int via_=0,
                                           int total_=0) const {
        return( MapMatrixConst(trajectory_.data()+(via_*n_cols)+1+n_state+n_angle+n_control+(dim_*ders),
                               (total_==0 ? viapoints : total_), ders,
                               Eigen::Stride<dynamic,dynamic>(1,n_cols)) );
    }



    MapMatrixConst const manipulability(const int via_=0, int total_=0) const {
        return( MapMatrixConst(trajectory_.data()+(via_*n_cols)+1+n_state+n_angle+n_control+n_joint,
                               (total_==0 ? viapoints : total_), 1,
                               Eigen::Stride<dynamic,dynamic>(1,n_cols)) );
    }


    /**
     * @brief Sets the mapping matrix between body-fixed quaternion rate
     * body-fixed angular rate
     * TO DO: this would fail if accessed at same time by 2 processes!
     */
    void get_Wp(const Vector4& q);


    /**
     * @brief Returns any instantaneous angular rate given pose and
     * corresponding quaternion rate
     */
    Vector3 get_angular_rate(const Vector4& q, const Vector4& q_rate);


    /**
     * @brief computes omegas from quaternion rates and stores values in
     * trajectory_
     */
    void compute_angular_rates();


    /**
     * @brief
     */
    void normalize_quaternions();


    /**
     * @brief check for antipodal switches in forward time direction and
     * correct values
     */
    void check_antipodal_quaternions();

    /**
     * @brief check for antipodal switch with one quaternion
     * hack for trajectories with 1 viapoint...
     */
    void check_antipodal_quaternion(const Array4& array_q0);


    /**
     * @brief sets values of trajectory for one dimension through given
     * derivative, starting at via_ to size of total_ is (dof x total_via)
     */
    void set_derivatives(const Matrix& values,
                         const int dim_=0,
                         const int der_=position,
                         const int nder_=1,
                         const int via_=0,
                         int total_=0);


    /**
     * @brief sets values of trajectory for given derivative for given dimension
     * range, starting at via_ to size of mapped is (dof x total_via))
     * defaults to positions (der_=0)
     */
    void set_dimensions(const Matrix& values,
                        const int der_=position,
                        const int dim_=0,
                        const int ndof_=1,
                        const int via_=0,
                        int total_=0);


    /**
     * @brief sets values angular rate vector for one or more points in time,
     * where values is a row vector or matrix of row vectors of angular rates,
     * e.g. |wx(t0) wy(t0) wz(t0)|
     *      | ...    ...    ...  |
     *      |wx(tf) wy(tf) wz(tf)|
     */
    void set_angular_rates(const Matrix& values,
                           const int der_=velocity,
                           const int via_=0,
                           int total_=0);


    /**
     * @brief sets control values for 1 or multiple dof with dim_ (starting) and
     * ndof_ (total)
     */
    void set_controls(const Matrix& values,
                      const int dim_=0,
                      const int ndof_=1,
                      const int via_=0,
                      int total_=0);


    /**
     * @brief sets values of trajectory for one joint dimension through given
     * derivative, starting at via_ to size of total_.
     * \a values has size (\a dof x \a total_)
     */
    void set_joint_derivatives(const Matrix& values,
                               const int dim_,
                               const int der_,
                               const int nder_,
                               const int via_,
                               int total_);


    /**
     * @brief sets values of trajectory for given joint derivative for given
     * dimension range, starting at via_ to size of mapped is (dof x total_via))
     * defaults to positions (der_=0)
     */
    void set_joint_dimensions(const Matrix& values,
                              const int der_,
                              const int dim_,
                              const int ndof_,
                              const int via_,
                              int total_);


    /**
     * @brief set_manipulability
     * @param values
     * @param via_
     * @param total_
     */
    void set_manipulability(const Matrix& values, const int via_, int total_);



    /**
     * @brief set_manipulability
     * @param values
     * @param via_
     */
    void set_manipulability(const scalar values, const int via_);



    /**
     * @brief sets trajectory values at time t_i given a node_t
     */
    void set_node(std::shared_ptr<const node_t> node_i, const int via);


    /**
     * @brief sets node_t object with data at time t_i in trajectory
     */
    void node(std::shared_ptr<node_t> node_i, const int via_) const;


private:
    scalar t0;
    scalar tf;
    scalar step;
    int viapoints;
    MatrixRowMajor trajectory_;

    // temporary storage variable
    Matrix34 Wp;
};


} // end amp namespace
#endif //TRAJECTORIES_HPP


