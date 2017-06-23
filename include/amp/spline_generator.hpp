/*
FILE amp/spline_generator.hpp

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

#ifndef SPLINE_GENERATOR_HPP
#define SPLINE_GENERATOR_HPP


#ifdef USE_EIGEN_PKGCONFIG_HEADERS
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <eigen3/unsupported/Eigen/Splines>
#else
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <unsupported/Eigen/Splines>
#endif

#include "EigenQuadProg/eiquadprog.hpp"


#include <cmath>
#include <vector>



namespace amp {

/**
 * \ingroup amp_core
 * \class spline_generator
 * \brief A wrapper to the Eigen Splined Module to represent robot trajectories.
 *
 * Provides an interface to the Eigen Splines Module for the
 * purpose of building bspline representations of variable order,
 * variable DOF and variable parameter size specifically for
 * motion planning applications.
 *
 * NOTE:
 * The general spline formula,
 *
 * \f{align*}
 *        C_k(u) = \sum  N_i,k(u) * P_i;  for i=0, ..., n_vrtx-1, k=order
 * \f}
 *
 * \tparam scalar The underlying data type (typically float or double)
 * \tparam order_ The order of the BSpline basis f(x)'s
 * \tparam n_fbcs_ The no. fixed boundary conditions
 * \tparam n_vert_ The no. of vertices in the spline control polygon
 */
template<typename scalar, int order_=Eigen::Dynamic, int n_fbcs_=Eigen::Dynamic, int n_vert_=Eigen::Dynamic>
class spline_generator
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // always 1 for time and 1 for spline value (e.g. time and position)
    static const int dof    = 2;
    static const int dyn    = Eigen::Dynamic;
    static const int n_via_ = Eigen::Dynamic;
    static const int degree     = order_!=dyn ? order_-1 : order_;
    static const int der_order_ = order_!=dyn ? order_+1 : order_;

    // Define Templated Spline Types
    //typedef typename spline_object::KnotVectorType knot_vector;
    //typedef typename spline_object::ControlPointVectorType vertex_vector;
    typedef typename Eigen::Spline<scalar,dof,degree> spline_object;
    typedef typename Eigen::SplineTraits<spline_object,order_>::BasisDerivativeType basis_derivative;

    // NOTE:  C_k(u) = \sum  N_i,k(u) * P_i;  for i = 0,...,n_vrtx-1 and k=order
    typedef Eigen::Matrix<scalar,n_fbcs_,1> bounds;
    typedef Eigen::Matrix<scalar,1,der_order_> point_on_curve;
    typedef Eigen::Array<scalar,1,dyn> knot_vector;
    typedef Eigen::Array<scalar,dof,n_vert_> vertex_vector;
    typedef Eigen::Matrix<scalar,1,n_vert_> vertex_row;
    typedef Eigen::Matrix<scalar,dyn,1> p_vector;

    // basis Matrix storage options
    typedef Eigen::Matrix<scalar,1,n_vert_> Nut;
    typedef Eigen::Matrix<scalar,n_via_,n_vert_> Nu;
    typedef Eigen::Matrix<scalar,order_,n_vert_> Nuts;
    typedef typename std::vector<Nu,Eigen::aligned_allocator<Nu> > Nus;

    typedef Eigen::Matrix<scalar,2,1> Vector2;
    typedef Eigen::Matrix<scalar,dyn,dyn> Matrix;
    typedef Eigen::Matrix<scalar,dyn,1> Vector;
    typedef Eigen::Matrix<scalar,1,dyn> RowArray;
    typedef Eigen::Matrix<scalar,1,dyn> RowVector;
    typedef Eigen::Matrix<scalar,n_via_,dyn> curve;


    ~spline_generator();


    /**
     * \brief Automatically sets time sequence and the knot vector to be
     * uniform.
     * TODO: make different constructors for different knot vector situations.
     *
     */
    spline_generator(const int order_in,
                     const int n_fbcs_in,
                     const int n_vert_in,
                     const knot_vector& knots_in,
                     const scalar& t_final_in);


    /** \brief Overload for template instantiations with order_ provided. */
    spline_generator(const int n_fbcs_in,
                     const int n_vert_in,
                     const knot_vector& knots_in,
                     const scalar& t_final_in);



    /** \brief Overload for template instantiations with order_ and n_vert_
     * provided. */
    spline_generator(const int n_vert_in,
                     const knot_vector& knots_in,
                     const scalar& t_final_in);


    /** \brief Overload for template instantiations with order_ and n_vert_
     * provided. */
    spline_generator(const knot_vector& knots_in,
                     const scalar& t_final_in);


    /**
     * \brief Returns copy of parameter array. This is a subset of the vertex
     * vector containing only the 'free' vertices, the vertices not fully
     * constrained by boundary conditions. These are typically the so-called
     * parameters of an optimization problem.
     *
     */
    inline p_vector parameters() const {
        return(vertices_.block(1,n_fbcs,1,n_p).transpose());
    }


    vertex_vector vertices() const {
        return(vertices_);
    }


    knot_vector knots() const {
        return(knots_);
    }


    /** \brief returns real-valued copy of initial boundary conditions */
    bounds initial_condition() const {
        return((time_pwrs.block(0,0,n_fbcs,n_fbcs)).inverse() * initl_bcs);
    }


    /** \brief returns real-valued copy of final boundary conditions */
    bounds final_condition() const {
        return((time_pwrs.block(0,0,n_fbcs,n_fbcs)).inverse() * final_bcs);
    }


    /** \brief returns copy of normalized bases */
    Nus get_bases() const {
        Nus Nu_;
        Matrix time_pwrs_mat(Matrix::Zero(num_viapts, num_viapts));
        for (int j=0; j<order; ++j) {
            time_pwrs_mat.diagonal() = Vector::Constant(num_viapts, 1, (1./time_pwrs(j,j)));
            Nu_.push_back( time_pwrs_mat * N_u_vector[j] );
        }
        return (Nu_);
    }

    scalar time_final() const {
        return(t_final);
    }


    inline int sparse_index(const int t_index,const int o_index) const {
        return( spline_diagonal_idx(t_index,o_index) );
    }


    /** \brief  Set new final time and recompute normalized step */
    void set_t_final(scalar t_final_in);


    /** \brief Define number of via points and step size for x-component (aka time) */
    void set_via_points(int via_in);


    /** \brief Set basis matrices for all via points! */
    void set_bases(Nus& N_u_in);


    /** \brief Initializes vertex polygon from boundary conditions */
    void init_generator(const bounds& initl_bcs_in,
                        const bounds& final_bcs_in,
                        const p_vector& known_vertices);


    /** \brief Return entrire trajectory for specified via points */
    void set_bcs(const bounds& initl_bcs_in, const bounds& final_bcs_in);


    /** \brief Fill in the free parameters of control polygon y-dim */
    void update_parameters(const p_vector& new_parameters);


    /** \brief Fill in all coefficients of the control polygon y-dim */
    void update_vertices(const vertex_row& new_vertices);


    /**
     * \brief Sets the new values and recomputes the fixed vertices such that the new
     * boundary conditions are satisfied.
     */
    void update_ics_and_vertices(const p_vector& new_parameters,
                                 const bounds& new_initl_bcs,
                                 const bounds& new_final_bcs);


    /**
     * \brief Finds vertices to minimize squared sum of \a der derivative
     * curve values requires that N_u_vector has been precomputed
     */
    void qp_initial_guess(const int der,
                          const int num_lin_constraints,
                          const Matrix& lims);


    /**
     * \brief Finds the parameters to fit \a der derivative curve values to
     * \a trth (column Eigen vector) requires that N_u_vector has been
     * precomputed
     */
    void llsq_initial_guess(const Vector& trth, const int der=0);


    /**
     * \brief Compute basis matrices for all via points. Defaults to minimal
     * computation only for 0th derivative bases.
     */
    void compute_bases(const int der_order=1);


    /**
     * @brief compute_bases_nonuniform
     * @param knot_seq
     * @param der_order
     */
    void compute_bases_nonuniform(const Vector& knot_seq, const int der_order);

    /**
     * \brief Compute basis matrices for one time instance via points. Defaults to minimal
     * computation only for 0th derivative bases.
     */
    Nuts compute_base(const scalar u_curr, const int der_order=1);


    /**
     * \brief Return trajectory values at time instance \a via_curr and all
     * corresponding derivatives. This method computes the reqd basis functions.
     */
    curve get_one_spline_value(const scalar via_curr, const int der_order=1);



    /**
     * \brief Return trajectory values at time instance \a via_curr and all
     * corresponding derivatives, given bases available in \a N_u_vector
     */
    point_on_curve get_one_spline_value(int via_curr, const int der_order=1);


    point_on_curve get_one_spline_value(Nus& N_u_in,
                                        int via_curr,
                                        const int der_order=1);

    /**
     * \brief Return entire trajectory for specified via points and all
     * existing time derivatives
     */
    curve get_spline_values(const int der_order=1);


    /**
     * \brief Return entire trajectory for specified via points for
     * \a der_order number of time derivatives.
     */
    curve get_spline_values(Nus& N_u_in, const int der_order=1);


private:
    /** \brief Build matrices for fitting functions */
    void linear_problem_setup(Matrix& A,
                              Vector& b,
                              const Matrix& Nspl,
                              const Vector& Pspl);

    /**
     * \brief Compute the unkknown polygon vertices for the case
     * \a n_fbcs = \a order.
     */
    void compute_unknown_vertices_symmetric();


    /**
     * \brief Compute the unkknown polygon vertices for the case
     * \a n_fbcs < \a order
     */
    void compute_unknown_vertices_asymmetric();


    const int order;
    const int n_vert;
    const int n_fbcs;
    const int n_p;

    /** \brief column vector of initial conditions */
    bounds initl_bcs;
    bounds final_bcs;

    /** \brief heap-allocated matrix of bases for entire curve */
    Nus N_u_vector;

    int num_viapts;
    const int num_knots;

    scalar t_final;
    scalar delta_tn;

    knot_vector knots_;
    vertex_vector vertices_;
    spline_object spline_;

    Eigen::Matrix<scalar,order_,order_> time_pwrs;
    Eigen::Matrix<int,dyn,2> spline_diagonal_idx;
};


} // end  namespace amp

#endif // SPLINE_GENERATOR_HPP
