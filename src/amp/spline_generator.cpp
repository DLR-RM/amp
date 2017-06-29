/*
FILE amp/spline_generator.cpp

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

#include "amp/spline_generator.hpp"


namespace amp {




template<typename scalar, int order_, int n_fbcs_, int n_vert_>
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::spline_generator(const int order_in,
                   const int n_fbcs_in,
                   const int n_vert_in,
                   const knot_vector& knots_in,
                   const scalar& t_final_in) :
    order(order_in),
    n_vert(n_vert_in),
    n_fbcs(n_fbcs_in),
    n_p(n_vert-2*n_fbcs),
    initl_bcs(bounds::Zero(n_fbcs)),
    final_bcs(bounds::Zero(n_fbcs)),
    num_viapts(0),
    num_knots(n_vert+order),
    t_final(t_final_in),
    delta_tn(0.0),
    knots_(knots_in),
    vertices_(vertex_vector::Zero(2,n_vert)),
    time_pwrs(Eigen::Matrix<scalar,order_,order_>::Zero(order,order))
{
    // Define the Time and Knot Sequences
    RowArray temp = RowArray::LinSpaced(n_vert, 0., t_final);
    vertices_.block(0,0,1,n_vert) = temp;

    // Define the time normalization matrix
    scalar bar = 0.0;
    for(int i=0; i<order; i++) {
        time_pwrs(i,i) = pow(t_final, bar);
        bar+=1.0;
    }

    // Create the Spline Object
    spline_ = spline_object(knots_, vertices_);
}



template<typename scalar, int order_, int n_fbcs_, int n_vert_>
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::spline_generator(const int n_fbcs_in,
                   const int n_vert_in,
                   const knot_vector& knots_in,
                   const scalar& t_final_in) :
    order(order_),
    n_vert(n_vert_in),
    n_fbcs(n_fbcs_in),
    n_p(n_vert-2*n_fbcs),
    initl_bcs(bounds::Zero(n_fbcs)),
    final_bcs(bounds::Zero(n_fbcs)),
    num_viapts(0),
    num_knots(n_vert+order),
    t_final(t_final_in),
    delta_tn(0.0),
    knots_(knots_in),
    vertices_(vertex_vector::Zero(2,n_vert)),
    time_pwrs(Eigen::Matrix<scalar,order_,order_>::Zero(order,order))
{
    // Define the Time and Knot Sequences
    RowArray temp = RowArray::LinSpaced(n_vert, 0., t_final);
    vertices_.block(0,0,1,n_vert) = temp;

    // Define the time normalization matrix
    scalar bar = 0.0;
    for(int i=0; i<order; i++) {
        time_pwrs(i,i) = pow(t_final, bar);
        bar+=1.0;
    }

    // Create the Spline Object
    spline_ = spline_object(knots_, vertices_);
}



template<typename scalar, int order_, int n_fbcs_, int n_vert_>
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::spline_generator(const int n_vert_in,
                   const knot_vector& knots_in,
                   const scalar& t_final_in) :
    order(order_),
    n_vert(n_vert_in),
    n_fbcs(n_fbcs_),
    n_p(n_vert-2*n_fbcs),
    initl_bcs(bounds::Zero(n_fbcs)),
    final_bcs(bounds::Zero(n_fbcs)),
    num_viapts(0),
    num_knots(n_vert+order),
    t_final(t_final_in),
    delta_tn(0.0),
    knots_(knots_in),
    vertices_(vertex_vector::Zero(2,n_vert)),
    time_pwrs(Eigen::Matrix<scalar,order_,order_>::Zero(order,order))
{
    // Define the Time and Knot Sequences
    RowArray temp = RowArray::LinSpaced(n_vert, 0., t_final);
    vertices_.block(0,0,1,n_vert) = temp;

    // Define the time normalization matrix
    scalar bar = 0.0;
    for(int i=0; i<order; i++) {
        time_pwrs(i,i) = pow(t_final, bar);
        bar+=1.0;
    }

    // Create the Spline Object
    spline_ = spline_object(knots_, vertices_);
}



template<typename scalar, int order_, int n_fbcs_, int n_vert_>
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::spline_generator(const knot_vector& knots_in,
                   const scalar& t_final_in) :
    order(order_),
    n_vert(n_vert_),
    n_fbcs(n_fbcs_),
    n_p(n_vert-2*n_fbcs),
    initl_bcs(bounds::Zero(n_fbcs)),
    final_bcs(bounds::Zero(n_fbcs)),
    num_viapts(0),
    num_knots(n_vert+order),
    t_final(t_final_in),
    delta_tn(0.0),
    knots_(knots_in),
    vertices_(vertex_vector::Zero(2,n_vert)),
    time_pwrs(Eigen::Matrix<scalar,order_,order_>::Zero(order,order))
{
    // Define the Time and Knot Sequences
    RowArray temp = RowArray::LinSpaced(n_vert, 0., t_final);
    vertices_.block(0,0,1,n_vert) = temp;

    // Define the time normalization matrix
    scalar bar = 0.0;
    for(int i=0; i<order; i++) {
        time_pwrs(i,i) = pow(t_final, bar);
        bar+=1.0;
    }

    // Create the Spline Object
    spline_ = spline_object(knots_, vertices_);
}




template<typename scalar, int order_, int n_fbcs_, int n_vert_>
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::~spline_generator()
{
}




template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::set_via_points(int via_in)
{
    // via points
    num_viapts = via_in;

    // normalized step in time dimesion
    delta_tn = 1.0/(static_cast<scalar>(num_viapts-1));

    spline_diagonal_idx = Eigen::Matrix<int,dyn,2>::Zero(num_viapts,2);
}







template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::set_t_final(scalar t_final_in)
{
    t_final = t_final_in;

    // 'x-dimension' of vertices
    vertices_.template block<1,n_vert_>(0,0,1,n_vert).setLinSpaced(0.0,t_final);

    // 2. Define the time normalization matrix
    scalar bar = 0.0;
    for(int i=0; i<order; i++) {
        time_pwrs(i,i) = pow(t_final, bar);
        bar+=1.0;
    }
}







template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::init_generator(const bounds& initl_bcs_in,
                 const bounds& final_bcs_in,
                 const p_vector& known_parameters)
{
    // Fill in the known Vertices (dim Y in spline world)
    update_parameters(known_parameters);

    // Set the Boundary Conditions scaled to the total trajectory time span
    set_bcs(initl_bcs_in, final_bcs_in);
}








template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::set_bcs(const bounds& initl_bcs_in, const bounds& final_bcs_in)
{
    // Set the Boundary Conditions and scale to the total trajectory time span
    initl_bcs = time_pwrs.block(0,0,n_fbcs,n_fbcs) * initl_bcs_in;
    final_bcs = time_pwrs.block(0,0,n_fbcs,n_fbcs) * final_bcs_in;

    // Set Unknown Vertices
    if ( n_fbcs == order ) {
        compute_unknown_vertices_symmetric();
    }
    else if ( n_fbcs < order ) {
        compute_unknown_vertices_asymmetric();
    }
}






template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::update_parameters(const p_vector& new_parameters)
{
    vertices_.block(1,n_fbcs,1,n_p) = new_parameters.transpose();
}





template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::update_vertices(const vertex_row& new_vertices)
{
    vertices_.block(1,0,1,n_vert) = new_vertices;
}






template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::update_ics_and_vertices(const p_vector& new_parameters,
                          const bounds& new_initl_bcs,
                          const bounds& new_final_bcs)
{
    update_parameters(new_parameters);
    set_bcs(new_initl_bcs, new_final_bcs);
}






template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::compute_unknown_vertices_symmetric()
{
    // local mem storage
    Vector vertices_tmp = Vector::Zero(n_fbcs,1);
    Matrix NMatPartial = Matrix::Zero(n_fbcs,n_fbcs);

    // Get the basis function derivatives for initial segments
    NMatPartial = spline_.basisFunctionDerivatives(knots_(0), n_fbcs);

    // compute the vertices
    vertices_tmp = (NMatPartial.inverse())*(initl_bcs);
    vertices_.block(1,0,1,n_fbcs) = vertices_tmp.transpose();

    // Get the basis function derivatives for final segments
    NMatPartial = spline_.basisFunctionDerivatives(knots_(num_knots-1), n_fbcs);

    // compute the vertices
    vertices_tmp = (NMatPartial.inverse())*(final_bcs);
    vertices_.block(1,(n_vert-n_fbcs),1,n_fbcs) = vertices_tmp.transpose();
}







template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::compute_unknown_vertices_asymmetric()
{

    // simplify names for clarity
    // no. of fixed conditions (== n_fbcs)
    const int m = n_fbcs;

    // no. of free conditions (== order - n_fbcs)
    const int h = order-n_fbcs; // no. of free conditions (== order - n_fbcs)

    // 1. local mem storage
    Matrix vertices_tmp =Matrix::Zero(m,1);
    Matrix NMatPartial = Matrix::Zero(m,m);
    Matrix NMatId = Matrix::Identity(m,m);

    Matrix Ncurr = Matrix::Zero(order,order);
    Matrix Ma = Matrix::Zero(m,m*h);
    Matrix Mb = Matrix::Zero(m*h,1);

    // 2. Get the basis function derivatives for initial segments
    Ncurr = spline_.basisFunctionDerivatives(knots_(0), order);
    NMatPartial = Ncurr.block(0,0,m,m);

    // 3. compute the vertices
    for (int i=0; i<h; i++) {
        Ma.block(0,i*m,m,m) = vertices_(1,i+m) * NMatId;
        Mb.block(i*m,0,m,1) = Ncurr.block(0,i+m,m,1);
    }
    vertices_tmp = (NMatPartial.inverse()) * (initl_bcs - (Ma * Mb));
    vertices_.block(1,0,1,m) = vertices_tmp.transpose();


    // 4. Get the basis function derivatives for final segments
    Ncurr = spline_.basisFunctionDerivatives(knots_(num_knots-1), order);
    NMatPartial = Ncurr.block(0,0+h,m,m);
    for (int i=0; i<h; i++) {
        Ma.block(0,i*m,m,m) = vertices_(1,n_vert-(m+h-i)) * NMatId;
        Mb.block(i*m,0,m,1) = Ncurr.block(0,i,m,1);
    }
    // 5. compute the vertices
    vertices_tmp = (NMatPartial.inverse()) * (final_bcs - (Ma * Mb));
    vertices_.block(1,n_vert-m,1,m) = vertices_tmp.transpose();
}




template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::qp_initial_guess(const int der,
                   const int num_lin_constraints,
                   const Matrix &lims)
{
    // set up matrices
    Matrix Ntotl = Matrix::Zero(num_viapts,n_vert);
    Vector Ptotl = Vector::Zero(n_vert);
    Matrix Nfree = Matrix::Zero(num_viapts,n_p);
    Vector bfree = Vector::Zero(num_viapts);

    Matrix rowv = Matrix::Zero(1,n_p);
    Matrix colv = Matrix::Zero(n_p,1);
    Vector f = Vector::Zero(n_p);
    Matrix H = Matrix::Zero(n_p,n_p);

    // spline vertices & basis matrix
    Ptotl = vertices_.block(1,0,1,n_vert).transpose();
    Ntotl = N_u_vector[der];
    Ntotl *= (1./time_pwrs(der,der));
    linear_problem_setup(Nfree, bfree, Ntotl, Ptotl);

    f = (bfree.transpose() * Nfree).transpose();
    for (int i=0; i<num_viapts; i++) {
        colv = (Nfree.block(i,0,1,n_p)).transpose();
        rowv = Nfree.block(i,0,1,n_p);
        H += colv * rowv;
    }

    // num inequality constraints
    int m = 0;
    m = num_viapts*(2*num_lin_constraints);
    p_vector x;
    typedef typename Eigen::QuadProg<scalar> qp_t;
    typename qp_t::e_constraints_matrix_t CE;
    typename qp_t::e_constraints_vector_t ce0;
    typename qp_t::ie_constraints_matrix_t CI;
    typename qp_t::ie_constraints_vector_t ci0;
    qp_t qp;

    // position envelope constraint
    CI = qp_t::ie_constraints_matrix_t::Zero(n_p, m);
    ci0 = qp_t::ie_constraints_vector_t::Zero(m);
    for (size_t i=0; i<num_lin_constraints; i++)
    {
        Ntotl = N_u_vector[i];
        Ntotl *= (1./time_pwrs(i,i));
        linear_problem_setup(Nfree, bfree, Ntotl, Ptotl);

        // lower bounds
        ci0.segment((i*2)*num_viapts,num_viapts)
                = -1.0*(Vector::Constant(num_viapts, lims(i,0)) - bfree);
        // upper bounds
        ci0.segment((i*2+1)*num_viapts, num_viapts)
                = Vector::Constant(num_viapts, lims(i,1)) - bfree;

        // lower bounds
        CI.block(0,(i*2)*num_viapts,n_p,num_viapts)
                = Nfree.transpose();
        // upper bounds
        CI.block(0,(i*2+1)*num_viapts,n_p,num_viapts)
                = -1.0*Nfree.transpose();
    }

    qp.solve_quadprog(H, f, CE, ce0, CI, ci0, x);

    // Set vertices equal to solution
    update_parameters(x);
}






template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::llsq_initial_guess(const Vector& trth, const int der)
{
    // Set Up Matrices
    p_vector x;
    Matrix Ntotl = Matrix::Zero(num_viapts,n_vert);
    Vector Ptotl = Vector::Zero(n_vert);
    Matrix Nfree = Matrix::Zero(num_viapts,n_p);
    Vector bfree = Vector::Zero(num_viapts);

    // spline data
    Ptotl = vertices_.block(1,0,1,n_vert).transpose();
    Ntotl = N_u_vector[der];
    Ntotl *= (1.0/time_pwrs(der,der));

    linear_problem_setup(Nfree, bfree, Ntotl, Ptotl);
    x = Nfree.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve((trth-bfree));

    // Set Vertices to solution
    update_parameters(x);
}







template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::linear_problem_setup(Matrix& A, Vector& b, const Matrix& Nspl, const Vector& Pspl)
{
    Matrix Nfixd = Matrix::Zero(num_viapts,2*n_fbcs);
    Vector Pfixd = Vector::Zero(2*n_fbcs);

    Nfixd.block(0,0,num_viapts,n_fbcs) = Nspl.block(0,0,num_viapts,n_fbcs);
    Nfixd.block(0,n_fbcs,num_viapts,n_fbcs) = Nspl.block(0,n_fbcs+(n_p),num_viapts,n_fbcs);
    Pfixd.head(n_fbcs) = Pspl.head(n_fbcs);
    Pfixd.tail(n_fbcs) = Pspl.tail(n_fbcs);

    A = Nspl.block(0,n_fbcs,num_viapts,n_p);
    b = Nfixd*Pfixd;
}






template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::set_bases(Nus& N_u_in)
{
    N_u_vector.clear();
    N_u_vector = N_u_in;
}





template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::compute_bases(const int der_order)
{
    N_u_vector.clear();
    Eigen::DenseIndex n = 0;
    scalar iter = 0.0;
    scalar u_curr = 0.0;

    // fill in with correct sized mats
    for (int j=0; j<der_order; j++) {
        N_u_vector.push_back(Nu::Zero(num_viapts,n_vert));
    }

    // The basis vectors
    basis_derivative NTmp(der_order,order);

    for(int i=0; i<num_viapts; i++) {
        // Get a point on the curve at every viapoint
        u_curr = iter*(delta_tn);
        n = spline_.span(u_curr);
        NTmp = spline_.basisFunctionDerivatives(u_curr,der_order);

        spline_diagonal_idx(i,0) = (n+1) - order - n_fbcs;
        spline_diagonal_idx(i,1) = n-n_fbcs;
        for (int j=0; j<der_order; j++) {
            N_u_vector[j].block(i,(n+1)-order,1,order) = NTmp.block(j,0,1,order);
        }
        iter+=1.0;
    }
}





template<typename scalar, int order_, int n_fbcs_, int n_vert_>
void spline_generator<scalar,order_,n_fbcs_,n_vert_>
::compute_bases_nonuniform(const Vector& knot_seq, const int der_order)
{
    N_u_vector.clear();
    Eigen::DenseIndex n = 0;
    scalar iter = 0.0;
    scalar u_curr = 0.0;

    // fill in with correct sized mats
    for (int j=0; j<der_order; j++) {
        N_u_vector.push_back(Nu::Zero(num_viapts,n_vert));
    }

    // The basis vectors
    basis_derivative NTmp(der_order,order);

    for(int i=0; i<knot_seq.rows(); i++) {
        // Get a point on the curve at every viapoint
        u_curr = knot_seq(i) * delta_tn;
        n = spline_.span(u_curr);
        NTmp = spline_.basisFunctionDerivatives(u_curr,der_order);

        spline_diagonal_idx(i,0) = (n+1) - order - n_fbcs;
        spline_diagonal_idx(i,1) = n-n_fbcs;
        for (int j=0; j<der_order; j++) {
            N_u_vector[j].block(i,(n+1)-order,1,order) = NTmp.block(j,0,1,order);
        }
        iter+=1.0;
    }
}


template<typename scalar, int order_, int n_fbcs_, int n_vert_>
typename spline_generator<scalar,order_,n_fbcs_,n_vert_>::Nuts
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::compute_base(const scalar u_curr, const int der_order)
{
    Nuts bases = Nuts::Zero(order,n_vert);
    Eigen::DenseIndex n = 0;
    n = spline_.span(u_curr);
    basis_derivative NTmp = spline_.basisFunctionDerivatives(u_curr,der_order);
    for (int j=0; j<der_order; j++) {
        bases.block(j,(n+1)-order,1,order) = NTmp.block(j,0,1,order);
    }
    return(bases);
}





template<typename scalar, int order_, int n_fbcs_, int n_vert_>
typename spline_generator<scalar,order_,n_fbcs_,n_vert_>::curve
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::get_one_spline_value(const scalar via_curr, const int der_order)
{
    // Compute the basis
    scalar normed_time = 0.;
    normed_time = via_curr / t_final;
    Nuts bases = compute_base(normed_time,der_order+1);

    // The Spline curve
    //point_on_curve CMat = point_on_curve::Zero(1,der_order+1);

    curve CMat(curve::Zero(1, der_order+1));


    // Get Curve data and scale to real values
    CMat(0,0) = via_curr;
    for (int j=0; j<der_order; ++j) {
        CMat.block(0,j+1,1,1) = bases.block(j,0,1,n_vert) *
                (vertices_.block(1,0,1,n_vert)).transpose().matrix();
        CMat.block(0,j+1,1,1) *= (1.0/time_pwrs(j,j));
    }
    return(CMat);
}




template<typename scalar, int order_, int n_fbcs_, int n_vert_>
typename spline_generator<scalar,order_,n_fbcs_,n_vert_>::point_on_curve
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::get_one_spline_value(int via_curr, const int der_order)
{
    // The Spline curve
    point_on_curve CMat = point_on_curve::Zero(der_order+1);

    // Get Curve data and scale to real values
    CMat(0,0) = static_cast<scalar>(via_curr)*(delta_tn*t_final);
    for (int j=0; j<der_order; ++j) {
        CMat.block(0,j+1,1,1) = (N_u_vector[j]).block(via_curr,0,1,n_vert) *
                (vertices_.block(1,0,1,n_vert)).matrix().transpose();
        CMat.block(0,j+1,1,1) *= (1.0/time_pwrs(j,j));
    }
    return(CMat);
}





template<typename scalar, int order_, int n_fbcs_, int n_vert_>
typename spline_generator<scalar,order_,n_fbcs_,n_vert_>::point_on_curve
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::get_one_spline_value(Nus& N_u_in,
                       int via_curr,
                       const int der_order)
{
    // The Spline curve
    point_on_curve CMat = point_on_curve::Zero(der_order+1);

    // Get Curve data and scale to real values
    CMat(0,0) = static_cast<scalar>(via_curr)*(delta_tn*t_final);
    for (int j=0; j<der_order; j++) {
        CMat.block(0,j+1,1,1) = N_u_in[j].block(via_curr,0,1,n_vert) *
                (vertices_.block(1,0,1,n_vert)).matrix().transpose();
        CMat.block(0,j+1,1,1) *= (1.0/time_pwrs(j,j));
    }
    return(CMat);
}





template<typename scalar, int order_, int n_fbcs_, int n_vert_>
typename spline_generator<scalar,order_,n_fbcs_,n_vert_>::curve
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::get_spline_values(const int der_order)
{
    // The Spline curve
    curve CMat(curve::Zero(num_viapts, der_order+1));

    // Get Curve data and scale to real values
    Matrix time_pwrs_mat(Matrix::Zero(num_viapts, num_viapts));
    scalar iter = 0.0;
    for(int i=0; i<num_viapts; i++) {
        CMat(i,0) = iter*(delta_tn)*t_final;
        iter+=1.0;
    }
    for (int j=0; j<der_order; j++) {
        time_pwrs_mat.diagonal() =
                Vector::Constant(num_viapts, 1, (1./time_pwrs(j,j)));
        CMat.block(0,j+1,num_viapts,1) = time_pwrs_mat *
                (N_u_vector[j]*(vertices_.block(1,0,1,n_vert)).matrix().transpose());
    }
    return(CMat);
}





template<typename scalar, int order_, int n_fbcs_, int n_vert_>
typename spline_generator<scalar,order_,n_fbcs_,n_vert_>::curve
spline_generator<scalar,order_,n_fbcs_,n_vert_>
::get_spline_values(Nus& N_u_in, const int der_order)
{
    // the spline curve
    curve CMat(curve::Zero(num_viapts, der_order+1));

    // get curve data and scale to real values
    Matrix time_pwrs_mat(Matrix::Zero(num_viapts, num_viapts));
    scalar iter = 0.0;
    for(int i=0; i<num_viapts; i++) {
        CMat(i,0) = iter*(delta_tn)*t_final;
        iter+=1.0;
    }
    for (int j=0; j<der_order; j++) {
        time_pwrs_mat.diagonal() = Vector::Constant(num_viapts, 1, (1./time_pwrs(j,j)));
        CMat.block(0,j+1,num_viapts,1) = time_pwrs_mat *
                (N_u_in[j] * (vertices_.block(1,0,1,n_vert)).matrix().transpose());
    }
    return(CMat);
}

// instantiate some templates
template class spline_generator<float>;
template class spline_generator<float,4>;
template class spline_generator<float,4,2>;

template class spline_generator<double>;
template class spline_generator<double,4>;
template class spline_generator<double,4,2>;

} // end namespace amp
