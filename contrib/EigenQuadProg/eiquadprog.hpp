#ifndef _EIGEN_QUADSOLVE_HPP_
#define _EIGEN_QUADSOLVE_HPP_


/*
 FILE eiquadprog.hh

 NOTE: this is a modified of uQuadProg++ package, working with Eigen data structures.
       uQuadProg++ is itself a port made by Angelo Furfaro of QuadProg++ originally developed by
       Luca Di Gaspero, working with ublas data structures.

 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani
 for the solution of a (convex) Quadratic Programming problem
by means of a dual method.

The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^T x + ce0 = 0
    CI^T x + ci0 >= 0

 The matrix and vectors dimensions are as follows:
     G: n * n
        g0: n

        CE: n * p
     ce0: p

      CI: n * m
   ci0: m

     x: n

 The function will return the cost of the solution written in the x vector or
 std::numeric_limits::infinity() if the problem is infeasible. In the latter case
 the value of the x vector is not correct.

 References: D. Goldfarb, A. utilsally stable dual method for solving
             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

 Notes:
  1. pay attention in setting up the vectors ce0 and ci0.
       If the constraints of your problem are specified in the form
       A^T x = b and C^T x >= d, then you should set ce0 = -b and ci0 = -d.
  2. The matrix G is modified within the function since it is used to compute
     the G = L^T L cholesky factorization for further computations inside the function.
     If you need the original matrix G you should make a copy of it and pass the copy
     to the function.


 The author will be grateful if the researchers using this software will
 acknowledge the contribution of this modified function and of Di Gaspero's
 original version in their research papers.


LICENSE

Copyright (2010) Gael Guennebaud
Copyright (2008) Angelo Furfaro
Copyright (2006) Luca Di Gaspero


This file is a porting of QuadProg++ routine, originally developed
by Luca Di Gaspero, exploiting uBlas data structures for vectors and
matrices instead of native C++ array.

uquadprog is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

uquadprog is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uquadprog; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#ifdef USE_EIGEN_PKGCONFIG_HEADERS
    #include <eigen3/Eigen/Dense>
#else
    #include <Eigen/Dense>
#endif


namespace Eigen {

namespace {
    template <int E_CONSTRAINTS, int IE_CONSTRAINTS>
    struct QuadProg_traits {
        static const int Constraints = E_CONSTRAINTS+IE_CONSTRAINTS;
    };

    template <>
    struct QuadProg_traits<-1,-1> {
        static const int Constraints = -1;
    };

    template <int E_CONSTRAINTS>
    struct QuadProg_traits<E_CONSTRAINTS,-1> {
        static const int Constraints = -1;
    };

    template <int IE_CONSTRAINTS>
    struct QuadProg_traits<-1,IE_CONSTRAINTS> {
        static const int Constraints = -1;
    };

    template <int E_CONSTRAINTS>
    struct update_np{
        template <typename T1, typename T2>
        static void run(MatrixBase<T1> & np,MatrixBase<T2> const & CE_col){
            np = CE_col;
        }
    };

    template <>
    struct update_np<0> {
        template <typename T1, typename T2>
        static void run(MatrixBase<T1> & np,MatrixBase<T2> const & CE_col){
            // do nothing, CE_col has shape size,0
        }
    };
}


// namespace internal {
template<typename Scalar, int SIZE=Dynamic,
    int E_CONSTRAINTS=Dynamic, int IE_CONSTRAINTS=Dynamic>
class QuadProg {
    // import template variables
    static const int Size = SIZE;
    static const int Constraints =
        QuadProg_traits<E_CONSTRAINTS,IE_CONSTRAINTS>::Constraints;

inline Scalar distance(Scalar a, Scalar b)
{
    // TODO: replaced references to outdated Eigen::internal::{abs,sqrt}
    // with references to std::{abs,sqrt}.
    // Will not work with non-standard data types.
    Scalar a1, b1, t;
    a1 = std::abs(a);
    b1 = std::abs(b);
    if (a1 > b1)
    {
        t = (b1 / a1);
        return a1 * std::sqrt(1.0 + t * t);
    }
    else if (b1 > a1)
    {
        t = (a1 / b1);
        return b1 * std::sqrt(1.0 + t * t);
    }
    return a1 * std::sqrt(2.0);
}

// }

// TODO: use types defined above here
template <typename T1, typename T2, typename T3>
inline void compute_d(MatrixBase<T1> &d,
        const MatrixBase<T2>& J, const MatrixBase<T3>& np)
{
    d = J.adjoint() * np;
}

template <typename T1, typename T2, typename T3>
inline void update_z(MatrixBase<T1>& z,
        const MatrixBase<T2>& J,
        const MatrixBase<T3>& d,  int iq)
{
    z = J.rightCols(z.size()-iq) * d.tail(d.size()-iq);
}

template <typename T1, typename T2, typename T3>
inline void update_r(const MatrixBase<T1>& R,
        MatrixBase<T2>& r,
        const MatrixBase<T3>& d, int iq)
{
    r.head(iq)= R.topLeftCorner(iq,iq).
        template triangularView<Upper>().solve(d.head(iq));
}

//bool add_constraint(MatrixXd& R, MatrixXd& J, VectorXd& d, int& iq, double& R_norm);
//void delete_constraint(MatrixXd& R, MatrixXd& J, VectorXi& A, VectorXd& u,  int p, int& iq, int l);

public:
//typedef -> later template for fixed size
typedef Eigen::Matrix<Scalar,Size,Size> system_matrix_t;

typedef Eigen::Matrix<Scalar,E_CONSTRAINTS,1> e_constraints_vector_t;
typedef Eigen::Matrix<Scalar,Size,E_CONSTRAINTS> e_constraints_matrix_t;
typedef Eigen::Matrix<Scalar,IE_CONSTRAINTS,1> ie_constraints_vector_t;
typedef Eigen::Matrix<Scalar,Size,IE_CONSTRAINTS> ie_constraints_matrix_t;

// only internal types?
typedef Eigen::Matrix<Scalar,Constraints,1> constraints_vector_t;
typedef Eigen::Matrix<Scalar,Constraints,Size> constraints_matrix_t;

// TODO:
// fix size for internal temporary variables
typedef Eigen::Matrix<Scalar,Size,1> system_vector_t;
typedef Eigen::Matrix<int,Constraints,1> constraints_activation_t;
typedef Eigen::Matrix<Scalar,Size,1> argument_t;

// Constructor to specify sizes for the variable problem size case
QuadProg(int _size, int _e_constraints, int _ie_constraints)
    :R(_size,_size),
     J(_size,_size),
     chol(_size),
     s(_e_constraints+_ie_constraints),
     r(_e_constraints+_ie_constraints),
     u(_e_constraints+_ie_constraints),
     z(_size),d(_size),np(_size),
     x_old(_size),
     u_old(_e_constraints+_ie_constraints),
     A(_e_constraints+_ie_constraints),
     A_old(_e_constraints+_ie_constraints),
     iai(_e_constraints+_ie_constraints){
         assert(Size==Dynamic);
         assert(Constraints==Dynamic);
     };

// constructor used
// * for the static size case
// * for 0 size instanciation with a dynamic resizing in solve_quadprog
QuadProg(){};

template <typename T1, typename T2, typename T3, typename T4, typename T5,
         typename T6, typename T7>
inline double solve_quadprog(MatrixBase<T1> & G,
                             MatrixBase<T2> & g0,
                             const MatrixBase<T3> & CE,
                             const MatrixBase<T4> & ce0,
                             const MatrixBase<T5> & CI,
                             const MatrixBase<T6> & ci0,
                             MatrixBase<T7> & x)
{
    int i, k, l; /* indices */
    int ip, me, mi;
    //int n = g0.size();
    int p = ce0.size();
    int m = ci0.size();
    if(Size!=Dynamic || Constraints!=Dynamic){
        // static size case, function arguments must match class member sizes
        assert(R.rows()==G.rows());
        assert(R.cols()==G.cols());
        assert(J.rows()==G.rows());
        assert(J.cols()==G.cols());
    }else{
        // Dynamic Size class: resize class variables if necessary
        int _size = G.rows();
        int _constraints = ce0.size() + ci0.size();
        R.resize(_size,_size);
        J.resize(_size,_size);
        //chol.resize(_size);
        s.resize(_constraints);
        r.resize(_constraints);
        u.resize(_constraints);
        z.resize(_size);
        d.resize(_size);
        np.resize(_size);
        x_old.resize(_size);
        u_old.resize(_constraints);
        A.resize(_constraints);
        A_old.resize(_constraints);
        iai.resize(_constraints);
    }


    // TODO more assertion on other class member sizes and input sizes
    // assert chol.size()==G.cols();

    double f_value, psi, c1, c2, sum, ss, R_norm;
    const double inf = std::numeric_limits<double>::infinity();
    double t, t1, t2; /* t is the step length, which is the minimum of the partial step length t1
    * and the full step length t2 */
    //int q;
    int iq, iter = 0;
    bool iaexcl[m + p];

    me = p; /* number of equality constraints */
    mi = m; /* number of inequality constraints */
    //q = 0;  /* size of the active set A (containing the indices of the active constraints) */

    /*
     * Preprocessing phase
     */

    /* compute the trace of the original matrix G */
    c1 = G.trace();

    /* decompose the matrix G in the form LL^T */
    chol.compute(G);

    /* initialize the matrix R */
    d.setZero();
    R.setZero();
    R_norm = 1.0; /* this variable will hold the norm of the matrix R */

    /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
    // J = L^-T
    J.setIdentity();
    J = chol.matrixU().solve(J);
    c2 = J.trace();
#ifdef TRACE_SOLVER
    print_matrix("J", J, n);
#endif

    /* c1 * c2 is an estimate for cond(G) */

    /*
    * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x
    * this is a feasible point in the dual space
     * x = G^-1 * g0
    */
    x = chol.solve(g0);
    x = -x;
    /* and compute the current solution value */
    f_value = 0.5 * g0.dot(x);
#ifdef TRACE_SOLVER
    std::cerr << "Unconstrained solution: " << f_value << std::endl;
    print_vector("x", x, n);
#endif

    /* Add equality constraints to the working set A */
    iq = 0;
    for (i = 0; i < me; i++)
    {
        // hack to avoid Eigen assertion in case of E_CONSTRAINT==0
        update_np<E_CONSTRAINTS>::run(np,CE.col(i));
        compute_d(d, J, np);
        update_z(z, J, d,  iq);
        update_r(R, r, d,  iq);
#ifdef TRACE_SOLVER
        print_matrix("R", R, iq);
        print_vector("z", z, n);
        print_vector("r", r, iq);
        print_vector("d", d, n);
#endif

        /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint
          becomes feasible */
        t2 = 0.0;
        if (std::abs(z.dot(z)) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
            t2 = (-np.dot(x) - ce0(i)) / z.dot(np);

        x += t2 * z;

        /* set u = u+ */
        u(iq) = t2;
        u.head(iq) -= t2 * r.head(iq);

        /* compute the new solution value */
        f_value += 0.5 * (t2 * t2) * z.dot(np);
        A(i) = -i - 1;

        if (!add_constraint(R, J, d, iq, R_norm))
        {
            // FIXME: it should raise an error
            // Equality constraints are linearly dependent
            return f_value;
        }
    }

    /* set iai = K \ A */
    for (i = 0; i < mi; i++)
        iai(i) = i;

l1:
    iter++;
#ifdef TRACE_SOLVER
    print_vector("x", x, n);
#endif
    /* step 1: choose a violated constraint */
    for (i = me; i < iq; i++)
    {
        ip = A(i);
        iai(ip) = -1;
    }

    /* compute s(x) = ci^T * x + ci0 for all elements of K \ A */
    ss = 0.0;
    psi = 0.0; /* this value will contain the sum of all infeasibilities */
    ip = 0; /* ip will be the index of the chosen violated constraint */
    for (i = 0; i < mi; i++)
    {
        iaexcl[i] = true;
        sum = CI.col(i).dot(x) + ci0(i);
        s(i) = sum;
        psi += std::min(0.0, sum);
    }
#ifdef TRACE_SOLVER
    print_vector("s", s, mi);
#endif


    if (std::abs(psi) <= mi * std::numeric_limits<double>::epsilon() * c1 * c2* 100.0)
    {
        /* numerically there are not infeasibilities anymore */
        //q = iq;
        return f_value;
    }

    /* save old values for u, x and A */
    u_old.head(iq) = u.head(iq);
    A_old.head(iq) = A.head(iq);
    x_old = x;

l2: /* Step 2: check for feasibility and determine a new S-pair */
    for (i = 0; i < mi; i++)
    {
        if (s(i) < ss && iai(i) != -1 && iaexcl[i])
        {
            ss = s(i);
            ip = i;
        }
    }
    if (ss >= 0.0)
    {
        //q = iq;
        return f_value;
    }

    /* set np = n(ip) */
    np = CI.col(ip);
    /* set u = (u 0)^T */
    u(iq) = 0.0;
    /* add ip to the active set A */
    A(iq) = ip;

#ifdef TRACE_SOLVER
    std::cerr << "Trying with constraint " << ip << std::endl;
    print_vector("np", np, n);
#endif

l2a:/* Step 2a: determine step direction */
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    compute_d(d, J, np);
    update_z(z, J, d, iq);
    /* compute N* np (if q > 0): the negative of the step direction in the dual space */
    update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
    std::cerr << "Step direction z" << std::endl;
    print_vector("z", z, n);
    print_vector("r", r, iq + 1);
    print_vector("u", u, iq + 1);
    print_vector("d", d, n);
    print_ivector("A", A, iq + 1);
#endif

    /* Step 2b: compute step length */
    l = 0;
    /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
    t1 = inf; /* +inf */
    /* find the index l s.t. it reaches the minimum of u+(x) / r */
    for (k = me; k < iq; k++)
    {
        double tmp;
        if (r(k) > 0.0 && ((tmp = u(k) / r(k)) < t1) )
        {
            t1 = tmp;
            l = A(k);
        }
    }
    /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
    if (std::abs(z.dot(z))  > std::numeric_limits<double>::epsilon()) // i.e. z != 0
        t2 = -s(ip) / z.dot(np);
    else
        t2 = inf; /* +inf */

    /* the step is chosen as the minimum of t1 and t2 */
    t = std::min(t1, t2);
#ifdef TRACE_SOLVER
    std::cerr << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ";
#endif

    /* Step 2c: determine new S-pair and take step: */

    /* case (i): no step in primal or dual space */
    if (t >= inf)
    {
        /* QPP is infeasible */
        // FIXME: unbounded to raise
        //q = iq;
        return inf;
    }
    /* case (ii): step in dual space */
    if (t2 >= inf)
    {
        /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
        u.head(iq) -= t * r.head(iq);
        u(iq) += t;
        iai(l) = l;
        delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
        std::cerr << " in dual space: "
                  << f_value << std::endl;
        print_vector("x", x, n);
        print_vector("z", z, n);
        print_ivector("A", A, iq + 1);
#endif
        goto l2a;
    }

    /* case (iii): step in primal and dual space */

    x += t * z;
    /* update the solution value */
    f_value += t * z.dot(np) * (0.5 * t + u(iq));

    u.head(iq) -= t * r.head(iq);
    u(iq) += t;
#ifdef TRACE_SOLVER
    std::cerr << " in both spaces: "
              << f_value << std::endl;
    print_vector("x", x, n);
    print_vector("u", u, iq + 1);
    print_vector("r", r, iq + 1);
    print_ivector("A", A, iq + 1);
#endif

    if (t == t2)
    {
#ifdef TRACE_SOLVER
        std::cerr << "Full step has taken " << t << std::endl;
        print_vector("x", x, n);
#endif
        /* full step has taken */
        /* add constraint ip to the active set*/
        if (!add_constraint(R, J, d, iq, R_norm))
        {
            iaexcl[ip] = false;
            delete_constraint(R, J, A, u, p, iq, ip);
#ifdef TRACE_SOLVER
            print_matrix("R", R, n);
            print_ivector("A", A, iq);
#endif
            for (i = 0; i < m; i++)
                iai(i) = i;
            for (i = 0; i < iq; i++)
            {
                A(i) = A_old(i);
                iai(A(i)) = -1;
                u(i) = u_old(i);
            }
            x = x_old;
            goto l2; /* go to step 2 */
        }
        else
            iai(ip) = -1;
#ifdef TRACE_SOLVER
        print_matrix("R", R, n);
        print_ivector("A", A, iq);
#endif
        goto l1;
    }

    /* a patial step has taken */
#ifdef TRACE_SOLVER
    std::cerr << "Partial step has taken " << t << std::endl;
    print_vector("x", x, n);
#endif
    /* drop constraint l */
    iai(l) = l;
    delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
    print_matrix("R", R, n);
    print_ivector("A", A, iq);
#endif

    s(ip) = CI.col(ip).dot(x) + ci0(ip);

#ifdef TRACE_SOLVER
    print_vector("s", s, mi);
#endif
    goto l2a;
}

private:
system_matrix_t R;
system_matrix_t J;
LLT<system_matrix_t,Lower> chol;

constraints_vector_t s,r,u;
system_vector_t z,d,np;
system_vector_t x_old, u_old;
constraints_activation_t A, A_old, iai;


template <typename T1, typename T2, typename T3>
bool add_constraint(MatrixBase<T1> & R, MatrixBase<T2>& J,
        MatrixBase<T3> & d, int& iq, double& R_norm)
{
    int n=J.rows();
#ifdef TRACE_SOLVER
    std::cerr << "Add constraint " << iq << '/';
#endif
    int j, k;
    double cc, ss, h, t1, t2, xny;

    /* we have to find the Givens rotation which will reduce the element
        d(j) to zero.
        if it is already zero we don't have to do anything, except of
        decreasing j */
    for (j = n - 1; j >= iq + 1; j--)
    {
        /* The Givens rotation is done with the matrix (cc cs, cs -cc).
                 If cc is one, then element (j) of d is zero compared with element
                 (j - 1). Hence we don't have to do anything.
                 If cc is zero, then we just have to switch column (j) and column (j - 1)
                 of J. Since we only switch columns in J, we have to be careful how we
                 update d depending on the sign of gs.
                 Otherwise we have to apply the Givens rotation to these columns.
                 The i - 1 element of d has to be updated to h. */
        cc = d(j - 1);
        ss = d(j);
        h = distance(cc, ss);
        if (h == 0.0)
            continue;
        d(j) = 0.0;
        ss = ss / h;
        cc = cc / h;
        if (cc < 0.0)
        {
            cc = -cc;
            ss = -ss;
            d(j - 1) = -h;
        }
        else
            d(j - 1) = h;
        xny = ss / (1.0 + cc);
        for (k = 0; k < n; k++)
        {
            t1 = J(k,j - 1);
            t2 = J(k,j);
            J(k,j - 1) = t1 * cc + t2 * ss;
            J(k,j) = xny * (t1 + J(k,j - 1)) - t2;
        }
    }
    /* update the number of constraints added*/
    iq++;
    /* To update R we have to put the iq components of the d vector
      into column iq - 1 of R
      */
    R.col(iq-1).head(iq) = d.head(iq);
#ifdef TRACE_SOLVER
    std::cerr << iq << std::endl;
#endif

    if (std::abs(d(iq - 1)) <= std::numeric_limits<double>::epsilon() * R_norm)
        // problem degenerate
        return false;
    R_norm = std::max<double>(R_norm, std::abs(d(iq - 1)));
    return true;
}


//TODO: these methods can access the class members directly
template <typename T1, typename T2, typename T3, typename T4>
void delete_constraint(MatrixBase<T1> & R,
        MatrixBase<T2>& J, MatrixBase<T3>& A, MatrixBase<T4>& u,
        int p, int& iq, int l)
{

    int n = R.rows();
#ifdef TRACE_SOLVER
    std::cerr << "Delete constraint " << l << ' ' << iq;
#endif
    int i, j, k, qq;
    double cc, ss, h, xny, t1, t2;
	
	// add this to prevent initialization compiler warning!
	qq=0;

    /* Find the index qq for active constraint l to be removed */
    for (i = p; i < iq; i++)
        if (A(i) == l)
        {
            qq = i;
            break;
        }

    /* remove the constraint from the active set and the duals */
    for (i = qq; i < iq - 1; i++)
    {
        A(i) = A(i + 1);
        u(i) = u(i + 1);
        R.col(i) = R.col(i+1);
    }

    A(iq - 1) = A(iq);
    u(iq - 1) = u(iq);
    A(iq) = 0;
    u(iq) = 0.0;
    for (j = 0; j < iq; j++)
        R(j,iq - 1) = 0.0;
    /* constraint has been fully removed */
    iq--;
#ifdef TRACE_SOLVER
    std::cerr << '/' << iq << std::endl;
#endif

    if (iq == 0)
        return;

    for (j = qq; j < iq; j++)
    {
        cc = R(j,j);
        ss = R(j + 1,j);
        h = distance(cc, ss);
        if (h == 0.0)
            continue;
        cc = cc / h;
        ss = ss / h;
        R(j + 1,j) = 0.0;
        if (cc < 0.0)
        {
            R(j,j) = -h;
            cc = -cc;
            ss = -ss;
        }
        else
            R(j,j) = h;

        xny = ss / (1.0 + cc);
        for (k = j + 1; k < iq; k++)
        {
            t1 = R(j,k);
            t2 = R(j + 1,k);
            R(j,k) = t1 * cc + t2 * ss;
            R(j + 1,k) = xny * (t1 + R(j,k)) - t2;
        }
        for (k = 0; k < n; k++)
        {
            t1 = J(k,j);
            t2 = J(k,j + 1);
            J(k,j) = t1 * cc + t2 * ss;
            J(k,j + 1) = xny * (J(k,j) + t1) - t2;
        }
    }
}

}; // class QuadProg

} // namespace Eigen

#endif
