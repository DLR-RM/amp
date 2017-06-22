#ifndef _INTERPOLATE_H_
#define _INTERPOLATE_H_

#include "amp/nodes.hpp"
#include "amp/trajectories.hpp"
#include "amp/spline_generator.hpp"

#ifdef USE_EIGEN_PKGCONFIG_HEADERS
    #include <eigen3/Eigen/Dense>
#else
    #include <Eigen/Dense>
#endif

#include <list>
#include <string>
#include <iostream>
#include <iomanip>
#include <memory>

#ifdef DEBUG
  #define _DEBUG_
#endif


namespace amp {


template <typename scalar>
class  interpolator
{
public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    static const int dynamic = Eigen::Dynamic;

    // flags for guess type
    enum guess_flag_t {
        DEFAULT_GUESS = 0, // set guess using default QP min. Jerk Solution
        NOT_SET = 0,       // set guess using default QP min. Jerk Solution
        DIRECT_PARAMETERS, // with directly passed parameters
        DIRECT_TRAJECTORY  // with directly passed trajectory
    };

    enum edge_flag_t {
        DEFAULT_EDGE = 0,
        IS_INTERMEDIATE = 0,
        IS_FINAL
    };

    typedef node<scalar> node_t;
    typedef trajectory<scalar> trajectory_t;

    //! TODO:
    //! this template instantiation is hardcoded in spline cpp file -
    //! how to make that more general??
    static const int order  = 4;
    typedef spline_generator<scalar,order> sgt;

    typedef Eigen::Matrix<bool,1,9> failure_modes_t;
    typedef Eigen::Matrix<scalar,2,1> Vector2;
    typedef Eigen::Matrix<scalar,3,1> Vector3;
    typedef Eigen::Matrix<scalar,3,3> Matrix33;
    typedef Eigen::Quaternion<scalar> Quaternion;
    typedef Eigen::Array<scalar,4,1> Array4;
    typedef Eigen::Matrix<scalar,dynamic,1> Vector;
    typedef Eigen::Matrix<scalar,1,dynamic> rowVector;
    typedef Eigen::Matrix<scalar,dynamic,dynamic> Matrix;
    typedef Eigen::Matrix<scalar,dynamic,dynamic,Eigen::RowMajor> MatrixRowMajor;
    typedef Eigen::Map<Matrix,0,Eigen::Stride<dynamic,dynamic> > MapMatrix;
    typedef Eigen::Map<MatrixRowMajor,0,Eigen::Stride<dynamic,dynamic> > MapRowMajorMatrix;
    typedef Eigen::Map<const Matrix,0,Eigen::Stride<dynamic,dynamic> > MapMatrixConst;

    //! holds detailed optimization results
    struct result_t {
        scalar cost;
        Vector cost_norm_constants;
        Matrix parameters;
        std::shared_ptr<const node_t> node_start;
        std::shared_ptr<const node_t> node_end;
        std::shared_ptr<trajectory_t> traj_result;

        result_t(const unsigned int cost_size,
                 const int dofj_,
                 const int ders_,
                 const unsigned int n_p_,
                 const unsigned int n_pd_) :
            cost(0.),
            cost_norm_constants(Vector::Zero(cost_size)),
            parameters(Matrix::Zero(n_p_,n_pd_)),
            node_start(std::make_shared<node_t>(dofj_,ders_)),
            node_end(std::make_shared<node_t>(dofj_,ders_))
        { }

        result_t(const unsigned int cost_size,
                 const unsigned int n_p_,
                 const unsigned int n_pd_,
                 std::shared_ptr<const node_t> n_start,
                 std::shared_ptr<const node_t> n_end) :
            cost(0.),
            cost_norm_constants(Vector::Zero(cost_size)),
            parameters(Matrix::Zero(n_p_,n_pd_)),
            node_start(n_start),
            node_end(n_end)
        { }
    };


    interpolator(const int via_points_,
                 const int dofj_,
                 const int ders_,
                 const guess_flag_t guess_type_in=DEFAULT_GUESS);

    ~interpolator();

    void set_num_fixed_bcs(const int n_fbcs_in) {
        n_fbcs = n_fbcs_in;
    }

    std::shared_ptr<trajectory_t> resample(Matrix& parameters_, const int via_points_in);

    void sample_one_time(std::shared_ptr<trajectory_t>& traj_, const scalar via_points_in);

    void set_edge_flag(const edge_flag_t& edge_flag_in) {
        edge_flag = edge_flag_in;
    }

    void set_guess_flag(const guess_flag_t& guess_flag_in) {
        guess_flag = guess_flag_in;
    }

	void set_num_params(const int n_p_in);
    
    void set_via_points(const int via_points_in);

    std::shared_ptr<result_t> run(const std::shared_ptr<const node_t>& state_from_in,
                                  const std::shared_ptr<const node_t>& state_to_in,
                                  const int n_p_);

    std::shared_ptr<result_t> run(std::shared_ptr<trajectory_t>& traj_guess,
                                  const std::shared_ptr<const node_t>& state_from_in,
                                  const std::shared_ptr<const node_t>& state_to_in,
                                  const int n_p_);

    scalar time_final() const { return(time_final_); }

    scalar time_initial() const { return(time_initial_); }

    scalar duration() const { return(duration_); }

private:
    edge_flag_t edge_flag;
    guess_flag_t guess_flag;

    int n_fbcs;
    const int dofj;
    const int ders;
    int n_p;
    int n_pd;
    int via_points;
    scalar duration_;
    scalar time_step;
    scalar time_final_;
    scalar time_initial_;
    std::shared_ptr<const node_t> state_to;
    std::shared_ptr<const node_t> state_from;

    Matrix guess_parameters;
    std::shared_ptr<trajectory_t> guess_trajectory;

    // keep (dof-1) spline objects in memory
    std::list<std::shared_ptr<sgt> > splines;


private:
    void set_knots(typename sgt::knot_vector& knots);

    void set_spline_parameterization(Matrix& parameters,
                                     std::shared_ptr<trajectory_t>& traj);

    void get_spline_parameterization(Matrix& coeff);

    void set_guess_with_qp(Matrix& parameters_,
                           std::shared_ptr<trajectory_t>& trajectory_,
                           const typename sgt::knot_vector& nuknots);

    void set_guess_from_parameters(const Matrix& parameters_,
                                   std::shared_ptr<trajectory_t>& trajectory_,
                                   const typename sgt::knot_vector& nuknots);

    void set_guess_from_trajectory(Matrix& parameters_,
                                   std::shared_ptr<trajectory_t>& trajectory_,
                                   const typename sgt::knot_vector& nuknots);
};


} // end namespace amp
#endif //_INTERPOLATE_H_
