#ifndef _AMP_SIMULATION_HPP_
#define _AMP_SIMULATION_HPP_

#include <string>

#ifdef USE_EIGEN_PKGCONFIG_HEADERS
    #include <eigen3/Eigen/Dense>
#else
    #include <Eigen/Dense>
#endif


namespace amp {


template<typename scalar>
class simulation 
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef Eigen::Matrix<scalar,3,1> vector3_t;
    typedef Eigen::Matrix<scalar,3,3> matrix3_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,1> vector_t;
    typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> matrix_t;
    typedef matrix_t jacobian_t;

    typedef Eigen::LLT<matrix_t> llt_t;
    typedef Eigen::LDLT<matrix_t> ldlt_t;
    typedef Eigen::HouseholderQR<matrix_t> qr_t;
    typedef Eigen::FullPivLU<matrix_t> lu_t;
    typedef Eigen::JacobiSVD<matrix_t> svd_t;

	/**
	 * @brief Eigen <a 
	 * href="http://eigen.tuxfamily.org/dox-devel/group__TutorialLinearAlgebra.html">
	 * linear solver</a> designation.
	 */
    enum lin_alg_t {
        DEFAULT_ALG = -1, LU = 0, QR, SVD, LLT, LDLT
    };

public:
    simulation() : 
        error_gains_(matrix_t::Zero(6,1)),
        tcp_offset_(vector_t::Zero(7)),
		tcp_id_(0) { }

	~simulation() { }

	/**
	 * @brief
	 */
	inline void set_error_gains(const scalar Kp_, const scalar Kq_) {
        error_gains_.block(0,0,3,1) = matrix_t::Constant(3,1,Kp_);
        error_gains_.block(3,0,3,1) = matrix_t::Constant(3,1,Kq_);
	}

	/**
	 * @brief
	 */ 
    inline matrix_t error_gains() const { return(error_gains_); }

    inline vector_t set_tcp_offset() const { return(tcp_offset_); }

    inline int set_tcp_id() const { return(tcp_id_); }

    inline void set_error_gains(const matrix_t& gains_in) { error_gains_ = gains_in; }

    inline void set_tcp_offset(const vector_t trafo_in) { tcp_offset_ = trafo_in; }

    inline void set_tcp_offset(double* trans_ptr, double* rot_ptr) {
        tcp_offset_.head(3) = Eigen::Map<vector_t>(trans_ptr, 3);
        tcp_offset_.tail(4) = Eigen::Map<vector_t>(rot_ptr,4);
    }
    
    inline void set_tcp_id(const int id_in) { tcp_id_ = id_in; }

	virtual const int get_model_dof() const = 0;

    /**
     * @brief Runs the forward kinematics, simplify accelerations to zero,
     * model accelerations, angular accelerations, corioli will not be
     * calculated.
     */
    virtual void fwd_kin(const vector_t& state_, const vector_t& dstate_) = 0;

    /**
     * @brief Runs the forward kinematics
     */
    virtual void fwd_dyn(vector_t& ddstate_, const vector_t& state_, 
			     const vector_t& dstate_, const vector_t& tau_) = 0;

    /**
     * @brief Runs the forward kinematics and sets current
     * tcp pose
     */
    virtual matrix_t tcp_pose(const vector_t& state_) = 0;

    virtual matrix3_t rotation_to_tcp_frame(const vector_t& state_) = 0;

	virtual matrix3_t rotation_to_base_frame(const vector_t& state_) = 0;

    /**
     * @brief Computes the manipulability
     * WARNING: performs the determinant on J
     */
    virtual scalar manipulability(const jacobian_t& J_) const = 0;

    /**
     * @brief Sets the current Jacobian
     */
    virtual jacobian_t get_J(const vector_t& state_) const = 0;
	

    virtual jacobian_t get_J_with_update(const vector_t& state_) const = 0;

    /**
     * @brief solves dstate: dtcp = J(state)*dstate
     * using one of the available linear solver
     */
    virtual vector_t inv_kin(const matrix_t& dtcp_, const vector_t& state_, 
			                 const lin_alg_t lin_=DEFAULT_ALG) = 0;

    /**
     * @brief Difference beteween desired tcp pose \a tcp_ and actual tcp
     * pose from forward kinematics, \a tcp_e
     */
    virtual matrix_t tcp_error(const matrix_t& tcp_, const vector_t& state_) = 0;

    /**
     * @brief Composes the 6x1 gain * tcp-position-error term used to
     * correct for inverse kinematics error (see Siciliano).
     */
    virtual matrix_t tcp_error_correction(const matrix_t& tcp_e) = 0;

    /**
     * @brief integrate states forward one time step
     */
    virtual vector_t diff_kin(const vector_t& state_, const vector_t& dstate_, 
			                  scalar dt_) = 0;

    /**
     * @brief runs the forward kinematics using numerically computed dstate
     */ 
	virtual matrix_t inv_kin_error_check(const vector_t& state_, 
			                             const vector_t& dstate_) = 0;

protected:
    matrix_t error_gains_;
	
    /** @brief tcp transformation from its parent joint, represented as a vector
     * where the first three elements are the relative translations and the last
     * four elements are the relative orientation expressed as the vector and
     * scalar parts of the quaterion. */
    vector_t tcp_offset_;
	
	/** @brief sequential id of tcp's parent joint */
    int tcp_id_;
};

} // namespace amp
#endif // _AMP_SIMULATION_HPP_

