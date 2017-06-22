#ifndef NODES_HPP
#define NODES_HPP


#ifdef USE_EIGEN_PKGCONFIG_HEADERS
  #include <eigen3/Eigen/Dense>
#else
  #include <Eigen/Dense>
#endif

#include <memory>


namespace amp {

/**
 * @brief The node class holds the state of a robot with dofj number of joints
 * and one defined tcp for one instant in time.
 *
 * \a doft -- degrees of translational freedom of tcp
 * \a dofr -- degrees of rotational freedom (+1 for quaternion rep.) of tcp
 * \a doa  -- degrees of freedom in actuation of tcp
 * \a dofj -- degrees of freedom of robot == number of joints
 * \a ders -- number of time derivatives in state (NOTE: counting from 1)
 */
template<typename scalar>
class  node {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    static const int doft = 3;
    static const int dofr = 4;
    static const int dof = doft+dofr;
    const int doa;
    const int dofj;
    const int ders;

    static const int dynamic = Eigen::Dynamic;
    typedef Eigen::Quaternion<scalar> Quaternion;
    typedef Eigen::Matrix<scalar,dynamic,dynamic> state_t; // (dof x ders)
    typedef Eigen::Matrix<scalar,dynamic,1> control_t;     // (doa x ders)
    typedef Eigen::Matrix<scalar,3,dynamic> angular_t;     // (3 x ders)
    typedef Eigen::Matrix<scalar,dynamic,dynamic> joint_state_t; // (dofj x ders)


    /**
     * @brief sets the runtime number of joints and desired time derivatives,
     * and all data members to zero
     */
    node(const int dofj_, const int ders_);


    /**
     * @brief instantiates new node with all the data and runtime sizes of
     * \a m_
     */
    node(const std::shared_ptr<const node<scalar> >& m_);


    /**
     * @brief Copy all data members of \a m_, assumes runtime sizes are equal!
     */
    void copy_data(const std::shared_ptr<const node<scalar> >& m_);

    scalar time_curr;


    /**
     * @brief tcp state
     */
    state_t state;


    /**
     * @brief tcp forces and torques
     */
    control_t control;


    /**
     * @brief tcp angular derivatives defined in the body frame
     */
    angular_t angular;


    /**
     * @brief joint states
     */
    joint_state_t joints;
};

} //end amp namespace

#endif //NODES_HPP


