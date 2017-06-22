/* Unit Test for RBDL urdf file compatibility */


#include <Eigen/Dense>
#include <rbdl/rbdl.h>
#include <rbdl/addons/urdfreader/urdfreader.h>

#include <iostream>
#include <cfenv>
#include <string>
#include <boost/program_options.hpp>

#include <iomanip>
#include <string>
#include <fenv.h>

#include "luca_dynamics/model.hpp"
#include "luca_dynamics/simulation/simulation.hpp"
#include "luca_dynamics/collision/collision_model.hpp"
#include <luca_dynamics/eigen_helpers.hpp>


using namespace Eigen;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;


struct simf {
    simf() {
    #ifdef EIGEN_RUNTIME_NO_MALLOC
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    }

    void init(std::string const & filename) {
        model = luca_dynamics::create_model_from_urdf(
                    std::string(filename));
        dyn = boost::shared_ptr<luca_dynamics::luca_dynamics>(
                    new luca_dynamics::luca_dynamics(model));

        // set gravity off
        // model->set_gravity_vector(Eigen::Vector3d::Zero());
        simulation = boost::shared_ptr<simulation_t>(new simulation_t(model));
        state.resize(simulation->state_size());
        dstate.resize(simulation->state_size());
        ddstate.resize(simulation->state_size());
        tau.resize(simulation->state_size());
        state.setConstant(1.12);
        dstate.setZero();
        tau.setZero();
    }

    boost::shared_ptr<luca_dynamics::model> model;
    typedef luca_dynamics::simulation::simulation simulation_t;
    boost::shared_ptr<simulation_t> simulation;
    boost::shared_ptr<luca_dynamics::luca_dynamics> dyn;

    simulation_t::vector_t state;
    simulation_t::vector_t dstate;
    simulation_t::vector_t ddstate;
    simulation_t::vector_t tau;
};



int main (int argc, char* argv[])
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
    std::string file_name_in;

    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
    ("help", "display help")
    ("urdf-model, F", boost::program_options::value<std::string>(&file_name_in),
     "urdf model of articulated robot");
    boost::program_options::variables_map vm;
    boost::program_options::store(
                boost::program_options::parse_command_line(argc, argv, desc),
                vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    simf fx;
    fx.init(file_name_in);

    std::cerr << "  model links: " << fx.model->get_fool() << std::endl;
    std::cerr << "      model n: " << fx.model->get_n() << std::endl;
    std::cerr << "model num dof: " << fx.model->get_dof() << std::endl;

    return 0;
}
