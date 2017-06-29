
Articulated-robot Motion Planner Library
===============================================================================

Amp provides a user-friendly API for loading robot models and generating robot 
joint and Tool Center Point (TCP) trajectories for point-to-point motion 
planning problems. Amp provides standalone trajectory parameterization and 
global search algorithms. Robot kinematic and dynamic algorithms are not a part 
of the library, but can be loaded at run-time allowing the user flexibility to 
suit their needs. The library includes general modules for: 

* Robot trajectory and state representation
    * Using Eigen matrix objects http://eigen.tuxfamily.org 
* Trajectory parameterization with general BSplines
    * Supports linear fitting
    * Supports quadratic minimization,  
* Interface to various robot kinematics and dynamics libraries:
    * RBDL https://bitbucket.org/rbdl/rbdl
    * Orocos KDL http://www.orocos.org/kdl 
* Forward and Inverse kinematics trajectory generation
    * Where IK solutions are computed numerically
* Sampling-based motion planning (RRT, RRT*)
    * Using a simple list-based tree structure and kd-tree nearest-neighbor search: https://github.com/jtsiomb/kdtree  
    * In any Tool Center Point (TCP) space
    * In the robot's joint space
    * With bound constraints on states and inputs, 
    * and collision constraints using the Bullet Collision Library http://bulletphysics.org  

The code is developed by Samantha Stoneman <samantha.stoneman@dlr.de> and 
Roberto Lampariello <roberto.lampariello@dlr.de> at the Institue of Robotics 
and Mechatronics http://www.dlr.de/rmc/rm/en at the German Aerospace Center. 
 

Installation
===============================================================================

Dependencies
-------------
Amp requires at least C++11 and is dependent on several other open-source 
libraries. The required dependencies are:

* Eigen http://eigen.tuxfamily.org (Core, Linear Algebra and Geometry modules)  
* Boost http://www.boost.org/ (odeint library)
* kdtree https://github.com/jtsiomb/kdtree 
* EigenQuadProg 

And the optional dependencies:
* RBDL https://bitbucket.org/rbdl/rbdl (rbdl, urdf/luamodel libraries)
* Orocos KDL http://www.orocos.org/kdl 
* Bullet Collision Library http://bulletphysics.org  (Collision library)

Amp will try to find the depencies using the native CMake `find_package()` tool as
well as with the `PkgConfig` module. Paths to the valid package paths for the
required dependencies and the desired optional ones should be set in the
`CMAKE_PREFIX_PATH`. Static versions of the kdtree and EigenQuadProg (header
only) are also provided in the ./contrib dir. 

Quick Installation
-------------------
Amp is built using CMake. From the package root directory, the simplest way to 
build amp using the default build options.

    > mkdir build
    > cd build && cmake ../ 
    > make 
    > make install

This will install the amp shared objects and binaries in ./lib and ./bin
directories directly in the package root directory. 

Detailed Installation
----------------------
Amp comes with several build options which affect what is built and how. The
standard cmake variable `CMAKE_INSTALL_PREFIX` allows changing the location of the
install dirs lib, bin, and incude. Amp can be built in Debug, RelWithDebInfo or
full Release mode with `CMAKE_BUILD_TYPE`.

All of the unique amp build options are detailed in the root CMakeLists.txt. The
most important of which are the options which turn on and off compilation of
code sections: 

* `DEBUG` - allows computations which are not strictly necessary e.g. matrix determinants
* `VERBOSE` - streams verbose output to std::cerr during runtime
* `LOG_ALL_DATA` - writes intermediate data to files during runtime 

The default is always OFF, setting these options to ON will greatly affect 
computation time, but can be helpful for debugging purposes. 


Recent Changes
===============================================================================



Performance
===============================================================================

The inverse kinematics -based RRT algorithm was tested on a difficult
point-to-point TCP motion which could not be solved in one trajectory. The
following table shows the performance results for several parameter settings.
This benchmark problem is described in more detail in the Examples section. 


Examples
===============================================================================

A sequence of examples are provided, increasing in complexity to demonstrate
the usage of the library. 

* A simple example for loading a robot model and environment from a .urdf file(s) and accessing model data members using the C++ API.
* An example demonstrating the inverse kinematics- based trajectory generation.
* An example demonstrating the use of the rrt planner module in a robot TCP space.


License
===============================================================================

The library is published under the GNU Lesser General Public License (LGPL).
Here is the full license text:

Amp - Articulated-robot Motion Planner library.
Copyright 2017 Samantha Stoneman
Research Associate, German Aerospace Center (DLR)

Amp is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Amp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Amp. If not, see http://www.gnu.org/licenses/.
