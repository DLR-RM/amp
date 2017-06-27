
Articulated-robot Motion Planner Library
===============================================================================

This is the documentation for Amp, the Articulated-robot Motion Planner
library. Amp provides a user-friendly API for loading robot models and 
generating robot joint and Tool Center Point (TCP) trajectories for point-
to-point motion planning problems. Amp provides standalone trajectory
parameterization and global search algorithms. Robot kinematic and
dynamic algorithms are not a part of the library, but can be loaded at
run-time allowing the user flexibility to suit their needs. The library 
includes general modules for: 

* Robot trajectory and state representation
    * Using Eigen matrix objects
* Trajectory parameterization with general BSplines
    * Supports linear fitting
    * Supports quadratic minimization 
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

Coming soon.


Recent Changes
===============================================================================

* 22 June 2017: New release 1.0.0:
    * Initial public release.


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
