#################################
#Build the rbdl interface library
#################################
if(rbdl_FOUND)
    link_directories(${rbdl_LIBRARY_DIRS})
    set(rbdl_itf_LIBS rbdl rbdl_urdfreader)
    set(rbdl_itf_SOURCES amp/simulation_rbdl.cpp)

    add_library(rbdl_interface SHARED ${rbdl_itf_SOURCES})
    target_include_directories(rbdl_interface SYSTEM PRIVATE ${Eigen_INCLUDE_DIRS})
    target_include_directories(rbdl_interface SYSTEM PRIVATE ${rbdl_INCLUDE_DIRS})
    target_include_directories(rbdl_interface SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/contrib)
    target_include_directories(rbdl_interface PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(rbdl_interface PRIVATE ${rbdl_itf_LIBS})
    install(TARGETS rbdl_interface DESTINATION lib)
endif(rbdl_FOUND)


############################
# Build the main amp library
############################
set(amp_SOURCES
    amp/config.cpp
    amp/nodes.cpp
    amp/trajectories.cpp
    amp/interpolate.cpp
    amp/integrator.cpp
    amp/planner.cpp
    amp/rrt.cpp
    amp/spline_generator.cpp)

link_directories(${kdtree_LIBRARY_DIRS})
link_directories(${yaml-cpp_LIBRARY_DIRS})

set(amp_LIBS rbdl_interface kdtree)
add_library(amp SHARED ${amp_SOURCES})

if(yaml_FOUND OR yaml-cpp_FOUND)
    target_include_directories(amp SYSTEM PRIVATE ${yaml-cpp_INCLUDE_DIRS})
    list(APPEND amp_LIBS yaml-cpp) # -Wl,-Bdynamic
endif(yaml_FOUND OR yaml-cpp_FOUND)

if(bullet_FOUND)
    target_include_directories(amp SYSTEM PRIVATE ${bullet_INCLUDE_DIRS})
    link_directories(${bullet_LIBRARY_DIRS})
    list(APPEND amp_LIBS BulletCollision LinearMath)
endif(bullet_FOUND)

if(nlopt_FOUND)
    target_include_directories(amp SYSTEM PRIVATE ${nlopt_INCLUDE_DIRS})
    link_directories(${nlopt_LIBRARY_DIRS})
    list(APPEND amp_LIBS nlopt m)
endif(nlopt_FOUND)

target_include_directories(amp SYSTEM PRIVATE ${Eigen_INCLUDE_DIRS})
target_include_directories(amp SYSTEM PRIVATE ${kdtree_INCLUDE_DIRS})
target_include_directories(amp SYSTEM PRIVATE ${boost_INCLUDE_DIRS})
target_include_directories(amp SYSTEM PRIVATE ${rbdl_INCLUDE_DIRS})
target_include_directories(amp SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/contrib)
target_include_directories(amp PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(amp PRIVATE ${amp_LIBS})
install(TARGETS amp DESTINATION lib)
