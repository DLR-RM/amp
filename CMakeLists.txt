cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(amp)
set(amp_VERSION_MAJOR 1)
set(amp_VERSION_MINOR 0)
set(amp_VERSION_PATCH 0)
set(amp_VERSION	${amp_VERSION_MAJOR}.${amp_VERSION_MINOR}.${amp_VERSION_PATCH})

option(BUILD_RBDL_INTERFACE "build shared library interface to RBDL." ON)
option(BUILD_KDL_INTERFACE "build shared library interface to KDL." ON)
option(BUILD_TESTS "build unit tests." OFF) 		
option(BUILD_EXAMPLES "build module examples." ON) 		
option(BUILD_KDTREE "build the kdtree dependency if not already exists on system." OFF)
option(DEBUG "activate debug-only code sections - this option impacts comp time." ON) 
option(VERBOSE "activate verbose output code sections - this option impacts comp time." ON)
option(LOG_ALL_DATA "activate agressive data logging - this option impacts comp time." OFF)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
		message(FATAL_ERROR "Insufficient GCC version!")
    endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.7)
		message(FATAL_ERROR "Insufficient Clang version!")
    endif()
endif()

########################
# Required Dependencies
########################
find_package(PkgConfig)

find_package(Eigen 3.3.2)
if(NOT Eigen_FOUND) #if no .cmake file try pkg-config
    pkg_search_module(Eigen REQUIRED eigen3>=3.3.2)
endif(NOT Eigen_FOUND)

find_package(boost 1.58.0)
if(NOT boost_FOUND) #if no .cmake file try pkg-config
    pkg_search_module(boost REQUIRED boost>=1.58.0)
endif(NOT boost_FOUND)

if(BUILD_KDTREE)
    #include(ExternalProject)
    #ExternalProject_Add(kdtree
    #    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/contrib/kdtree
    #    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/contrib/kdtree
    #    CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/contrib/kdtree/configure --prefix=share
    #    BUILD_COMMAND ${MAKE}})
    #ExternalProject_Get_Property(lib_kdtree install_dir)
    #add_library(kdtree SHARED IMPORTED)
    #set_property(TARGET kdtree PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libkdtree.so)
    #add_dependencies(kdtree project_kdtree)
else(BUILD_KDTREE)
    find_package(kdtree)
    if(NOT kdtree_FOUND) #if no .cmake file try pkg-config
        pkg_search_module(kdtree REQUIRED kdtree)
    endif(NOT kdtree_FOUND)
endif(BUILD_KDTREE)

#######################
# Optional Dependencies
#######################
pkg_search_module(rbdl REQUIRED rbdl>=2.5.0)
find_package(yaml-cpp 0.5.3)
#find_package(nlopt 2.4.2 REQUIRED nlopt)
#find_package(bullet LinearMath BulletCollision)
find_package(tinyxml 2.6.2)

##################
# compiler options
##################
set(CMAKE_CXX_FLAGS "-Wall -fPIC -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb") # -DEIGEN_DONT_VECTORIZE")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG ")

if(DEBUG)
    add_definitions(-DDEBUG)
endif(DEBUG)

if(VERBOSE)
    add_definitions(-DVERBOSE)
endif(VERBOSE)

if(LOG_ALL_DATA)
	add_definitions(-DLOG_ALL_DATA)
endif(LOG_ALL_DATA)

#############################
# Build Libraries and targets
#############################
add_subdirectory(src)

if(BUILD_TESTS)
	add_subdirectory(tests)
endif(BUILD_TESTS)

if(BUILD_EXAMPLES)
	add_subdirectory(examples)
endif(BUILD_EXAMPLES)

