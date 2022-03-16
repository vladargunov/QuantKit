cmake_minimum_required(VERSION 3.20)

project(QuantKit VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(${PROJECT_NAME} main.cpp Regression.cpp Regression.h FileManagement.cpp FileManagement.h Option.cpp Option.h nonuniform_grid.cpp nonuniform_grid.h Simulation.cpp Simulation.h HestonOption.cpp HestonOption.h OptimalExecution.cpp OptimalExecution.h)

# Add Eigen library
SET( EIGEN3_INCLUDE_DIR ../external/eigen)
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ("${EIGEN3_INCLUDE_DIR}")
