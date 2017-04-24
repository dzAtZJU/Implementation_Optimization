cmake_minimum_required(VERSION 3.7)
project(Implementation_Optimization)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES GaussNewtonSolver/ExactLineSearch.cpp GaussNewtonSolver/ExactLineSearch.h GaussNewtonSolver/LineSearch.cpp GaussNewtonSolver/LineSearch.h GaussNewtonSolver/SoftLineSearch.cpp GaussNewtonSolver/SoftLineSearch.h GaussNewtonSolver/hw3_gn.h GaussNewtonSolver/EllipsoidFitting.cpp GaussNewtonSolver/EllipsoidFitting.h GaussNewtonSolver/Solver2776.cpp GaussNewtonSolver/Solver2776.h GaussNewtonSolver/MatrixUtility.cpp GaussNewtonSolver/MatrixUtility.h main.cpp)

add_executable(Implementation_Optimization ${SOURCE_FILES})
target_link_libraries(Implementation_Optimization ${OpenCV_LIBS} )
