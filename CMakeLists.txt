cmake_minimum_required(VERSION 3.16)
project(Hermite-TENO)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -ffloat-store -O3")

find_package (Eigen3 3.3 REQUIRED)
include_directories(${Eigen_INCLUDE_DIR})
find_package(OpenMP)

add_executable(Strawberry main.cpp Equations_Euler1D.cpp Equations_Euler1D.h
        Equations_Euler2D.cpp Equations_Euler2D.h
        Equations_Scalar1D_Burgers.cpp Equations_Scalar1D_Burgers.h
        Equations_Scalar1D_LinearAdvection.cpp Equations_Scalar1D_LinearAdvection.h
        Equations_Scalar2D_LinearAdvection.cpp Equations_Scalar2D_LinearAdvection.h
        FV_Scalar1D.cpp FV_Scalar1D.h
        FV_Scalar2D.cpp FV_Scalar2D.h
        FV_System1D.cpp FV_System1D.h
        FV_System2D.cpp FV_System2D.h
        FVHermite_Scalar1D.cpp FVHermite_Scalar1D.h
        FVHermite_Scalar2D.cpp FVHermite_Scalar2D.h
        FVHermite_System1D.cpp FVHermite_System1D.h
        FVHermite_System2D.cpp FVHermite_System2D.h
        Schemes_Alternative.cpp Schemes_Alternative.h
        Schemes_Hermite.cpp Schemes_Hermite.h
        Schemes_Lagrange.cpp Schemes_Lagrange.h
        Schemes_SSPRK.cpp Schemes_SSPRK.h)
target_link_libraries(Strawberry OpenMP::OpenMP_CXX Eigen3::Eigen)