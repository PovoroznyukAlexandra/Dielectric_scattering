cmake_minimum_required(VERSION 3.20)
project(Dielectric_scattering)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapack")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(LAPACK REQUIRED)


add_executable(Dielectric_scattering dielectric_scattering.cpp diel_scat_func.h diel_scat_func.cpp)
