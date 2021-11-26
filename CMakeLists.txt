cmake_minimum_required(VERSION 3.10)
project(sim-soa)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -DCMAKE_BUILD_TYPE=Release -fopenmp)
add_executable(sim-psoa sim-soa.cpp)
target_link_libraries (sim-psoa PUBLIC OpenMP::OpenMP_CXX)
