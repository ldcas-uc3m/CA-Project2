cmake_minimum_required(VERSION 3.10)
project(sim-aos)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -DCMAKE_BUILD_TYPE=Release -fopenmp)
add_executable(sim-paos sim-aos.cpp)
target_link_libraries (sim-paos PUBLIC OpenMP::OpenMP_CXX)
