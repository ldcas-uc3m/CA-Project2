cmake_minimum_required(VERSION 3.10)
project(sim-aos)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -DCMAKE_BUILD_TYPE=Release)
add_executable(sim-paos sim-aos.cpp)
