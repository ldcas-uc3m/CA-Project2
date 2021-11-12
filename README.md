# CA-Project2
By team 89-02: Luis Daniel Casais Mezquida, Ivan Darío Cersósimo & Pablo Ruiz Fernández

## Project description
This project has as a fundamental goal to help students to acquire concerns on performance of parallel programs and to get familiarized with parallel programming techniques. In addition, they will get introduced to performance evaluation techniques for parallel programs.  
More specifically, the project focuses in the development of parallel software in the C++ programming language (including the latest versions), starting from sequential software previously developed.  
  
In this project you will generate parallel versions from the sequential versions of the gravitational simulation application developed in the [first project](https://github.com/ldcas-uc3m/CA-Project1). For code parallelization you will be using OpenMP.  
Two versions of the program will be implemented using the arrays of structures/objects (``paos``) and structures/objects of arrays (``psoa``) techniques.  
The program must accurately produce the same results than the sequential versions of the corresponding applications.

## Parallel version development
This task consists in the development of the parallel version of the described application in C++17 (or C++20).  
  
All your source files must compile without problems and shall not emit any compiler warning. In particular, the following compiler warning specific flags must be enabled:
```
-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors
```
Keep also in mind that you will have to perform all evaluations with compiler optimizations enabled (compiler options ``-O3`` and ``-DNDEBUG``). You can easily get this with ``-DCMAKE_BUILD_TYPE=Release``.  
You are allowed to use additional compiler flags as long as you document them in the project report and justify its use. Such flags must be properly included in the CMake configuration file.  

### Libraries
You are allowed to use any element from the C++ standard library that is included in your compiler distribution.  
You are also allowed to make use of the OpenMP library (header file ``<omp.h>``).  
Additional libraries are not allowed.  

## Performance evaluation
This task consists of the performance evaluation of the parallel application in its two versions (``paos`` and ``psoa``).  
To carry out the performance evaluation you must measure the application execution time. You are expected to represent graphically your results. Keep in mind the following considerations:
- All evaluations must be performed in the same machine that was used for the [first project](https://github.com/ldcas-uc3m/CA-Project1) of this course.
- You must include in your report all the relevant parameters of the machine where you have run the experiments (processor model, number of physical cores, main memory size, cache memory hierarchy, ...) as well as system software (operating system version, compiler version, ...).
    -  The used computer must have, as a minimum, four physical cores.
    - The operating system shall be GNU/Linux natively installed (using a virtual machine is not allowed).
    - The compiler shall be ``gcc`` version 9 or higher.
- Run each experiment a given number of times and take the average value. You are recommended to perform each experiment at least 10 or more executions so that you can give a confidence interval.
- Study results for different object populations. Consider cases with 4000 and 8000.
- Study results for different number of iterations: Consider cases with 250 and 500.

Represent graphically total execution times. Represent graphically average time per iteration.  
You must perform a performance evaluation considering different numbers of threads, from 1 to 16 (1, 2, 4, 8 and 16). You are expected to represent the obtained speedup over the sequential version.  
Include in your report the conclusions you may infer from results. Do not limit to describing data. You must search for a convincing explanation of results that includes the impact of the computer architecture.

