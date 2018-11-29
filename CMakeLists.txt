cmake_minimum_required(VERSION 3.6)
project(REO)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE FILEPATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE FILEPATH "C++ compiler" FORCE)

#FIND_PACKAGE(OpenMP REQUIRED)
#if(OPENMP_FOUND)
#    message("OPENMP FOUND")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#    include_directories(${OpenMP_INCLUDE_DIRS})
#endif()

set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
SET(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
set(OpenMP_EXE_LINKER_LFLAGS "-lomp")


set(SOURCE_FILES main.cpp Maths/Algebra.h Common/Solver.cpp Common/Solver.h Maths/Algebra.cpp REO/REOSolver.cpp REO/REOSolver.h Common/SpaceMesh.h Common/SpaceMesh.cpp ExactRiemannSolver/ExactRiemannSolver.h HLLC/HLLCSolver.cpp HLLC/HLLCSolver.h WENO/WENOSolver.h)
add_executable(REO ${SOURCE_FILES})
target_link_libraries(REO ${OpenMP_omp_LIBRARY})

