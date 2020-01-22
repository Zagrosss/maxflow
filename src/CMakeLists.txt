set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE "RELEASE")
#set(CMAKE_BUILD_TYPE "DEBUG")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CXX_COVERAGE_COMPILE_FLAGS "-Wall -pedantic -Wextra -g -O3 -march=native -funroll-loops -ftree-vectorize -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_COVERAGE_COMPILE_FLAGS}")

add_executable(maxflow maxflow.cpp)
target_include_directories(maxflow PRIVATE lib)
