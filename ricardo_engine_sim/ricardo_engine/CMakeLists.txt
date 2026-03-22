cmake_minimum_required(VERSION 3.16)
project(RicardoEngineSim VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -fsanitize=address")

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Find optional dependencies
find_package(OpenMP QUIET)
find_package(Threads REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)

# Collect all source files
file(GLOB_RECURSE SOURCES 
    src/mesh/*.cpp
    src/physics/*.cpp
    src/chemistry/*.cpp
    src/combustion/*.cpp
    src/engine/*.cpp
    src/analysis/*.cpp
    src/utils/*.cpp
)

# Main library
add_library(ricardo_engine_lib STATIC ${SOURCES})
target_include_directories(ricardo_engine_lib PUBLIC include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(ricardo_engine_lib PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(ricardo_engine_lib PUBLIC USE_OPENMP)
endif()
target_link_libraries(ricardo_engine_lib PUBLIC Threads::Threads m)

# Main executable
add_executable(ricardo_sim src/main.cpp)
target_link_libraries(ricardo_sim PRIVATE ricardo_engine_lib)

# Analysis executable
add_executable(combustion_analysis src/combustion_analysis_main.cpp)
target_link_libraries(combustion_analysis PRIVATE ricardo_engine_lib)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
