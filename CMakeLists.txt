# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# Start configuration
cmake_minimum_required(VERSION 2.8)
project(HDF5EXTRACT CXX)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# Linking and include info
link_libraries(boost_program_options hdf5 hdf5_cpp ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL})

# Set the standard so that the compiler doesnt complain
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

# Create executable
add_executable(h5extract h5extract.cpp)
triqs_set_rpath_for_target(h5extract)
