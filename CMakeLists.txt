# Set minimum CMake version
cmake_minimum_required(VERSION 3.10)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()


# Allow user to specify NVHPC and GCC installation paths
set(NVHPC_HOME "/home/Aaron/nvhpc/Linux_x86_64/21.3" CACHE STRING "Path to NVHPC installation")
set(GCC_HOME "/home/Aaron/gcc13" CACHE STRING "Path to GCC installation")

# Set compilers explicitly before project()
# set(CMAKE_C_COMPILER "${NVHPC_HOME}/comm_libs/mpi/bin/mpicc" CACHE PATH "C compiler" FORCE)
# set(CMAKE_CXX_COMPILER "${NVHPC_HOME}/comm_libs/mpi/bin/mpicxx" CACHE PATH "C++ compiler" FORCE)
set(CMAKE_C_COMPILER "/home/Aaron/openmpi/bin/mpicc" CACHE PATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER "/home/Aaron/openmpi/bin/mpicxx" CACHE PATH "C++ compiler" FORCE)
set(CMAKE_Fortran_COMPILER "${GCC_HOME}/bin/gfortran" CACHE PATH "Fortran compiler" FORCE)
set(CMAKE_CUDA_COMPILER "${NVHPC_HOME}/compilers/bin/nvcc" CACHE PATH "CUDA compiler" FORCE)


# Disable dependency file generation (-MF issue) for NVHPC
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    set(CMAKE_CXX_DEPENDS "")
    set(CMAKE_C_DEPENDS "")
    set(CMAKE_DEPFILE_FLAGS_CXX "")
    set(CMAKE_DEPFILE_FLAGS_C "")
endif()

# Define project after compiler settings
project(VVM2D VERSION 0.1.0 LANGUAGES CXX Fortran)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set HOME_DIR
set(HOME_DIR "/home/Aaron")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Option to enable or disable PETSc
option(USE_PETSC "Enable PETSc support" OFF)

if(USE_PETSC)
    # Set root location to find PETSc's pkg-config
    set(PETSC $ENV{PETSC_DIR})
    set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

    # Find PETSc using pkg-config
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
endif()

# Remove OpenMP since we're using MPI
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Set Fortran compiler flags
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fdefault-real-8")

# Enable exporting compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/2DVVM/include
    ${HOME_DIR}/local/include
    ${MPI_INCLUDE_PATH}
    ${CMAKE_SOURCE_DIR}/2DVVM/external/rte-rrtmgp-cpp/include
    ${CMAKE_SOURCE_DIR}/2DVVM/external/rte-rrtmgp-cpp/include_test
)

# Add PETSc include directories only if PETSc is enabled
if(USE_PETSC)
    include_directories(${HOME_DIR}/local/petsc/include)
endif()

# Source files (C++ files)
file(GLOB_RECURSE SRC_CXX CSSWM/src/*.cpp CSSWM/src/*.hpp 2DVVM/src/*.cpp 2DVVM/src/*.hpp src/*.cpp src/*.hpp)

# Delete the main.cpp files from the 2DVVM and CSSWM directories
list(REMOVE_ITEM SRC_CXX ${CMAKE_CURRENT_SOURCE_DIR}/2DVVM/src/main.cpp)
list(REMOVE_ITEM SRC_CXX ${CMAKE_CURRENT_SOURCE_DIR}/CSSWM/src/main.cpp)

# Fortran source files
file(GLOB_RECURSE SRC_Fortran 2DVVM/external/P3-microphysics/*.f90)

# Create a shared library for Fortran code
if(SRC_Fortran)
    add_library(fortran_lib SHARED ${SRC_Fortran})
    set_target_properties(fortran_lib PROPERTIES OUTPUT_NAME "fortran_lib")
endif()

# Add rte-rrtmgp-cpp as a subdirectory
add_subdirectory(2DVVM/external/rte-rrtmgp-cpp)

# Force rte-rrtmgp to use only _deps/install headers
set_target_properties(rte_rrtmgp PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/2DVVM/include;${CMAKE_SOURCE_DIR}/2DVVM/external/rte-rrtmgp-cpp/include"
)

# Add the executable target
add_executable(TMIF ${SRC_CXX})

# Link Fortran library and runtime to the executable
if(SRC_Fortran)
    target_link_libraries(TMIF PRIVATE fortran_lib gfortran)
    add_dependencies(TMIF fortran_lib)
endif()

# Link MPI libraries
target_link_libraries(TMIF PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# Link PETSc libraries if PETSc is enabled
if(USE_PETSC)
    target_link_libraries(TMIF PRIVATE PkgConfig::PETSC)
endif()

# Find and link other libraries
find_library(libncPath netcdf "${HOME_DIR}/local/lib")
find_library(libncxxPath netcdf_c++4 "${HOME_DIR}/local/lib")
find_library(libpetscPath petsc "${HOME_DIR}/local/petsc/lib")

target_link_libraries(TMIF PRIVATE ${libncPath} ${libncxxPath} rte_rrtmgp)

# Link PETSc libraries conditionally
if(USE_PETSC)
    target_link_libraries(TMIF PRIVATE ${libncPath} ${libncxxPath} ${libpetscPath} )
endif()

# Add MPI compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

# CPack settings
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
