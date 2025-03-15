cmake_minimum_required(VERSION 3.5.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(TMIF VERSION 1.0.0 LANGUAGES CXX Fortran)  # Added Fortran language

set(HOME_DIR $ENV{HOME})

# Set compilers
set(CMAKE_C_COMPILER "${HOME_DIR}/gcc13/bin/gcc")
set(CMAKE_CXX_COMPILER "${HOME_DIR}/gcc13/bin/g++")
set(CMAKE_Fortran_COMPILER "${HOME_DIR}/gcc13/bin/gfortran")  # Added Fortran compiler

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Option to enable or disable PETSc
option(USE_PETSC "Enable PETSc support" OFF)

if(USE_PETSC)
    # Set root location to find PETSc's pkg-config
    set(PETSC $ENV{PETSC_DIR})
    set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

    # Uncomment these lines if you want PETSc to determine the compilers
    execute_process(COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
    SET(CMAKE_C_COMPILER ${C_COMPILER})
    
    execute_process(COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
    if (CXX_COMPILER)
        SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
    endif()

    execute_process(COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
    if (FORTRAN_COMPILER)
        SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
    endif()

    # Find PETSc using pkg-config
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
endif()

# Find and enable OpenMP if available
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")  # Added for Fortran
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Set Fortran compiler flags (similar to successful case)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8")

# Enable exporting compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/2DVVM/include
    ${HOME_DIR}/local/include
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

# Fortran source files (adjust path as needed)
file(GLOB_RECURSE SRC_Fortran 2DVVM/src/*.f90)  # Add your Fortran files here

# Create a shared library for Fortran code (if applicable, adjust name and files)
if(SRC_Fortran)
    add_library(fortran_lib SHARED ${SRC_Fortran})
    set_target_properties(fortran_lib PROPERTIES OUTPUT_NAME "fortran_lib")
endif()

# Add the executable target
add_executable(TMIF ${SRC_CXX})

# Link Fortran library and runtime to the executable
if(SRC_Fortran)
    target_link_libraries(TMIF PRIVATE fortran_lib gfortran)  # Link Fortran runtime
    add_dependencies(TMIF fortran_lib)  # Ensure Fortran library is built first
endif()

# Link PETSc libraries if PETSc is enabled
if(USE_PETSC)
    target_link_libraries(TMIF PRIVATE PkgConfig::PETSC)
endif()

# Find and link other libraries
find_library(libncPath netcdf "${HOME_DIR}/local/lib")
find_library(libncxxPath netcdf_c++4 "${HOME_DIR}/local/lib")
find_library(libpetscPath petsc "${HOME_DIR}/local/petsc/lib")

target_link_libraries(TMIF PRIVATE ${libncPath} ${libncxxPath})

# Link PETSc libraries conditionally
if(USE_PETSC)
    target_link_libraries(TMIF PRIVATE ${libncPath} ${libncxxPath} ${libpetscPath})
endif()

# CPack settings
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
