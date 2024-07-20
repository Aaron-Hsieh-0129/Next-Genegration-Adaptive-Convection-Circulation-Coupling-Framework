cmake_minimum_required(VERSION 3.5.0)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(TMIF VERSION 0.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# set root of location to find PETSc's pkg-config
set(PETSC $ENV{PETSC_DIR})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

# Remove the lines below if you do not wish to have PETSc determine the compilers
# execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
# SET(CMAKE_C_COMPILER ${C_COMPILER})
# execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
# if (CXX_COMPILER)
#     SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
# endif (CXX_COMPILER)
# execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
# if (FORTRAN_COMPILER)
#     SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
#     enable_language(Fortran)
# endif (FORTRAN_COMPILER)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

set(CMAKE_EXPORT_COMPILE_COMMANDS True)


include_directories(
    ./2DVVM/include
    /home/Aaron/local/include
    /home/Aaron/local/petsc/include
)

file(GLOB_RECURSE SRC CSSWM/src/*.cpp CSSWM/src/*.hpp 2DVVM/src/*.cpp 2DVVM/src/*.hpp main.cpp)
# delete the main.cpp files from the 2DVVM and CSSWM directories
list(REMOVE_ITEM SRC ${CMAKE_CURRENT_SOURCE_DIR}/2DVVM/src/main.cpp)
list(REMOVE_ITEM SRC ${CMAKE_CURRENT_SOURCE_DIR}/CSSWM/src/main.cpp)

add_executable(TMIF ${SRC})

target_link_libraries(TMIF PkgConfig::PETSC)

find_library(libncPath netcdf "/home/Aaron/local/lib")
find_library(libncxxPath netcdf_c++4 "/home/Aaron/local/lib")
find_library(libpetscPath petsc "/home/Aaron/local/petsc/lib")
target_link_libraries(TMIF ${libncPath} ${libncxxPath} ${petsc})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
