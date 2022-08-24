#
# This is a modified/original file distributed in RESPACK code under GNU GPL ver.3.
# https://sites.google.com/view/kazuma7k6r
#

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(WAN2RESPACK NONE)

if(CONFIG)    
  message(STATUS "Loading configration: " ${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
  include(${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
endif(CONFIG)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
enable_language(Fortran C)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

option(DEBUG "Debug" OFF)

find_package(LAPACK)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
#CMAKE_Fortran_COMPILER_VERSION cannot be detected
if("${CMAKE_C_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
     set(OMP_FLAG_Intel "-openmp")
  else()
     set(OMP_FLAG_Intel "-qopenmp")
  endif()
  find_package(OpenMP)
  if(OPENMP_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OMP_FLAG_Intel}")
  endif(OPENMP_FOUND)
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
  endif(OPENMP_FOUND)
endif()

add_subdirectory(util/wan2respack)

