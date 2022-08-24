#
# This is a modified/original file distributed in RESPACK code under GNU GPL ver.3.
# https://sites.google.com/view/kazuma7k6r
#

# for Intel Compiler
# to detect omp flag
set(CMAKE_C_COMPILER "icc" CACHE STRING "" FORCE) 
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -qopenmp -xHost -g -traceback" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -qopenmp -g -traceback" CACHE STRING "" FORCE)

# for Intel MKL
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "" FORCE)
