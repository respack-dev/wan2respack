#!/bin/bash

QE_DIR=/path/to/qe-6.6/
WAN90_DIR=/path/to/wannier90-3.1.0/
RES=/path/to/RESPACK-20200113/
PATH_to_INSTALL=/path/to/wan2respack/build

# scf
mpirun -n 16 $QE_DIR/bin/pw.x -nk 4 < La2CuO4.scf.in > La2CuO4.scf.out

# nscf IBZ calc
mpirun -n 16 $QE_DIR/bin/pw.x -nk 4 < La2CuO4.nscf.in > La2CuO4.nscf.out

# pre-process
python $PATH_to_INSTALL/bin/wan2respack.py -pp conf.toml

# wannier90 run
mpirun -n 16 $QE_DIR/bin/pw.x -nk 4 < La2CuO4.nscf_wannier.in > La2CuO4.nscf_wannier.out
$WAN90_DIR/wannier90.x -pp La2CuO4
mpirun -n 16 $QE_DIR/bin/pw2wannier90.x < La2CuO4.pw2wan.in > La2CuO4.pw2wan.out
$WAN90_DIR/wannier90.x La2CuO4

# wannier90 results to RESPACK inputs
python $PATH_to_INSTALL/bin/wan2respack.py conf.toml

# RESPACK
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
mpirun -np 4 $RES/src/chiqw/calc_chiqw < respack.in > LOG.chiqw

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
$RES/src/calc_int/calc_w3d < respack.in > LOG.w3d
$RES/src/calc_int/calc_j3d < respack.in > LOG.j3d
