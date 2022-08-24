#!/bin/bash

QE_DIR=/path/to/qe-6.6/
WAN90_DIR=/path/to/wannier90-3.1.0/
RES=/path/to/RESPACK-20200113/


# prepare(make dir wfn)
mpirun -n 16 $QE_DIR/bin/pw.x -nk 4 < La2CuO4.scf.in > La2CuO4.scf.out
mpirun -n 16 $QE_DIR/bin/pw.x -nk 4 < La2CuO4.nscf.in > La2CuO4.nscf.out

python $RES/util/qe2respack/qe2respack.py ./work/La2CuO4.save/

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16

# wannierize
$RES/src/wannier/calc_wannier < input.in > LOG.wannier
$RES/src/chiqw/calc_chiqw < input.in > LOG.chiqw
$RES/src/calc_int/calc_w3d < input.in > LOG.w3d
$RES/src/calc_int/calc_j3d < input.in > LOG.j3d

