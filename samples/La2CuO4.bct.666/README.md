# La2CuO4.bct.666 Sample

This sample demonstrates the calculation of Wannier functions and Coulomb interactions for La2CuO4 in body-centered tetragonal structure.

## Directory Structure
```
La2CuO4.bct.666/
├── inputs/
│   ├── wan2respack/
│   │   └── conf.toml
│   ├── QE/
│   │   ├── La2CuO4.scf.in
│   │   ├── La2CuO4.nscf.in
│   │   └── La2CuO4.pw2wan.in
│   ├── Wannier90/
│   │   └── La2CuO4.win.ref
│   ├── RESPACK/
│   │   └── respack.in
│   └── Script/
│       └── submit.sh
└── README.md
```

## Calculation Steps

1. Navigate to the inputs directory:
```bash
cd inputs
```

2. Perform DFT calculations:
```bash
$QE/bin/pw.x < QE/La2CuO4.scf.in > La2CuO4.scf.out
$QE/bin/pw.x < QE/La2CuO4.nscf.in > La2CuO4.nscf.out
```

3. Prepare Wannier90 input files:
```bash
python $PATH_to_Install/bin/wan2respack.py -pp wan2respack/conf.toml
```

4. Generate Wannier functions:
```bash
$QE/bin/pw.x < La2CuO4.nscf_wannier.in > La2CuO4.nscf_wannier.out
$Wannier90/wannier90.x -pp La2CuO4
$QE/bin/pw2wannier90.x < QE/La2CuO4.pw2wan.in > La2CuO4.pw2wan.out
$Wannier90/wannier90.x La2CuO4
```

5. Convert to RESPACK format:
```bash
python $PATH_to_Install/bin/wan2respack.py wan2respack/conf.toml
```

6. Calculate Coulomb interactions:
```bash
$RESPACK/bin/calc_chiqw < RESPACK/respack.in > LOG.chiqw
$RESPACK/bin/calc_w3d < RESPACK/respack.in > LOG.W3d
$RESPACK/bin/calc_j3d < RESPACK/respack.in > LOG.J3d
```

## Input Files Description

### Quantum ESPRESSO Inputs
- `La2CuO4.scf.in`: Self-consistent field calculation input
- `La2CuO4.nscf.in`: Non-self-consistent field calculation input
- `La2CuO4.pw2wan.in`: Input for pw2wannier90.x

### Wannier90 Input
- `La2CuO4.win.ref`: Reference input file for Wannier90

### RESPACK Input
- `respack.in`: Input file for RESPACK calculations

### Configuration
- `conf.toml`: Configuration file for wan2respack

## Expected Output
- `dir-wfn/`: Directory containing wave functions in RESPACK format
- `dir-wan/`: Directory containing Wannier functions and additional information
- `LOG.*`: Log files for RESPACK calculations

## Notes
- This calculation requires more computational resources than the Al.fcc.666 example due to the larger unit cell and more complex electronic structure.
- The Wannier functions are centered on the Cu 3d orbitals and O 2p orbitals. 