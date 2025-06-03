# SrVO3.sc.666 Sample

This sample demonstrates the calculation of Wannier functions and Coulomb interactions for SrVO3 in simple cubic structure.

## Directory Structure
```
SrVO3.sc.666/
├── inputs/
│   ├── wan2respack/
│   │   └── conf.toml
│   ├── QE/
│   │   ├── SrVO3.scf.in
│   │   ├── SrVO3.nscf.in
│   │   └── SrVO3.pw2wan.in
│   ├── Wannier90/
│   │   └── SrVO3.win.ref
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
$QE/bin/pw.x < QE/SrVO3.scf.in > SrVO3.scf.out
$QE/bin/pw.x < QE/SrVO3.nscf.in > SrVO3.nscf.out
```

3. Prepare Wannier90 input files:
```bash
python $PATH_to_Install/bin/wan2respack.py -pp wan2respack/conf.toml
```

4. Generate Wannier functions:
```bash
$QE/bin/pw.x < SrVO3.nscf_wannier.in > SrVO3.nscf_wannier.out
$Wannier90/wannier90.x -pp SrVO3
$QE/bin/pw2wannier90.x < QE/SrVO3.pw2wan.in > SrVO3.pw2wan.out
$Wannier90/wannier90.x SrVO3
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
- `SrVO3.scf.in`: Self-consistent field calculation input
- `SrVO3.nscf.in`: Non-self-consistent field calculation input
- `SrVO3.pw2wan.in`: Input for pw2wannier90.x

### Wannier90 Input
- `SrVO3.win.ref`: Reference input file for Wannier90

### RESPACK Input
- `respack.in`: Input file for RESPACK calculations

### Configuration
- `conf.toml`: Configuration file for wan2respack

## Expected Output
- `dir-wfn/`: Directory containing wave functions in RESPACK format
- `dir-wan/`: Directory containing Wannier functions and additional information
- `LOG.*`: Log files for RESPACK calculations

## Notes
- This calculation focuses on the V 3d orbitals, which are the main contributors to the electronic properties near the Fermi level.
- The calculation is intermediate in complexity between Al.fcc.666 and La2CuO4.bct.666 examples. 