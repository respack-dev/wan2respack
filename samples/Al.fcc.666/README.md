# Al.fcc.666 Sample

This sample demonstrates the calculation of Wannier functions and Coulomb interactions for aluminum in FCC structure.

## Directory Structure
```
Al.fcc.666/
├── inputs/
│   ├── wan2respack/
│   │   └── conf.toml
│   ├── QE/
│   │   ├── Al.scf.in
│   │   ├── Al.nscf.in
│   │   └── Al.pw2wan.in
│   ├── Wannier90/
│   │   └── Al.win.ref
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
$QE/bin/pw.x < QE/Al.scf.in > Al.scf.out
$QE/bin/pw.x < QE/Al.nscf.in > Al.nscf.out
```

3. Prepare Wannier90 input files:
```bash
python $PATH_to_Install/bin/wan2respack.py -pp wan2respack/conf.toml
```

4. Generate Wannier functions:
```bash
$QE/bin/pw.x < Al.nscf_wannier.in > Al.nscf_wannier.out
$Wannier90/wannier90.x -pp Al
$QE/bin/pw2wannier90.x < QE/Al.pw2wan.in > Al.pw2wan.out
$Wannier90/wannier90.x Al
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
- `Al.scf.in`: Self-consistent field calculation input
- `Al.nscf.in`: Non-self-consistent field calculation input
- `Al.pw2wan.in`: Input for pw2wannier90.x

### Wannier90 Input
- `Al.win.ref`: Reference input file for Wannier90

### RESPACK Input
- `respack.in`: Input file for RESPACK calculations

### Configuration
- `conf.toml`: Configuration file for wan2respack

## Expected Output
- `dir-wfn/`: Directory containing wave functions in RESPACK format
- `dir-wan/`: Directory containing Wannier functions and additional information
- `LOG.*`: Log files for RESPACK calculations 