# Sample Calculations

This directory contains example calculations demonstrating the use of wan2respack for different materials and structures. Each sample includes complete input files and step-by-step instructions for calculating Wannier functions and Coulomb interactions.

## Available Samples

### 1. Al.fcc.666
- **Structure**: Face-centered cubic (FCC) aluminum
- **Complexity**: Basic
- **Features**:
  - Simple metallic system
  - Small unit cell
  - Good starting point for learning the workflow
- **Key orbitals**: Al 3s and 3p orbitals
- **Recommended for**: Beginners and testing the installation

### 2. SrVO3.sc.666
- **Structure**: Simple cubic SrVO3
- **Complexity**: Intermediate
- **Features**:
  - Transition metal oxide
  - Correlated electron system
  - Intermediate unit cell size
- **Key orbitals**: V 3d orbitals
- **Recommended for**: Learning about correlated systems

### 3. La2CuO4.bct.666
- **Structure**: Body-centered tetragonal La2CuO4
- **Complexity**: Advanced
- **Features**:
  - High-temperature superconductor parent compound
  - Large unit cell
  - Complex electronic structure
- **Key orbitals**: Cu 3d and O 2p orbitals
- **Recommended for**: Advanced users and complex systems

## Directory Structure
Each sample directory contains:
```
sample_name/
├── inputs/
│   ├── wan2respack/
│   │   └── conf.toml
│   ├── QE/
│   │   ├── [material].scf.in
│   │   ├── [material].nscf.in
│   │   └── [material].pw2wan.in
│   ├── Wannier90/
│   │   └── [material].win.ref
│   ├── RESPACK/
│   │   └── respack.in
│   └── Script/
│       └── submit.sh
└── README.md
```

## Pseudopotentials
Each sample directory also contains a `PP` directory with pseudopotential files from the [SG15 library](http://www.quantum-simulation.org/potentials/sg15_oncv/). These pseudopotentials are used in the Quantum ESPRESSO calculations.

## How to Use
1. Choose a sample based on your needs and experience level
2. Navigate to the sample directory
3. Follow the instructions in the sample's README.md file
4. Each sample includes complete input files and step-by-step calculation procedures

## Computational Requirements
- **Al.fcc.666**: Minimal computational resources required
- **SrVO3.sc.666**: Moderate computational resources required
- **La2CuO4.bct.666**: Significant computational resources required

## Notes
- All samples use the same basic workflow but with different input parameters
- The complexity increases from Al.fcc.666 to La2CuO4.bct.666
- Each sample's README.md contains detailed instructions specific to that material
- Some input files are modified/original files distributed in RESPACK code under GNU GPL ver.3
  (https://sites.google.com/view/kazuma7k6r)
