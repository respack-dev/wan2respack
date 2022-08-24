***********************************
How to use wan2respack
***********************************

Prerequisite
============
**Required packages**

Install the following three packages

1. `Quantum Espresso <https://www.quantum-espresso.org>`_ (6.6)

2. `Wannier90 <http://www.wannier.org>`_ (3.0.0)

3. `RESPACK <https://sites.google.com/view/kazuma7k6r>`_ (20200113)

The versions in parentheses have been tested.

**About Python**

Python 3.6 or higher version is required.

Python requires tomli library. 
This tomli library is the standard library for python 3.11.
Execute the following command. ::

  pip install tomli

Structure
=========

*wan2respack* has the following directory structure ::

  |--CMakeLists.txt
  |--LICENSE         
  |--README.md       
  |
  |--config
  |    |
  |    |--intel.cmake
  |    |--gcc.cmake
  |
  |--docs
  |
  |--samples
  |    |
  |    |--README.md
  |    |--Al.fcc.666
  |    |--La2CuO4.bct.666        
  |    |--SrVO3.sc.666
  |
  |--util
       |
       |--wan2respack

Installation
============
*wan2respack* can be downloaded from the following GitHub page.
https://github.com/respack-dev/wan2respack

Users can complie *wan2respack* by using CMake.
An example for the installtion of *wan2respack* is as follows ::
  cd $PATH_to_wan2respack
  mkdir build 
  cd build
  cmake ../ -DCONFIG=$Type_of_Configure -DCMAKE_INSTALL_PREFIX=$PATH_to_Install
  make
  make install

, where ``$PATH_to_wan2respack`` is the path to *wan2respack* directory and ``$PATH_to_Install`` is the path to the directory for the installation.
By replacing ``$Type_of_Configure`` for a name of CMake configure files, users can specify compilers they want to use.
In :math:`\alpha` version, the following ``$Type_of_Configure`` are available ::
  intel: Intel compiler + MKL
  gcc: GCC compiler

The details of compiler options can be found in CMake configure files in ``$PATH_to_wan2respack/config`` directory.

All the binary files and the Python script are installed to ``$PATH_to_Install/bin``.
Their details are as follows:

  ``wan2respack.py``
   - Main Python script including two modes: pre-process and core-process. 
     This script requires a configuration toml file as an argument.

  ``gen_mk.x``
   - Fortran90 code for calculating k-points mesh for Wannier90. This binary is called by wan2respack_pre.py.
  
  ``gen_wan.x``
   - Fortran90 code for converting the Wannier90 results into the RESPACK Wannier format. This binary is called by wan2respack_core.py.
  
  ``wan2respack_pre.py``
   - Python script for saving the QE results and exporting k-points with gen_mk.x and qe2respack.py.
  
  ``wan2respack_core.py``
   - Python script for preparing files about Wannier functions in RESPACK format with gen_wan.x and qe2respack.py.

  ``qe2respack.py``
   - Python script for generating input files of RESPACK from QE band calculations. This script is originally distributed under GNU GPL ver.3 by open-source software RESPACK ver. 20200113.
  
  ``init.py``
   - Python module in which common functions are defined.

Basic usage
===========

1. Perform first principles calculations by QE

  - Only norm-conserving pseudopotentials can work in RESPACK.
  - Perform the calculation at the irreducible k-points. 

2. Running wan2respack.py in pre-process mode

  - Export k-points to be calculated to nscf-input and wannier90-input.
  
3. Generate Wannier functions by QE & Wannier90

  - Use the input files made by the previous step.

4. Convert Wannier functions to RESPACK format by running wan2respack.py 


Prepare input files
-------------------
Prepare the following files. See :ref:`file_spec` for details.

- QE scf input.

- QE nscf input.

  -  Be sure to use *{automatic}* to set the k-point.

- wannier90 input file.

  - Do not write kpoints-block.

- pw2wannier90 input file.

- RESPACK input file.

- configuration toml file. 


Run *wan2respack*
----------------
After the calculations at the irreducible k-points, ::

  python wan2respack.py conf.toml -pp

The above command generates ``new_nscf`` and ``new_win`` files with the k-points list to be calculated.

After the Wannier functions are generated, ::

  python wan2respack.py conf.toml

The ``dir-wan`` directory and four files inside it are generated.
