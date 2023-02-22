******************
File specification
******************

.. _file_spec:

Input files
===========
This section explains all input files from SCF calculation to Coulomb interaction calculation.

- QE scf input file

- QE nscf input file

- Wannier90 input file
    Do not write a k-point block.
    This file is used only as a reference.

- pw2wannier90 input file

- RESPACK input file

- conf.toml 
    The format is shown below. ::

        [base]
        QE_output_dir = "./work/prefix.save"
        seedname = "seedname"
        selfk = false            #(optional) Default: false

        [pre.ref]
        nscf = "nscf.in"
        win = "seedname.win.ref"

        [pre.output]
        nscf = "nscf_wannier.in"
        win = "seedname.win"


    1. ``base``
        - ``QE_output_dir``: *QE* output directory
        - ``seedname``: Same string as the seedname used in *Wannier90*
        - ``selfk`` (optional, Default: ``false``): Flag to set k-points manually at the pre-process mode

    2. ``pre.ref``
        - ``nscf``: File name of the *QE* nscf input file prepared by the user
        - ``win``: File name of the *Wannier90* input file prepared by the user

    3. ``pre.output``
        - ``nscf``: File name of the new *QE* nscf input file that is automatically generated based on ``[pre.ref]nscf``
        - ``win``: File name of the new *Wannier90* input that is automatically generated based on ``[pre.ref]win``

Output files
============
The details of the output files are explained.

Preprocess
-------------

- ``[pre.output]nscf``
    The *QE* input file having the name determined by ``[pre.output]nscf`` in conf.toml.
    This file is automatically made based on the reference file: ``[pre.ref]nscf``.

- ``[pre.output]win``
    The *Wannier90* input file having the name determined by ``[pre.output]win`` in conf.toml.
    This file is automatically made based on the reference file: ``[pre.ref]win``.

- dat.sample_mk
    The intermediate file for making input files of *QE* and *Wannier90*, including k-points.
    The first line gives the total number of k-points.
    The next block gives k-points in terms of the reciprocal lattice vectors.

- dat.kg_respack
    The intermediate file for making ``dat.wan``, including G-vectors.
    This file consists of the number of blocks equal to the total number of k-points.
    The first line of each block gives the number of G-vectors.
    The remaining lines of each block give the G-vectors in terms of the reciprocal lattice vectors.

- LOG.mk
    Log file.

Coreprocess
-------------

- dat.ns-nb, dat.umat, dat.wan, dat.wan-center
    These files include information related to Wannier functions.
    The file format is the same as that of *RESPACK*.
    See :ref:`file_expression-label` for details.

- LOG.genwan
    Log file.
