********
Tutorial
********

In this tutorial, we demonstrate the following calculations:
(i) performing the first principle calculations using *QE*,
(ii) calculating Wannier functions using *Wannier90*,
(iii) generating input files for *RESPACK* using *wan2respack*, and
(iv) calculating Coulomb interactions using *RESPACK*.
The sample files are located in ``samples/Al.fcc.666``.
In this directory, there are the following four directories:

1. PP

  - Including files of pseudopotentials.

2. inputs

  - Including input files.

3. inputs_selfk

  - Including input files when setting k-points manually by the user.

4. reference

  - Including reference input files for generating Wannier functions using ``calc_wannier`` in *RESPACK* .

Al.fcc
========

The sample files of this tutorial are located in ``samples/Al.fcc.666/inputs`` .
First, change the directory to ``samples/Al.fcc.666/inputs``::

    $cd samples/Al.fcc.666/inputs

In this directory, the following input files are included:

- QE

  - Al.scf.in: Input file for scf calculation.
  - Al.nscf.in: Input file for nscf calculation.
  - Al.pw2wan.in: Input file for generating ``mmn`` and ``amn`` files.

- Wannier90

  - Al.win.ref: Reference file of *wan2respack* for generating an input file of *Wannier90* .

- wan2respack

  - conf.toml: Input file of *wan2respack* for setting a path to the output directory of *QE* and seed names, etc.

- RESPACK

  - respack.in : Input file *RESPACK* .

First principle calculations for the irreducible k-points using *QE*.
-----------------------------------------------------------------------------

By typing the following commands, first principles calculations of scf and nscf using *QE* will be performed::

    $QE/bin/pw.x < Al.scf.in > Al.scf.out
    $QE/bin/pw.x < Al.nscf.in > Al.nscf.out

Here, ``$QE`` indicates a path to a directory where *QE* is installed.
**Note that target k-points must be irreducible.**

Export k-points to be calculated by *Wannier90*.
-----------------------------------------------------------------------------------------

Next, as a preprocess, export k-points to be calculated to input files of nscf calculation and Wannier90 by typing the following command: ::

    $python $PATH_to_Install/bin/wan2respack.py -pp conf.toml

The contents of ``conf.toml`` are shown below: ::

    [base]
    QE_output_dir = "./work/Al.save"
    seedname = "Al"

    [pre.ref]
    nscf = "Al.nscf.in"
    win = "Al.win.ref"

    [pre.output]
    nscf = "Al.nscf_wannier.in"
    win = "Al.win"

In the ``[base]`` section, an output directory of *QE* and seed name are indicated by ``QE_output_dir`` and a ``seedname`` , respectively.
In the ``[pre.ref]`` section, reference files for generating input files are indicated by ``nscf`` and ``win``, respectively.
In the ``[pre.output]`` section, names of output files generated by *wan2respack* are indicated by ``nscf`` and ``win``, respectively.

After finishing a calculation, ``dir-wfn`` directory, ``Al.nscf_wannier.in`` and ``Al.win`` files will be generated
(k-points will be  added to ``Al.nscf_wannier.in`` and ``Al.win`` ).

Generate Wannier functions with *QE* and *Wannier90*
-----------------------------------------------------------------------------------------

Using ``Al.nscf_wannier.in`` and ``Al.win``,
Wannier functions are generated using *QE* and *Wannier90* by typing the following commands: ::

    $QE/bin/pw.x < Al.nscf_wannier.in > Al.nscf_wannier.out
    $Wanier90/wannier90.x -pp Al
    $QE/bin/pw2wannier90.x < Al.pw2wan.in > Al.pw2wan.out
    $Wanier90/wannier90.x Al

Convert Wannier functions to *RESPACK* format
-----------------------------------------------------------------------------------------

Converting Wannier functions to *RESPACK* format can be performed using ``wan2respack.py``
The execution command is described below: ::

    $python $PATH_to_Install/bin/wan2respack.py conf.toml

After finishing calculations,  four files are generated in the ``dir-wan`` directory.

Calculation of Coulomb interactions using *RESPACK*
-----------------------------------------------------------------------------------------

The input file of *RESPACK* is prepared as ``respack.in`` .
Using this file, we can calculate Coulomb interactions using constrained Random Phase Approximation by *RESPACK* .

The execution command is described below. ::

    $RESPACK/bin/calc_chiqw < respack.in > LOG.chiqw
    $RESPACK/bin/calc_w3d < respack.in > LOG.W3d
    $RESPACK/bin/calc_j3d < respack.in > LOG.J3d

The obtained results are shown in :numref:`WvsR`.
The horizontal axis indicates the distance, and the vertical axis indicates the screened Coulomb interaction.
In these figures, we also plotted the numerical results obtained using *RESPACK*.
In this case, the calculation of Wannier functions was performed by ``calc_wannier`` in *RESPACK*.
The input file ``respack.in`` of ``calc_wannier`` is located in the ``reference`` directory.
Qualitatively, almost the same trend is obtained, indicating that the tool is working well.
The difference shown in the inset of :numref:`WvsR` is due to the fact
that the Wannier function obtained with Wannier90 is not maximally localized (num_iter =0).

.. figure:: ../../figs/wvsr-Al-1.pdf
    :name: WvsR
    :scale: 40%

    Coulomb interactions obtained by constrained Random Phase Approximation.
    W90-Wannier indicates the numerical results obtained by this tutorial.
    RESAPCK-Wannier indicates the numerical results obtained using *RESPACK* (the calculation of Wannier functions was performed by ``calc_wannier`` in *RESPACK* ).
    The inset shows an enlarged view of Coulomb interactions obtained by constrained Random Phase Approximation.

We also prepare sample files of La2CuO4 and SrVO3 in the ``samples`` directory.

[Optional] Set k-points manually by the user
--------------------------------------------

In the above tutorial, k-points are automatically exported.
Here, k-points can be set manually by the user, if desired, by typing the following command in the ``inputs_selfk`` directory: ::

    $python $PATH_to_Install/bin/wan2respack.py -pp conf.toml

The contents of ``conf.toml`` are shown below: ::

    [base]
    QE_output_dir = "./work/Al.save"
    seedname = "Al"
    selfk = true

The ``selfk`` flag in the ``[base]`` section must be ``true`` in this mode.
The k-point list is written in ``dat.sample_mk``.
A k-point list in ``Al.nscf_wannier.in`` and ``Al.win`` is determined based on ``dat.sample_mk``.

