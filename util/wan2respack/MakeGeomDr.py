#!/usr/bin/env python
"""
RESPACK interface code for Wannier90.

This script converts Wannier90 output files to RESPACK format, specifically creating
geometry and distance-related files (zvo_geom.dat and zvo_dr.dat) needed for RESPACK
calculations.

The script performs the following main tasks:
1. Reads Wannier90 output files (.wout and _hr.dat)
2. Extracts lattice vectors and Wannier function centers
3. Converts coordinates to fractional coordinates
4. Generates RESPACK input files
5. Calculates Green's functions using hwave
6. Creates RESPACK-compatible output files

Dependencies
-----------
- numpy: For numerical computations
- hwave: For Green's function calculations
  - hwave is a Python package for quantum many-body calculations
  - It must be installed and available in the Python environment
  - The hwave.qlms module is used for the UHFk calculation

Citation
--------
If you use hwave in your research, please cite:
    T. Aoyama, K. Yoshimi, K. Ido, Y. Motoyama, T. Kawamura, T. Misawa, T. Kato and A. Kobayashi,
    "H-wave – A Python package for the Hartree-Fock approximation and the random phase approximation",
    Computer Physics Communications, 2024, 298, 109087.
    https://doi.org/10.1016/j.cpc.2024.109087

Examples
--------
Basic usage:
    $ python MakeGeomDr.py conf.toml

With preprocessing flag (only generates zvo_geom.dat):
    $ python MakeGeomDr.py conf.toml --pp

Note
----
The script requires hwave to be properly installed and configured. If hwave is not
available, the script will fail when trying to calculate Green's functions.
"""

import argparse
import os
import numpy as np
try:
    import hwave.qlms
except ImportError:
    raise ImportError(
        "hwave package is required but not found. "
        "Please install hwave using pip or your package manager."
    )

def main():
    """
    Main function to execute the RESPACK interface code for Wannier90.

    This function parses command-line arguments, reads configuration from a
    TOML file, and generates necessary data files for further processing.
    When the preprocessing flag is not set, it also performs Green's function
    calculations using hwave.qlms.

    Examples
    --------
    Run the script with a configuration file:
        $ python MakeGeomDr.py conf.toml

    Run the script with preprocessing flag:
        $ python MakeGeomDr.py conf.toml --pp

    Notes
    -----
    The script requires hwave to be properly installed and configured. The hwave
    calculation is skipped if the preprocessing flag (-pp) is set.
    """
    parser = argparse.ArgumentParser(
        description="RESPACK interface code for Wannier90"
    )
    parser.add_argument(
        'conf', type=str, metavar="conf.toml", help='Configure TOML file.'
    )
    parser.add_argument(
        "-pp", "--pp", action="store_true", help="Flag for preprocessing"
    )
    args = parser.parse_args()
    conf = simple_variables(args)
    seedname = conf["base"]["seedname"]

    # [s] make zvo_geom.dat
    file_name = f"{seedname}.wout"
    vec_lat = GetLatticeVec(file_name)
    vec_wan = GetWanCentre(file_name)
    frac_wan = Convert2Frac(vec_lat, vec_wan)
    with open("zvo_geom.dat", 'w') as f:
        for i in range(3):
            print(f"  {vec_lat[i][0]} {vec_lat[i][1]} {vec_lat[i][2]} ", file=f)
        print(f"  {frac_wan.shape[0]} ", file=f)
        for cnt in range(frac_wan.shape[0]):
            print(f"  {frac_wan[cnt][0]} {frac_wan[cnt][1]} {frac_wan[cnt][2]} ", file=f)
    # [e] make zvo_geom.dat

    # [s] make zvo_dr.dat
    name_hr = f"{seedname}_hr.dat"
    Lx, Ly, Lz, orb_num = read_w90(name_hr)
    print(Lx, Ly, Lz, orb_num)

    Ncond = int(1 * Lx * Ly * Lz * orb_num / 3)

    MakeInputToml("input.toml", "zvo_geom.dat", name_hr, Lx, Ly, Lz, Ncond)
    # os.system("hwave input.toml")
    hwave.qlms.run(input_file="input.toml")

    data = np.load("output/green.npz")
    green = data["green"]

    name_out = "zvo_dr.dat"
    with open(name_out, "w") as fw:
        print("# zvo_dr.dat by wannier90 format", file=fw)
        print(f"{orb_num} ", file=fw)
        Nlattice = Lx * Ly * Lz
        print(f"{Nlattice} ", file=fw)
        for idx in range(Nlattice):
            print(" 1 ", end="", file=fw)
            if (idx + 1) % 15 == 0:
                print(" ", file=fw)
        if (idx + 1) % 15 != 0:
            print(" ", file=fw)
        
        from itertools import product
        for tmp_cnt_x, tmp_cnt_y, tmp_cnt_z in product(
            range(-Lx//2, Lx//2 + 1),
            range(-Ly//2, Ly//2 + 1),
            range(-Lz//2, Lz//2 + 1)
        ):
            cnt_x = ConvertCnt(Lx, tmp_cnt_x)
            cnt_y = ConvertCnt(Ly, tmp_cnt_y)
            cnt_z = ConvertCnt(Lz, tmp_cnt_z)
            cnt_tot = cnt_z + Lz * cnt_y + Lz * Ly * cnt_x
            for orb_j, orb_i in product(range(orb_num), range(orb_num)):
                tmp_val = green[cnt_tot][0][orb_i][0][orb_j] + green[cnt_tot][1][orb_i][1][orb_j]
                print(
                    f"{tmp_cnt_x} {tmp_cnt_y} {tmp_cnt_z} {orb_i + 1} {orb_j + 1} "
                    f"{tmp_val.real} {tmp_val.imag}",
                    file=fw,
                )

        # for tmp_cnt_x in range(-int(Lx / 2), int(Lx / 2) + 1, 1):
        #     cnt_x = ConvertCnt(Lx, tmp_cnt_x)
        #     for tmp_cnt_y in range(-int(Ly / 2), int(Ly / 2) + 1, 1):
        #         cnt_y = ConvertCnt(Ly, tmp_cnt_y)
        #         for tmp_cnt_z in range(-int(Lz / 2), int(Lz / 2) + 1, 1):
        #             cnt_z = ConvertCnt(Lz, tmp_cnt_z)
        #             cnt_tot = cnt_z + Lz * cnt_y + Lz * Ly * cnt_x
        #             for orb_j in range(orb_num):
        #                 for orb_i in range(orb_num):
        #                     tmp_val = green[cnt_tot][0][orb_i][0][orb_j]
        #                     tmp_val += green[cnt_tot][1][orb_i][1][orb_j]
        #                     print(
        #                         f" {tmp_cnt_x} {tmp_cnt_y} {tmp_cnt_z} {orb_i + 1} {orb_j + 1} {tmp_val.real} {tmp_val.imag} ",
        #                         file=fw,
        #                     )

def MakeInputToml(name_in, name_geom, name_hr, Lx, Ly, Lz, Ncond):
    """
    Create an input TOML file for the hwave.qlms module.

    This function generates a TOML configuration file for hwave.qlms, which is used
    to calculate Green's functions. The configuration includes:
    - Log settings for output control
    - Calculation mode (UHFk)
    - Parameters for the UHFk calculation
    - File paths for input and output

    Parameters
    ----------
    name_in : str
        Name of the input TOML file to be created.
    name_geom : str
        Name of the geometry file.
    name_hr : str
        Name of the Hamiltonian file.
    Lx, Ly, Lz : int
        Dimensions of the lattice.
    Ncond : int
        Number of conduction electrons.

    Examples
    --------
    >>> MakeInputToml("input.toml", "zvo_geom.dat", "seed_hr.dat", 10, 10, 10, 100)

    Notes
    -----
    The generated TOML file is used by hwave.qlms.run() to perform the UHFk
    calculation. The calculation parameters are set to default values that are
    typically suitable for most cases.
    """
    with open(name_in, 'w') as fw:
        print("[log]", file=fw)
        print("  print_level = 1", file=fw)
        print("  print_step = 10", file=fw)
        print("[mode]", file=fw)
        print("  mode = \"UHFk\"  ", file=fw)
        print("[mode.param]", file=fw)
        print(" # 2Sz = 0", file=fw)
        print(f"  Ncond = {Ncond}", file=fw)
        print("  IterationMax = 1000", file=fw)
        print("  EPS = 8", file=fw)
        print("  Mix = 0.5", file=fw)
        print("  T = 0.0", file=fw)
        print(f"  CellShape = [ {Lx}, {Ly}, {Lz} ]", file=fw)
        print("  SubShape = [ 1, 1, 1 ]", file=fw)
        print("[file]", file=fw)
        print("[file.input.interaction]", file=fw)
        print("  path_to_input = \"./\"", file=fw)
        print(f"  Geometry = \"{name_geom}\"", file=fw)
        print(f"  Transfer = \"{name_hr}\"", file=fw)
        print("[file.output]", file=fw)
        print("  path_to_output = \"output\"", file=fw)
        print("  energy = \"energy.dat\"", file=fw)
        print("  eigen = \"eigen\"", file=fw)
        print("  green = \"green\"", file=fw)

def read_w90(name_in):
    """
    Read the Wannier90 Hamiltonian file and extract lattice dimensions.

    Parameters
    ----------
    name_in : str
        Name of the Wannier90 Hamiltonian file.

    Returns
    -------
    tuple
        Dimensions of the lattice (Lx, Ly, Lz) and number of orbitals.

    Examples
    --------
    >>> Lx, Ly, Lz, orb_num = read_w90("seed_hr.dat")
    """
    with open(name_in, 'r') as f:
        l_strip = [s.strip() for s in f.readlines()[1:]]

    nr = int(l_strip[1])
    nints_per_line = 15
    skip_line = nr // nints_per_line
    if nr % nints_per_line != 0:
        skip_line += 1
    L_max = np.zeros((4), dtype=np.int64)
    L_min = np.zeros((4), dtype=np.int64)
    for idx, line in enumerate(l_strip[2 + skip_line:]):
        values = line.split()
        if len(values) == 0:
            break
        for cnt in range(len(L_max)):
            if int(values[cnt]) > L_max[cnt]:
                L_max[cnt] = int(values[cnt])
            if int(values[cnt]) < L_min[cnt]:
                L_min[cnt] = int(values[cnt])

    return (
        L_max[0] - L_min[0] + 1,
        L_max[1] - L_min[1] + 1,
        L_max[2] - L_min[2] + 1,
        L_max[3] - L_min[3],
    )

def ConvertCnt(max_L, cnt):
    """
    Convert a negative index to a positive one based on the maximum length.

    Parameters
    ----------
    max_L : int
        Maximum length of the dimension.
    cnt : int
        Index to be converted.

    Returns
    -------
    int
        Converted index.

    Examples
    --------
    >>> ConvertCnt(10, -1)
    9
    """
    if cnt < 0:
        cnt = cnt + max_L
    return cnt

def Convert2Frac(vec_lat, vec_wan):
    """
    Convert Cartesian coordinates to fractional coordinates.

    Parameters
    ----------
    vec_lat : ndarray
        Lattice vectors.
    vec_wan : ndarray
        Wannier center vectors.

    Returns
    -------
    ndarray
        Fractional coordinates of Wannier centers.

    Examples
    --------
    >>> vec_lat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> vec_wan = np.array([[0.5, 0.5, 0.5]])
    >>> Convert2Frac(vec_lat, vec_wan)
    array([[0.5, 0.5, 0.5]])
    """
    inv_vec_lat = np.linalg.inv(vec_lat)
    frac_wan = vec_wan @ inv_vec_lat
    # num_wan = vec_wan.shape[0]
    # frac_wan = np.zeros((num_wan, 3), dtype=np.float64)
    # for cnt in range(num_wan):
    #     for idx in range(3):
    #         tmp = 0.0
    #         for cnt_i in range(3):
    #             tmp += vec_wan[cnt][cnt_i] * inv_vec_lat[cnt_i][idx]
    #         frac_wan[cnt][idx] = tmp
    return frac_wan

def GetLatticeVec(file_name):
    """
    Extract lattice vectors from a file.

    Parameters
    ----------
    file_name : str
        Name of the file containing lattice vectors.

    Returns
    -------
    ndarray
        Lattice vectors.

    Examples
    --------
    >>> GetLatticeVec("seed.wout")
    array([[1.0, 0.0, 0.0],
           [0.0, 1.0, 0.0],
           [0.0, 0.0, 1.0]])
    """
    with open(file_name) as f:
        tmp = f.read().split("\n")
    
    # Find the line index where "Lattice Vectors" appears
    for i, line in enumerate(tmp):
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "Lattice" and parts[1] == "Vectors":
            start_idx = i + 1  # Next line contains the first lattice vector
            break
    else:
        # If "Lattice Vectors" is not found, return a zero matrix or handle as needed
        return np.zeros((3, 3), dtype=np.float64)

    # Read the next three lines and extract columns 1–3 as floats to form a 3×3 array
    vec_lat = np.array([
        [float(x) for x in tmp[j].split()[1:4]]
        for j in range(start_idx, start_idx + 3)
    ], dtype=np.float64)

    # for cnt in range(len(tmp)):
    #     tmp_2 = tmp[cnt].split()
    #     if len(tmp_2) > 0 and tmp_2[0] == "Lattice" and tmp_2[1] == "Vectors":
    #         cnt_s = cnt + 1
    #         break
    # vec_lat = np.zeros((3, 3), dtype=np.float64)
    # idx = 0
    # for cnt in range(cnt_s, cnt_s + 3):
    #     tmp_2 = tmp[cnt].split()
    #     vec_lat[idx][0] = float(tmp_2[1])
    #     vec_lat[idx][1] = float(tmp_2[2])
    #     vec_lat[idx][2] = float(tmp_2[3])
    #     idx += 1
    return vec_lat

def GetWanCentre(file_name):
    """
    Extract Wannier center vectors from a file.

    Parameters
    ----------
    file_name : str
        Name of the file containing Wannier center information.

    Returns
    -------
    ndarray
        Wannier center vectors.

    Examples
    --------
    >>> GetWanCentre("seed.wout")
    array([[0.5, 0.5, 0.5]])
    """
    with open(file_name) as f:
        lines = f.read().splitlines()

    # 1) Find the index of the line containing "Final State" using enumerate 
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "Final" and parts[1] == "State":
            start_idx = i + 1  # Next line is the first potential "WF centre" line
            break
    else:
        # If "Final State" not found, return an empty array
        return np.zeros((0, 3), dtype=np.float64)

    # 2) Collect all consecutive lines beginning with "WF centre" from start_idx 
    wan_parts = []
    for line in lines[start_idx:]:
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "WF" and parts[1] == "centre":
            wan_parts.append(parts)
        else:
            break

    # Number of wannier centres found
    num_wan = len(wan_parts)

    # 3) Initialize array to hold coordinates
    vec_wan = np.zeros((num_wan, 3), dtype=np.float64)

    # 4) Fill vec_wan based on the index given in the 5th token of each "WF centre" line
    for parts in wan_parts:
        # parts[4] holds the 1-based index of this Wannier function
        idx = int(parts[4]) - 1  # Convert to zero-based index
        # parts[6:9] hold the x,y,z coordinates (with trailing commas)
        x = float(parts[6].rstrip(','))
        y = float(parts[7].rstrip(','))
        z = float(parts[8].rstrip(','))
        vec_wan[idx] = (x, y, z)
    # with open(file_name) as f:
    #     tmp = f.read().split("\n")        
    # for cnt in range(len(tmp)):
    #     tmp_2 = tmp[cnt].split()
    #     if len(tmp_2) > 0 and tmp_2[0] == "Final" and tmp_2[1] == "State":
    #         cnt_s = cnt + 1
    #         break
    # cnt_e = cnt_s
    # for cnt in range(cnt_s, len(tmp)):
    #     tmp_2 = tmp[cnt].split()
    #     if len(tmp_2) > 0 and tmp_2[0] == "WF" and tmp_2[1] == "centre":
    #         cnt_e += 1
    #     else:
    #         break
    # num_wan = cnt_e - cnt_s
    # vec_wan = np.zeros((num_wan, 3), dtype=np.float64)
    # for cnt in range(cnt_s, cnt_e):
    #     tmp_2 = tmp[cnt].split()
    #     cnt_wan = int(tmp_2[4]) - 1
    #     vec_wan[cnt_wan][0] = float(tmp_2[6].rstrip(','))
    #     vec_wan[cnt_wan][1] = float(tmp_2[7].rstrip(','))
    #     vec_wan[cnt_wan][2] = float(tmp_2[8].rstrip(','))
    return vec_wan

def simple_variables(args):
    """
    Read input arguments and a TOML file (if exists) and set variables.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments.

    Returns
    -------
    dict
        Configuration dictionary.

    Examples
    --------
    >>> args = parser.parse_args(["conf.toml"])
    >>> simple_variables(args)
    {'base': {'seedname': 'example'}}
    """
    if os.path.exists(args.conf):
        try:
            import tomli
            print(f"Reading {args.conf}")
            with open(args.conf, "rb") as f:
                conf = tomli.load(f)
        except Exception as e:
            raise RuntimeError("Failed to read TOML file.") from e
    else:
        raise RuntimeError("TOML file does not exist.")

    return conf

if __name__ == "__main__":
    main()
