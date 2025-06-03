#!/usr/bin/env python
"""
RESPACK interface code for Wannier90.

This script processes Wannier90 output files to extract lattice vectors and Wannier
function centers, converting them to fractional coordinates for RESPACK input.

Usage
-----
1. Basic usage:
    python GetWannireCentre.py config.toml

2. With preprocessing flag:
    python GetWannireCentre.py config.toml --pp

The script requires a TOML configuration file with the following structure:
    [base]
    seedname = "your_seedname"  # Base name for Wannier90 files

The script will generate a zvo_geom.dat file containing:
    - Lattice vectors (3x3 matrix)
    - Number of Wannier centers
    - Fractional coordinates of Wannier centers

Examples
--------
1. Process a Wannier90 output file:
    $ python GetWannireCentre.py config.toml
    Reading config.toml
    # Output: zvo_geom.dat file is created

2. Process with preprocessing:
    $ python GetWannireCentre.py config.toml --pp
    Reading config.toml
    # Output: zvo_geom.dat file is created with preprocessing
"""

import argparse
import os
import numpy as np


def main():
    """
    Main function to process Wannier90 output and generate RESPACK geometry file.

    Reads configuration from TOML file, processes Wannier90 output (.wout file),
    and generates zvo_geom.dat containing lattice vectors and Wannier centers
    in fractional coordinates.

    Examples
    --------
    >>> # Example config.toml content:
    >>> # [base]
    >>> # seedname = "silicon"
    >>> # Running the script:
    >>> # $ python GetWannireCentre.py config.toml
    >>> # Output: zvo_geom.dat file is created with:
    >>> # - Lattice vectors from silicon.wout
    >>> # - Wannier centers in fractional coordinates
    """
    parser = argparse.ArgumentParser(
        description="RESPACK interface code for Wannier90"
    )
    parser.add_argument(
        'conf',
        type=str,
        metavar="conf.toml",
        help='Configuration TOML file containing input parameters'
    )
    parser.add_argument(
        "-pp",
        "--pp",
        action="store_true",
        help="Flag to enable preprocessing mode"
    )
    args = parser.parse_args()
    conf = simple_variables(args)
    seedname = conf["base"]["seedname"]

    # Process Wannier90 output file
    file_name = f"{seedname}.wout"
    vec_lat = GetLatticeVec(file_name)  # Get lattice vectors
    vec_wan = GetWanCentre(file_name)   # Get Wannier centers
    frac_wan = Convert2Frac(vec_lat, vec_wan)  # Convert to fractional coordinates

    # Write geometry data to zvo_geom.dat
    with open("zvo_geom.dat", 'w') as f:
        # Write lattice vectors
        for i in range(3):
            print(f"  {vec_lat[i][0]:f} {vec_lat[i][1]:f} {vec_lat[i][2]:f}", file=f)
        # Write number of Wannier centers
        print(f"  {frac_wan.shape[0]:d}", file=f)
        # Write fractional coordinates of Wannier centers
        for cnt in range(frac_wan.shape[0]):
            print(
                f"  {frac_wan[cnt][0]:f} {frac_wan[cnt][1]:f} {frac_wan[cnt][2]:f}",
                file=f
            )


def Convert2Frac(vec_lat: np.ndarray, vec_wan: np.ndarray) -> np.ndarray:
    """
    Convert Wannier centers from Cartesian to fractional coordinates.

    Parameters
    ----------
    vec_lat : np.ndarray
        Lattice vectors matrix (3x3) in Cartesian coordinates.
    vec_wan : np.ndarray
        Wannier centers coordinates (Nx3) in Cartesian coordinates.

    Returns
    -------
    np.ndarray
        Wannier centers coordinates (Nx3) in fractional coordinates.

    Notes
    -----
    The conversion is performed using the inverse of the lattice vector matrix.

    Examples
    --------
    >>> # Example lattice vectors (3x3 matrix)
    >>> vec_lat = np.array([
    ...     [5.0, 0.0, 0.0],
    ...     [0.0, 5.0, 0.0],
    ...     [0.0, 0.0, 5.0]
    ... ])
    >>> # Example Wannier centers (2x3 matrix)
    >>> vec_wan = np.array([
    ...     [2.5, 2.5, 2.5],
    ...     [1.0, 1.0, 1.0]
    ... ])
    >>> frac_coords = Convert2Frac(vec_lat, vec_wan)
    >>> print(frac_coords)
    [[0.5 0.5 0.5]
     [0.2 0.2 0.2]]
    """
    inv_vec_lat = np.linalg.inv(vec_lat)
    num_wan = vec_wan.shape[0]
    frac_wan = np.zeros((num_wan, 3), dtype=np.float64)

    # Convert each Wannier center to fractional coordinates
    for cnt in range(num_wan):
        for idx in range(3):
            tmp = 0.0
            for cnt_i in range(3):
                tmp += vec_wan[cnt][cnt_i] * inv_vec_lat[cnt_i][idx]
            frac_wan[cnt][idx] = tmp
    return frac_wan


def GetLatticeVec(file_name: str) -> np.ndarray:
    """
    Extract lattice vectors from Wannier90 output file.

    Parameters
    ----------
    file_name : str
        Path to the Wannier90 .wout file.

    Returns
    -------
    np.ndarray
        3x3 matrix containing lattice vectors in Cartesian coordinates.

    Raises
    ------
    RuntimeError
        If lattice vectors cannot be found in the file.

    Examples
    --------
    >>> # Example .wout file content:
    >>> # Lattice Vectors (Ang)
    >>> #    5.00000000    0.00000000    0.00000000
    >>> #    0.00000000    5.00000000    0.00000000
    >>> #    0.00000000    0.00000000    5.00000000
    >>> vec_lat = GetLatticeVec("silicon.wout")
    >>> print(vec_lat)
    [[5.0 0.0 0.0]
     [0.0 5.0 0.0]
     [0.0 0.0 5.0]]
    """
    with open(file_name) as f:
        tmp = f.read()
        tmp = tmp.split("\n")

    # Find the section containing lattice vectors
    cnt_s = None
    for cnt in range(len(tmp)):
        tmp_2 = tmp[cnt].split()
        if len(tmp_2) > 0 and tmp_2[0] == "Lattice" and tmp_2[1] == "Vectors":
            cnt_s = cnt + 1
            break

    if cnt_s is None:
        raise RuntimeError("Lattice vectors not found in Wannier90 output file")

    # Extract lattice vectors
    vec_lat = np.zeros((3, 3), dtype=np.float64)
    for idx, cnt in enumerate(range(cnt_s, cnt_s + 3)):
        tmp_2 = tmp[cnt].split()
        vec_lat[idx] = [float(tmp_2[1]), float(tmp_2[2]), float(tmp_2[3])]
    return vec_lat


def GetWanCentre(file_name: str) -> np.ndarray:
    """
    Extract Wannier function centers from Wannier90 output file.

    Parameters
    ----------
    file_name : str
        Path to the Wannier90 .wout file.

    Returns
    -------
    np.ndarray
        Nx3 matrix containing Wannier function centers in Cartesian coordinates.

    Raises
    ------
    RuntimeError
        If Wannier centers cannot be found in the file.

    Examples
    --------
    >>> # Example .wout file content:
    >>> # Final State
    >>> # WF centre and spread    1  (  2.50000000,  2.50000000,  2.50000000 )
    >>> # WF centre and spread    2  (  1.00000000,  1.00000000,  1.00000000 )
    >>> vec_wan = GetWanCentre("silicon.wout")
    >>> print(vec_wan)
    [[2.5 2.5 2.5]
     [1.0 1.0 1.0]]
    """
    with open(file_name) as f:
        tmp = f.read()
        tmp = tmp.split("\n")

    # Find the section containing Wannier centers
    cnt_s = None
    for cnt in range(len(tmp)):
        tmp_2 = tmp[cnt].split()
        if len(tmp_2) > 0 and tmp_2[0] == "Final" and tmp_2[1] == "State":
            cnt_s = cnt + 1
            break

    if cnt_s is None:
        raise RuntimeError("Wannier centers not found in Wannier90 output file")

    # Count number of Wannier centers
    cnt_e = cnt_s
    for cnt in range(cnt_s, len(tmp)):
        tmp_2 = tmp[cnt].split()
        if len(tmp_2) > 0 and tmp_2[0] == "WF" and tmp_2[1] == "centre":
            cnt_e += 1
        else:
            break

    num_wan = cnt_e - cnt_s
    vec_wan = np.zeros((num_wan, 3), dtype=np.float64)

    # Extract Wannier center coordinates
    for cnt in range(cnt_s, cnt_e):
        tmp_2 = tmp[cnt].split()
        cnt_wan = int(tmp_2[4]) - 1
        vec_wan[cnt_wan] = [
            float(tmp_2[6].rstrip(',')),
            float(tmp_2[7].rstrip(',')),
            float(tmp_2[8].rstrip(','))
        ]
    return vec_wan


def simple_variables(args: argparse.Namespace) -> dict:
    """
    Read and parse configuration from TOML file.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments containing the path to the TOML file.

    Returns
    -------
    dict
        Configuration dictionary loaded from TOML file.

    Raises
    ------
    RuntimeError
        If TOML file cannot be read or does not exist.

    Examples
    --------
    >>> # Example config.toml content:
    >>> # [base]
    >>> # seedname = "silicon"
    >>> args = argparse.Namespace(conf="config.toml")
    >>> conf = simple_variables(args)
    >>> print(conf["base"]["seedname"])
    'silicon'
    """
    if os.path.exists(args.conf):
        try:
            import tomli
            print(f"Reading {args.conf}")
            with open(args.conf, "rb") as f:
                conf = tomli.load(f)
        except Exception as e:
            raise RuntimeError(f"Failed to read TOML file: {str(e)}")
    else:
        raise RuntimeError(f"TOML file not found: {args.conf}")

    return conf


if __name__ == "__main__":
    main()
