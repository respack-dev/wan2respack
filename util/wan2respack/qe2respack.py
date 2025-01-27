#!/usr/bin/env python
# Copyright (c) 2017 Yusuke Nomura, 2018 Yuichi Motoyama, 2019 Terumasa Tadano, 2021 Jean-Baptiste Mor\'ee

"""
qe2respack.py -- convert Quantum ESPRESSO output into RESPACK input.
See `qe2respack.py --help`
"""

from __future__ import print_function

import argparse
import os
import os.path
import shutil
import struct
import sys
import xml.etree.ElementTree as ET

import numpy as np

if sys.version_info[0:2] < (2, 7):
    ET._ElementInterface.iter = ET._ElementInterface.getiterator


class Iotk_dat():
    def __init__(self, filename, endian=sys.byteorder):
        self.f = open(filename, 'rb')
        self.integer_nptype = {2: np.int16, 4: np.int32, 8: np.int64}
        self.integer_fmt = {2: 'h', 4: 'i', 8: 'q'}

        # Some environments do not have float128 dtype
        if hasattr(np, 'float128') and hasattr(np, 'complex256'):
            self.float_nptype = {4: np.float32, 8: np.float64, 16: np.float128}
            self.complex_nptype = {4: np.complex64, 8: np.complex128, 16: np.complex256}
        else:
            self.float_nptype = {4: np.float32, 8: np.float64}
            self.complex_nptype = {4: np.complex64, 8: np.complex128}

        self.float_fmt = {4: 'f', 8: 'd'}

        if endian == 'little':
            self.endian_fmt = '<'
        elif endian == 'big':
            self.endian_fmt = '>'
        else:
            raise RuntimeError('unknown endian, {0}'.format(endian))

    def __del__(self):
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, ex_type, ex_value, trace):
        self.f.close()

    def getattrs(self, line):
        ret = {}
        for word in line.split():
            if b'=' in word:
                w = word.split(b'=')
                ret[w[0]] = w[1].split(b'"')[1]
        return ret

    def load(self, name, rewind=False, raw=False):
        if type(name) is str:
            name = name.encode()
        if rewind:
            self.f.seek(0)
        s = self.f.readline().strip()
        while not s.startswith(b'<' + name):
            s = self.f.readline().strip()
        attrs = self.getattrs(s)
        nc = int(attrs.get(b'columns', '1'))
        nr = int(attrs[b'size']) // nc
        kind = int(attrs[b'kind'])

        self.f.read(4)  # size of the last line (write)

        if raw:
            b = struct.unpack(self.endian_fmt + 'i', self.f.read(4))[0] - 4
            self.f.read(4)  # some flag (complex or not???)
            return b, self.f.read(b)
        else:
            if attrs[b'type'] == b'integer':
                nptype = self.integer_nptype[kind]
                strfmt = self.endian_fmt + self.integer_fmt[kind]
                bcomplex = False
            elif attrs[b'type'] == b'real':
                nptype = self.float_nptype[kind]
                strfmt = self.endian_fmt + self.float_fmt[kind]
                bcomplex = False
            elif attrs[b'type'] == b'complex':
                nptype = self.complex_nptype[kind]
                strfmt = self.endian_fmt + self.float_fmt[kind]
                bcomplex = True

            self.f.read(8)  # data length and some flag (complex or not???)
            ret = np.zeros((nr, nc), nptype)
            for r in range(nr):
                for c in range(nc):
                    x = struct.unpack(strfmt, self.f.read(kind))[0]
                    if bcomplex:
                        y = struct.unpack(strfmt, self.f.read(kind))[0]
                        z = complex(x, y)
                        ret[r, c] = z
                    else:
                        ret[r, c] = x
            return ret


def calc_fermienergy_insulator(root):
    highest = -float('inf')
    lowest = float('inf')
    for kse in root.find('output').find('band_structure').iter('ks_energies'):
        es = map(float, kse.find('eigenvalues').text.strip().split())
        os = map(float, kse.find('occupations').text.strip().split())
        for e, o in zip(es, os):
            if o == 0.0:
                lowest = min(lowest, e)
                break
            else:
                highest = max(highest, e)
    return 0.5 * (highest + lowest)


def band_structure_info(root, oldxml=False, lsda=False):
    if oldxml:
        child = root.find('BAND_STRUCTURE_INFO')
        num_k = int(child.find('NUMBER_OF_K-POINTS').text)
        eFermi = float(child.find('FERMI_ENERGY').text)
        if lsda:
            num_b = [int(child.find('NUMBER_OF_BANDS').attrib['UP']),
                     int(child.find('NUMBER_OF_BANDS').attrib['DW'])]
        else:
            num_b = int(child.find('NUMBER_OF_BANDS').text)
    else:
        child = root.find('output').find('band_structure')
        num_k = int(child.find('nks').text)
        if lsda:
            num_b_up = int(child.find('nbnd_up').text)
            num_b_dw = int(child.find('nbnd_dw').text)
            num_b = [num_b_up, num_b_dw]
        else:
            num_b = int(child.find('nbnd').text)
        if child.find('fermi_energy') == None:
            # eFermi = float(child.find('highestOccupiedLevel').text)
            eFermi = calc_fermienergy_insulator(root)
        else:
            eFermi = float(child.find('fermi_energy').text)
    return num_k, num_b, eFermi


def kpoint_coords(root, num_k, oldxml=False):
    k_vec = np.zeros((3, num_k))
    if oldxml:
        child = root.find('EIGENVALUES')
        for i in range(num_k):
            ev = child.find('K-POINT.{0}'.format(i + 1))
            k_vec[:, i] = [float(x) for x in ev.find('K-POINT_COORDS').text.strip().split()]
    else:
        for i, kse in enumerate(root.find('output').iter('ks_energies')):
            k_vec[:, i] = [float(x) for x in kse.find('k_point').text.strip().split()]
    return k_vec


def number_of_GK_vectors(root, num_k, oldxml=False):
    if oldxml:
        child = root.find('EIGENVECTORS')
        num_Gk = [int(child.find('K-POINT.{0}'.format(i + 1))
                      .find('NUMBER_OF_GK-VECTORS').text) for i in range(num_k)]
    else:
        num_Gk = [int(kse.find('npw').text) for kse in root.find('output').iter('ks_energies')]
    return num_Gk


def latvectors(root, oldxml=False):
    A = np.zeros((3, 3))
    if oldxml:
        child = root.find('CELL')
        celldm = float(child.find('LATTICE_PARAMETER').text)
        cc = child.find('DIRECT_LATTICE_VECTORS')
        for i in range(3):
            A[i, :] = [float(x) for x in cc.find('a{0}'.format(i + 1)).text.split()]
    else:
        input_tag = root.find('input')
        output_tag = root.find('output')
        if input_tag is not None:
            child = input_tag.find('atomic_structure')
            celldm = float(child.attrib['alat'])
            cell = child.find('cell')
            # The celldm should be alat when "CELL_PARAMETERS alat" is used, but it isn't
            # in the <input> tag probably due to a bug in QE. The value in the 'alat' attribute
            # of the <atomic_structure> tag inside <output> is alat (= celldm(1)).
            # So, let's update celldm as follows:
            celldm = float(output_tag.find('atomic_structure').attrib['alat'])
        else:
            # In some cases, <input> tag may be absent. 
            # If so, let's parse data from <output> instead.
            child = output_tag.find('atomic_structure')
            celldm = float(child.attrib['alat'])
            cell = child.find('cell')
        for i in range(3):
            A[i, :] = [float(x) for x in cell.find('a{0}'.format(i + 1)).text.split()]
    return A, celldm


def atoms_information(root, oldxml=False):
    if oldxml:
        child = root.find('IONS')
        n_atoms = int(child.find('NUMBER_OF_ATOMS').text)
        atom_symbs = ['' for i in range(n_atoms)]
        atom_positions = np.zeros((3, n_atoms))
        for i in range(n_atoms):
            atom = child.find('ATOM.{0}'.format(i + 1)).attrib
            atom_symbs[i] = atom['SPECIES'].strip()
            atom_positions[:, i] = [float(x) for x in atom['tau'].split()]
    else:
        child = root.find('output').find('atomic_structure')
        n_atoms = int(child.attrib['nat'])
        atom_symbs = ['' for i in range(n_atoms)]
        atom_positions = np.zeros((3, n_atoms))
        for i, atom in enumerate(child.find('atomic_positions').iter('atom')):
            atom_symbs[i] = atom.attrib['name']
            atom_positions[:, i] = [float(x) for x in atom.text.split()]
    return atom_positions, atom_symbs


def wfc_cutoff(root, oldxml=False):
    if oldxml:
        child = root.find('PLANE_WAVES')
        Ecut_for_psi = float(child.find('WFC_CUTOFF').text)
    else:
        input_tag = root.find('input')
        if input_tag is not None:
            child = root.find('input').find('basis')
            Ecut_for_psi = float(child.find('ecutwfc').text)
        else:
            child = root.find('output').find('basis_set')
            Ecut_for_psi = float(child.find('ecutwfc').text)
    return Ecut_for_psi


def symmetry(root, oldxml=False):
    if oldxml:
        child = root.find('SYMMETRIES')
        n_sym = int(child.find('NUMBER_OF_SYMMETRIES').text)
        ftau = np.zeros((3, n_sym))
        mat_sym = np.zeros((3, 3, n_sym), dtype=int)
        for i in range(n_sym):
            sym = child.find('SYMM.{0}'.format(i + 1))
            rot = sym.find('ROTATION').text.strip().split('\n')
            for j in range(3):
                mat_sym[j, :, i] = [int(x) for x in rot[j].split()]
            ftau[:, i] = [float(x) for x in sym.find('FRACTIONAL_TRANSLATION').text.split()]
        pass
    else:
        child = root.find('output').find('symmetries')
        n_sym = int(child.find('nsym').text)
        ftau = np.zeros((3, n_sym))
        mat_sym = np.zeros((3, 3, n_sym), dtype=int)
        for i, sym in enumerate(child.iter('symmetry')):
            ftau[:, i] = [float(x) for x in sym.find('fractional_translation').text.split()]
            rot = sym.find('rotation').text.strip().split('\n')
            for j in range(3):
                mat_sym[j, :, i] = [int(round(float(x))) for x in rot[j].split()]
            if i == n_sym - 1:
                break
    return mat_sym, ftau


def eigenvalues(dirname, num_k, num_b, oldxml=False, lsda=False):
    """
     Parse and return eigenvalues.
     If lsda = True, the first half of evs[k,:] corresponds to the spin 'up' state
     and the last half of env[k,:] corresponds to the 'down' state.
     If lsda = False (non spin-polarized calculation and noncollinear calculation),
     the eigenvalues are sorted in the ascending order.
    """
    if lsda:
        evs = np.zeros((num_k, np.sum(num_b)))
    else:
        evs = np.zeros((num_k, num_b))

    if oldxml:
        if lsda:
            for k in range(num_k):
                evs_tmp = []
                for ispin in range(2):
                    tree = ET.parse(os.path.join(dirname, 'K{0:0>5}/eigenval{1}.xml'.format(k + 1, ispin + 1)))
                    root = tree.getroot()
                    child = root.find('EIGENVALUES')
                    evs_tmp.extend([float(x) for x in child.text.strip().split()])
                evs[k, :] = np.array(evs_tmp)
        else:
            for k in range(num_k):
                tree = ET.parse(os.path.join(dirname, 'K{0:0>5}/eigenval.xml'.format(k + 1)))
                root = tree.getroot()
                child = root.find('EIGENVALUES')
                evs[k, :] = [float(x) for x in child.text.strip().split()]
    else:
        tree = ET.parse(os.path.join(dirname, 'data-file-schema.xml'))
        root = tree.getroot()
        child = root.find('output').find('band_structure')
        for k, kse in enumerate(child.iter('ks_energies')):
            evs[k, :] = [float(x) for x in kse.find('eigenvalues').text.strip().split()]

    return evs


# def occupations(dirname, num_k, num_b, oldxml=False):
#    occs = np.zeros((num_k,num_b))
#    if oldxml:
#        for k in range(num_k):
#            tree = ET.parse(os.path.join(dirname, 'K{0:0>5}/eigenval.xml'.format(k+1)))
#            root = tree.getroot()
#            child = root.find('OCCUPATIONS')
#            occs[k,:] = [float(x) for x in child.text.strip().split()]
#    else:
#        tree = ET.parse(os.path.join(dirname, 'data-file-schema.xml'))
#        root = tree.getroot()
#        child = root.find('output').find('band_structure')
#        for k,kse in enumerate(child.iter('ks_energies')):
#            occs[k,:] = [float(x) for x in kse.find('occupations').text.strip().split()]
#    return occs

def spin_orbit_info(root, oldxml=False):
    if oldxml:
        child = root.find('SPIN')
        is_SpinOrbit = (child.find('SPIN-ORBIT_CALCULATION').text.strip() == 'T')
        is_noncolin = (child.find('NON-COLINEAR_CALCULATION').text.strip() == 'T')
        is_LSDA = (child.find('LSDA').text.strip() == 'T')
    else:
        child = root.find('output').find('magnetization')
        is_noncolin = (child.find('noncolin').text == 'true')
        is_SpinOrbit = (child.find('spinorbit').text == 'true')
        is_LSDA = (child.find('lsda').text == 'true')
    return is_noncolin, is_SpinOrbit, is_LSDA


def qe2respack(dirname, endian=sys.byteorder):
    if endian == 'little':
        endian_fmt = '<'
    elif endian == 'big':
        endian_fmt = '>'
    else:
        raise RuntimeError('unknown endian, {0}'.format(endian))

    xmlfile = os.path.join(dirname, 'data-file-schema.xml')
    if os.path.exists(xmlfile):
        print('New style QE output XML file (data-file-schema.xml) is found.')
        oldxml = False
    else:
        xmlfile = os.path.join(dirname, 'data-file.xml')
        if os.path.exists(xmlfile):
            print('Old style QE output XML file (data-file.xml) is found.')
            oldxml = True
        else:
            raise RuntimeError('QE output XML file is NOT found in {0}.'.format(dirname))

    print('loading {0}'.format(xmlfile))
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    is_noncolin, is_SpinOrbit, is_LSDA = spin_orbit_info(root, oldxml=oldxml)
    num_k, num_b, eFermi = band_structure_info(root, oldxml=oldxml, lsda=is_LSDA)
    print('num_k = {0}'.format(num_k))
    print('num_b = {0}'.format(num_b))
    print('eFermi = {0}'.format(eFermi))

    k_vec = kpoint_coords(root, num_k, oldxml=oldxml)
    num_Gk = number_of_GK_vectors(root, num_k, oldxml=oldxml)

    A, celldm = latvectors(root, oldxml=oldxml)
    Ainv = np.linalg.inv(A.transpose())

    atom_positions, atom_symbs = atoms_information(root, oldxml=oldxml)
    n_atoms = len(atom_symbs)

    Ecut_for_psi = wfc_cutoff(root, oldxml=oldxml)

    mat_sym, ftau = symmetry(root, oldxml=oldxml)
    n_sym = ftau.shape[1]
    print('n_sym = {0}'.format(n_sym))

    mat_sym, ftau = symmetry(root, oldxml=oldxml)
    print('is_LSDA = {0}'.format(is_LSDA))
    print('is_noncolin = {0}'.format(is_noncolin))
    print('is_SpinOrbit = {0}'.format(is_SpinOrbit))
    ncomp = 1
    if is_noncolin or is_SpinOrbit: ncomp = 2

    # end of read XML file

    k_vec = np.dot(A, k_vec) / celldm

    if not os.path.exists('dir-wfn'):
        os.mkdir('dir-wfn')

    print('generating dir-wfn/dat.sample-k')
    with open('./dir-wfn/dat.sample-k', 'w') as f:
        f.write('{0}\n'.format(num_k))
        for i in range(num_k):
            f.write('{0} {1} {2}\n'.format(k_vec[0, i], k_vec[1, i], k_vec[2, i]))

    print('generating dir-wfn/dat.lattice')
    with open('./dir-wfn/dat.lattice', 'w') as f:
        for i in range(3):
            f.write('{0:20.15f} {1:20.15f} {2:20.15f}\n'.format(A[i, 0], A[i, 1], A[i, 2]))

    print('generating dir-wfn/dat.bandcalc')
    with open('./dir-wfn/dat.bandcalc', 'w') as f:
        f.write('{0}\n'.format(Ecut_for_psi * 2.0))  # Hartree => Rydberg
        f.write('{0}\n'.format(eFermi))
        f.write('0.0\n')

    print('generating dir-wfn/dat.nkm')
    with open('./dir-wfn/dat.nkm', 'w') as f:
        for i in range(num_k):
            f.write('{0}\n'.format(num_Gk[i]))

    print('generating dir-wfn/dat.eigenvalue')
    eigvals = eigenvalues(dirname, num_k, num_b, oldxml=oldxml, lsda=is_LSDA)
    if is_LSDA:
        with open('./dir-wfn/dat.eigenvalue.up', 'w') as f:
            f.write('{0}\n'.format(num_b[0]))
            for k in range(num_k):
                for i in range(num_b[0]):
                    f.write(str(eigvals[k, i]))
                    f.write('\n')

        with open('./dir-wfn/dat.eigenvalue.dn', 'w') as f:
            f.write('{0}\n'.format(num_b[1]))
            for k in range(num_k):
                for i in range(num_b[1]):
                    f.write(str(eigvals[k, i + num_b[0]]))
                    f.write('\n')
    else:
        with open('./dir-wfn/dat.eigenvalue', 'w') as f:
            f.write('{0}\n'.format(num_b))
            for k in range(num_k):
                for i in range(num_b):
                    f.write(str(eigvals[k, i]))
                    f.write('\n')

    #    print('generating dir-wfn/dat.occ')
    #    occs = occupations(dirname, num_k, num_b, oldxml=oldxml)
    #    with open('./dir-wfn/dat.occ', 'w') as f:
    #        for k in range(num_k):
    #            for i in range(num_b):
    #                f.write(str(occs[k,i]))
    #                f.write('\n')

    print('generating dir-wfn/dat.atom_position')
    with open('./dir-wfn/dat.atom_position', 'w') as f:
        f.write('{0}\n'.format(n_atoms))
        X = np.dot(Ainv, atom_positions)
        for i in range(n_atoms):
            f.write('{0} {1} {2} {3}\n'.format(atom_symbs[i], X[0, i], X[1, i], X[2, i]))

    print('generating dir-wfn/dat.symmetry')
    with open('./dir-wfn/dat.symmetry', 'w') as f:
        for i in range(1, 1001):
            tmp = (ftau * i) - (ftau * i).round()
            if (abs(tmp) <= 1.0e-7).all():
                f.write('{0}\n'.format(n_sym))
                f.write('{0}\n'.format(i))
                tau = -((ftau * i).round())
                break
        for i in range(n_sym):
            mat = np.zeros((3, 3))
            mat[:, :] = mat_sym[:, :, i]
            mat[:, :] = np.linalg.inv(mat)

            for j in range(3):
                for k in range(3):
                    f.write('{0} '.format(int(mat[j, k])))
            f.write('\n')
            f.write('{0} {1} {2}\n'.format(int(tau[0, i]), int(tau[1, i]), int(tau[2, i])))

    print('generating dir-wfn/dat.kg')
    with open('./dir-wfn/dat.kg', 'w') as f:
        if oldxml:
            for k in range(num_k):
                f.write('{0}\n'.format(num_Gk[k]))
                with Iotk_dat(os.path.join(dirname, 'K{0:0>5}/gkvectors.dat'.format(k + 1)), endian=endian) as inp:
                    dat = inp.load('GRID')
                    nr, nc = dat.shape
                    for r in range(nr):
                        f.write('{0} {1} {2}\n'.format(dat[r, 0], dat[r, 1], dat[r, 2]))
        else:
            for k in range(num_k):
                f.write('{0}\n'.format(num_Gk[k]))
                if is_LSDA:
                    fname_to_parse = os.path.join(dirname, 'wfcup{0}.dat'.format(k + 1))
                else:
                    fname_to_parse = os.path.join(dirname, 'wfc{0}.dat'.format(k + 1))

                with open(fname_to_parse, 'rb') as inp:
                    fmt = endian_fmt + 'i'
                    fmt2 = endian_fmt + 'iii'
                    # skip three blocks
                    for i in range(3):
                        # Get the record length from the heading 4 bytes
                        n = struct.unpack(fmt, inp.read(4))[0]
                        # Read the body record + tailing 4 bytes
                        inp.read(n + 4)

                    n = struct.unpack(fmt, inp.read(4))[0]
                    nr = n // 12  # three integers
                    for r in range(nr):
                        kg = struct.unpack(fmt2, inp.read(12))
                        f.write('{0} {1} {2}\n'.format(kg[0], kg[1], kg[2]))

    print(ncomp)

    if oldxml:
        if is_LSDA:
            for ispin, spin in enumerate(['up', 'dn']):
                print('generating dir-wfn/dat.wfn.%s' % spin)
                with open('./dir-wfn/dat.wfn.%s' % spin, 'wb') as f:
                    f.write(struct.pack(endian_fmt + 'i', 4))
                    f.write(struct.pack(endian_fmt + 'i', ncomp))
                    f.write(struct.pack(endian_fmt + 'i', 4))
                    for k in range(num_k):
                        with Iotk_dat(os.path.join(dirname, 'K{0:0>5}/evc{1}.dat'.format(k + 1, ispin + 1)),
                                      endian=endian) as inp:
                            for ib in range(num_b[ispin]):
                                size, dat = inp.load('evc.{0}'.format(ib + 1), raw=True)
                                f.write(struct.pack(endian_fmt + 'i', size))
                                f.write(dat)
                                f.write(struct.pack(endian_fmt + 'i', size))
        else:
            print('generating dir-wfn/dat.wfn')
            with open('./dir-wfn/dat.wfn', 'wb') as f:
                # with open('./dir-wfn/dat.wfn_ascii', 'w') as f2:
                f.write(struct.pack(endian_fmt + 'i', 4))
                f.write(struct.pack(endian_fmt + 'i', ncomp))
                f.write(struct.pack(endian_fmt + 'i', 4))
                for k in range(num_k):
                    if not (is_SpinOrbit or is_noncolin):
                        with Iotk_dat(os.path.join(dirname, 'K{0:0>5}/evc.dat'.format(k + 1)), endian=endian) as inp:
                            for ib in range(num_b):
                                size, dat = inp.load('evc.{0}'.format(ib + 1), raw=True)
                                f.write(struct.pack(endian_fmt + 'i', size))
                                f.write(dat)
                                f.write(struct.pack(endian_fmt + 'i', size))
                    else:
                        # noncollinear case
                        with Iotk_dat(os.path.join(dirname, 'K{0:0>5}/evc1.dat'.format(k + 1)),
                                      endian=endian) as inp1:
                            with Iotk_dat(os.path.join(dirname, 'K{0:0>5}/evc2.dat'.format(k + 1)),
                                          endian=endian) as inp2:
                                for ib in range(num_b):
                                    size1, dat1 = inp1.load('evc.{0}'.format(ib + 1), raw=True)
                                    size2, dat2 = inp2.load('evc.{0}'.format(ib + 1), raw=True)
                                    f.write(struct.pack(endian_fmt + 'i', size1 + size2))
                                    f.write(dat1)
                                    #                            f.write(struct.pack(endian_fmt+'i', size1))
                                    #                            f.write(struct.pack(endian_fmt+'i', size2))
                                    f.write(dat2)
                                    f.write(struct.pack(endian_fmt + 'i', size1 + size2))

    else:
        # New format
        if is_LSDA:
            slabel_qe = ['up', 'dw']
            for ispin, spin in enumerate(['up', 'dn']):
                print('generating dir-wfn/dat.wfn.%s' % spin)
                with open('./dir-wfn/dat.wfn.%s' % spin, 'wb') as f:
                    f.write(struct.pack(endian_fmt + 'i', 4))
                    f.write(struct.pack(endian_fmt + 'i', ncomp))
                    f.write(struct.pack(endian_fmt + 'i', 4))
                    for k in range(num_k):
                        with open(os.path.join(dirname, 'wfc{0}{1}.dat'.format(slabel_qe[ispin], k + 1)), 'rb') as inp:
                            # skip four blocks
                            fmt = endian_fmt + 'i'
                            n = struct.unpack(fmt, inp.read(4))[0]
                            inp.read(n + 4)
                            n = struct.unpack(fmt, inp.read(4))[0]
                            inp.read(n + 4)
                            n = struct.unpack(fmt, inp.read(4))[0]
                            inp.read(n + 4)
                            n = struct.unpack(fmt, inp.read(4))[0]
                            inp.read(n + 4)
                            for ib in range(num_b[ispin]):
                                n = struct.unpack(fmt, inp.read(4))[0]
                                dat = inp.read(n)
                                inp.read(4)
                                f.write(struct.pack(endian_fmt + 'i', n))
                                f.write(dat)
                                f.write(struct.pack(endian_fmt + 'i', n))
        else:
            print('generating dir-wfn/dat.wfn')
            with open('./dir-wfn/dat.wfn', 'wb') as f:
                # with open('./dir-wfn/dat.wfn_ascii', 'w') as f2:
                f.write(struct.pack(endian_fmt + 'i', 4))
                f.write(struct.pack(endian_fmt + 'i', ncomp))
                f.write(struct.pack(endian_fmt + 'i', 4))
                for k in range(num_k):
                    with open(os.path.join(dirname, 'wfc{0}.dat'.format(k + 1)), 'rb') as inp:
                        # skip four blocks
                        fmt = endian_fmt + 'i'
                        n = struct.unpack(fmt, inp.read(4))[0]
                        inp.read(n + 4)
                        n = struct.unpack(fmt, inp.read(4))[0]
                        inp.read(n + 4)
                        n = struct.unpack(fmt, inp.read(4))[0]
                        inp.read(n + 4)
                        n = struct.unpack(fmt, inp.read(4))[0]
                        inp.read(n + 4)
                        for ib in range(num_b):
                            n = struct.unpack(fmt, inp.read(4))[0]
                            dat = inp.read(n)
                            # print(ib,n)
                            # print(dat)
                            inp.read(4)
                            f.write(struct.pack(endian_fmt + 'i', n))
                            f.write(dat)
                            f.write(struct.pack(endian_fmt + 'i', n))
                            # f2.write(int(dat))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='convert Quantum Espresso output into RESPACK input.')
    parser.add_argument('QE_output_dir',
                        help='output directory of Quantum Espresso (where data-file-schema.xml exists).')
    args = parser.parse_args()

    if os.path.exists('dir-wfn'):
        print('dir-wfn/ already exists. qe2respack.py moves this to dir-wfn.backup/ as a backup.')
        if os.path.exists('dir-wfn.backup'):
            shutil.rmtree('dir-wfn.backup')
        shutil.move('dir-wfn', 'dir-wfn.backup')

    qe2respack(args.QE_output_dir)
