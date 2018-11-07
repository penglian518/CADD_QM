#!/usr/bin/env python
#
# @Author: Peng Lian
# @Date: 05/24/2018
# @Function: This script extracts coordinates and parameters from the output file of Quantum Espresso.
#

import argparse
import numpy as np

def find_flag_idx(fcon, startstring, findall=False):
    '''
    turn on findall will return a list, otherwise a number.
    '''
    # find all or not
    if findall:
        idx = []
    else:
        idx = 0

    conter = 0
    while conter < len(fcon):
        line = fcon[conter]
        if line.startswith(startstring):
            if findall:
                idx.append(conter)
                conter += 1
            else:
                idx = conter
                conter += len(fcon)
        else:
            conter += 1
    return idx

def get_natoms(fcon):
    idx = find_flag_idx(fcon, '     number of atoms')
    line = fcon[idx]
    natoms = int(line.split()[-1])
    return natoms

def get_cell(fcon):
    ang2au = 1.8897259885789

    # get lattice parameter
    idx = find_flag_idx(fcon, '     lattice parameter')
    line = fcon[idx]
    lattice = float(line.split()[-2])

    # get lattice vectors and convert to angstrom
    idx = find_flag_idx(fcon, '     crystal axes:')
    line = fcon[idx+1]
    a = np.array([ lattice/ang2au * float(i) for i in line.split()[3:6]])
    line = fcon[idx+2]
    b = np.array([ lattice/ang2au * float(i) for i in line.split()[3:6]])
    line = fcon[idx+3]
    c = np.array([ lattice/ang2au * float(i) for i in line.split()[3:6]])

    # compute the angles. (not necessary)
    cosAB = np.dot(a, b)/np.linalg.norm(a)/np.linalg.norm(b)
    gamma = np.degrees(np.arccos(cosAB))
    cosAC = np.dot(a, c)/np.linalg.norm(a)/np.linalg.norm(c)
    beta = np.degrees(np.arccos(cosAC))
    cosBC = np.dot(b, c)/np.linalg.norm(b)/np.linalg.norm(c)
    alpha = np.degrees(np.arccos(cosBC))

    return a, b, c, alpha, beta, gamma

def convert_2_xyz(qe_output, fout):
    fcon = open(qe_output).readlines()
    natoms = get_natoms(fcon)
    a, b, c, alpha, beta, gamma = get_cell(fcon)
    comments = 'Lattice="%s %s %s" unit_cell=conventional pbc="T T T"' % (' '.join(list(a.astype(str))), ' '.join(list(b.astype(str))), ' '.join(list(c.astype(str))))

    # find the coordinates
    output_list = []
    idx = find_flag_idx(fcon, 'ATOMIC_POSITIONS', findall=True)
    for i in idx:
        output_list.append(str(natoms) + '\n')
        output_list.append(comments + '\n')
        output_list += fcon[i+1:i+1+natoms]

    # write the xyz file
    open(fout, 'w').writelines(output_list)
    return

if __name__ == '__main__':
    # parser
    parser = argparse.ArgumentParser(description='get coordinates from the output file of Quantum Espresso.')
    parser.add_argument("qe_output", type=str, default="", help="Output file from QE. (with file type e.g. 'mol.out')")
    # all arguments
    args = parser.parse_args()
    outxyz = '%s.xyz' % args.qe_output[:-4]
    convert_2_xyz(args.qe_output, outxyz)

