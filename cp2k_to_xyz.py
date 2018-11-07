#!/usr/bin/env python
#
# @Author: Peng Lian
# @Date: 05/24/2018
# @Function: This script adds PBC lattices to the trajectory of CP2K output so that it can be viewed by ase gui.
# @Version: 0.0.1
# @Known Bugs: The Cell vectors are extraced from input file instead of from output file. 
#              Therefore, for cell optimization, the cell will be fixed instead of dyanamic.
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

def get_project_name(fcon):
    '''get project name from output file'''
    idx = find_flag_idx(fcon, ' GLOBAL| Project name')
    line = fcon[idx]
    proName = line.split()[-1]
    return proName

def neb_or_not(fcon):
    '''determine it neb output or not'''
    idx = find_flag_idx(fcon, ' NEB| Replica_env Setup. START')
    if idx > 0:
        return True
    else:
        return False

def get_Nreplica(fcon):
    '''get number of NEB replica from output file'''
    idx = find_flag_idx(fcon, ' NUMBER OF NEB REPLICA')
    line = fcon[idx]
    Nreplica = line.split()[-1]
    return int(Nreplica)

def get_trjfile(fcon):
    '''get trjfiel(s) from output file. If NEB, return a list of filenames, else a filename'''
    proName = get_project_name(fcon)

    if neb_or_not(fcon):      
        Nreplica = get_Nreplica(fcon)
        trjfile = ['%s-pos-Replica_nr_%d-1.xyz' % (proName, i) for i in range(1, 1+Nreplica)]
    else:
        trjfile = '%s-pos-1.xyz' % proName
    return trjfile

def get_inputfile(fcon):
    '''get from output file'''
    idx = find_flag_idx(fcon, ' CP2K| Input file name')
    line = fcon[idx]
    inp = line.split()[-1]
    return inp



def get_cell(fcon, inp=''):
    '''get inp from input file, and then get find cell info form inp'''

    # read inp
    if inp == '':
        inp = get_inputfile(fcon)
    
    # if faild to find inp file, throw an IOError
    try:
        fcon_inp = open(inp).readlines()
    except IOError:
        raise IOError
        exit()

    # remove the blanks in the front
    fcon_inp = [i.strip().upper() for i in fcon_inp]

    # get idx for different signs
    all_sign = find_flag_idx(fcon_inp, '&', findall=True)
    cell_sign = find_flag_idx(fcon_inp, '&CELL')
    cell_end_sign = all_sign[all_sign.index(cell_sign)+1]
    cell_section = fcon_inp[cell_sign:cell_end_sign]

    # get lattice parameter from cell section
    idx_ABC = find_flag_idx(cell_section, 'ABC ')
    idx_A = find_flag_idx(cell_section, 'A ')
    idx_B = find_flag_idx(cell_section, 'B ')
    idx_C = find_flag_idx(cell_section, 'C ')

    if idx_ABC == 0:
        a_list = cell_section[idx_A].split()[1:4]
        b_list = cell_section[idx_B].split()[1:4]
        c_list = cell_section[idx_C].split()[1:4]
    else:
        abc_list = cell_section[idx_ABC].split()[1:4]
        a_list = [abc_list[0], '0.0', '0.0']
        b_list = ['0.0', abc_list[1], '0.0']
        c_list = ['0.0', '0.0', abc_list[2]]

    return a_list, b_list, c_list

def convert_2_xyz(cp2k_output, fout, inp='', trj=''):

    fcon = open(cp2k_output).readlines()

    a, b, c = get_cell(fcon, inp=inp)
    comments = 'Lattice="%s %s %s" unit_cell=conventional pbc="T T T"' % (' '.join(a), ' '.join(b), ' '.join(c))

    # read trj file(s)
    if trj == '':
        trj = get_trjfile(fcon)

    # if NEB get the coordinates on the path
    if type(trj) == list:
        # get number of atoms
        fcon_trj = open(trj[0]).readlines()
        natoms = int(fcon_trj[0].strip())
        output_list = []
        for fname in trj:
            fcon_trj = open(fname).readlines()
            # find the coordinates
            idx = find_flag_idx(fcon_trj, ' i =', findall=True)
            for i in idx[-1:]:
                output_list.append(str(natoms) + '\n')
                output_list.append(comments + '\n')
                output_list += fcon_trj[i+1:i+1+natoms]
    else:
        fcon_trj = open(trj).readlines()
        natoms = int(fcon_trj[0].strip())
        
        # find the coordinates
        output_list = []
        idx = find_flag_idx(fcon_trj, ' i =', findall=True)
        for i in idx:
            output_list.append(str(natoms) + '\n')
            output_list.append(comments + '\n')
            output_list += fcon_trj[i+1:i+1+natoms]

    # write the xyz file
    open(fout, 'w').writelines(output_list)
    return

if __name__ == '__main__':
    # parser
    parser = argparse.ArgumentParser(description='adds PBC lattices to the trajectory of CP2K output.')
    parser.add_argument("--inp", nargs='?', type=str, default="", help="input file for cp2k.")
    parser.add_argument("--trj", nargs='?', type=str, default="", help="trajectory file for cp2k.")

    parser.add_argument("output", type=str, default="", help="OUTPUT file from cp2k. (with file type e.g. 'mol.out')")
    # all arguments
    args = parser.parse_args()
    outxyz = '%s.xyz' % args.output[:-4]
    
    # if catch the IOError, use the filename of xxx.out as the input file name.
    try:
        convert_2_xyz(args.output, outxyz, inp=args.inp, trj=args.trj)
    except IOError:
        inp = '%s.inp' % args.output[:-4]
        print('Using %s as input file.' % inp)
        convert_2_xyz(args.output, outxyz, inp=inp, trj=args.trj)

