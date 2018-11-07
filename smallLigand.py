#!/usr/bin/env python
#
# @purpose
#   manipulate small ligands
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Mar 15 2017
#


import os, shutil, logging, argparse, random, subprocess, math
import numpy as np
import pandas as pd

try:
    import vmd
except:
    pass

import VMD

class smallLigand():
    def __init__(self):
        self.element_charge = {
            'H' :  1,  'H1':  1,  'H2'  :  1,
            'C' : 6,
            'N' : 7,
            'O': 8, 'OH2': 8,
        }

        self.switch = 10.0
        self.cutoff = 12.0
        self.dielectronic = 1.0
        self.tempname = 'tmp'
        self.atomselection = 'all'
        self.obminimize = '/home/p6n/anaconda2/bin/obminimize'
        self.obabel = '/home/p6n/anaconda2/bin/obabel'
        self.namd2 = '/opt/NAMD_2.12b1_Linux-x86_64-multicore/namd2'
        self.namd2ff = '/home/p6n/tools/vmd-1.9.2/lib/vmd/plugins/noarch/tcl/readcharmmpar1.2/par_all27_prot_lipid_na.inp'



    def make_rand_vector(self, dims=3):
        vec = [random.gauss(0, 1) for i in range(dims)]
        mag = sum(x**2 for x in vec)**.5
        return [x / mag for x in vec]


    def coor2cmat(self, coormatrix):
        '''
        convert a coordinate matrix to Coulomb matrix file
        :param coormatrix: input coormatrix, in the following format:

         [nuclear charge,     x  ,     y  ,     z   ]
        matrix([[ 16.   ,   1.32 ,  -0.561,   1.693],
                [  1.   ,   1.018,  -0.737,   0.795],
                [  1.   ,   2.218,  -0.263,   1.533],
                [ 16.   ,   0.278,   1.922,  -0.8  ],
                [  1.   ,   0.701,   1.392,  -1.44 ],
                [  1.   ,  -0.444,   1.406,  -0.585]])

        '''

        # get xyz matrix
        xyz = coormatrix[:, 1:]
        m, n = xyz.shape

        # cal diagonal element
        z = coormatrix[:, 0]
        zz = z * z.T
        #cmat_diagonal = np.power(np.multiply(zz, np.eye(m)), 2.4) * 0.5
        cmat_diagonal = np.power(np.multiply(zz, np.eye(m)), 0.5)

        # cal the distance matrix
        dist_list = []
        for i in range(m):
            # cal dist between atom i and others
            dist_i = np.power(np.sum(np.power(xyz - xyz[i, :], 2), axis=1), 0.5)
            dist_list.append(dist_i)
        dist_matrix = np.concatenate(dist_list, axis=1)

        # cal off diagonal element
        cmat_offdiagonal = np.multiply(zz, 1 / (dist_matrix + np.eye(m)) - np.eye(m))

        # combine
        cmat_combine = cmat_diagonal + cmat_offdiagonal

        # keep lower corner
        cmat = np.tril(cmat_combine, k=0)

        # convert to vector
        ind_tril = np.tril_indices(m)
        cmat_vector = cmat[ind_tril]

        '''
        cmat_vector = []
        i = 0
        while i < m:
            j = 0
            while j <= i:
                cmat_vector.append(cmat[i][j])
                j += 1
            i += 1
        '''

        return cmat, cmat_vector

    def perturb_coor(self, coorarray, offset=0.1):
        '''
        add a rondomly generated offset to the coordinate

        :param coorarray:

         array([[ 0.278,  1.922, -0.8  ],
                [ 0.701,  1.392, -1.44 ],
                [-0.444,  1.406, -0.585]])

        :param offset: better not larger than 0.5, as may generate overlap atoms
        :return: perturbed coordinates array
        '''

        offset_matrix = np.random.random(coorarray.shape) * offset
        # combine
        coorarray_off = coorarray.astype(float) + offset_matrix

        return coorarray_off

    def rotate_coor(self, coorarray):
        '''
        rotate the molecule in X, Y, Z axis by a rondom degree

        :param coorarray:

         array([[ 0.278,  1.922, -0.8  ],
                [ 0.701,  1.392, -1.44 ],
                [-0.444,  1.406, -0.585]])

        :return: perturbed coordinates array
        '''

        # gene random theta (-pi to pi)
        theta = random.uniform(-1.0, 1.0) * math.pi

        # gene random axis
        x, y, z = np.array(self.make_rand_vector(dims=3))

        Rxyz = np.matrix([
            [x*x*(1-math.cos(theta))+math.cos(theta), x*y*(1-math.cos(theta))-z*math.sin(theta), x*z*(1-math.cos(theta))+y*math.sin(theta)],
            [x*y*(1-math.cos(theta))+z*math.sin(theta), y*y*(1-math.cos(theta))+math.cos(theta), y*z*(1-math.cos(theta))-x*math.sin(theta)],
            [x*z*(1-math.cos(theta))-y*math.sin(theta), y*z*(1-math.cos(theta))+x*math.sin(theta), z*z*(1-math.cos(theta))+math.cos(theta)]
        ])

        ## Rx(theta), Ry(theta), Rz(theta)
        #Rx = np.matrix([[1., 0., 0.], [0., math.cos(theta), -1. * math.sin(theta)], [0., math.sin(theta), math.cos(theta)]])
        #Ry = np.matrix([[math.cos(theta), 0., math.sin(theta)], [0., 1., 0.], [-1. * math.sin(theta), 0., math.cos(theta)]])
        #Rz = np.matrix([[math.cos(theta), -1. * math.sin(theta), 0.], [math.sin(theta), math.cos(theta), 0.], [0., 0., 1.]])

        # rotate
        coor_rotated = np.array(np.matrix(coorarray) * Rxyz)

        # tanslate back to original point (according to the 1st atom)
        vector = coorarray[0] - coor_rotated[0]
        coor_rotated = self.translate_coor(coor_rotated, distance=1.0, vector=vector, random=False)

        return coor_rotated

    def translate_coor(self, coorarray, distance=1.0, vector=(0.,0.,0.), random=True):
        '''
        translate the molecule in X, Y, Z axis by a rondom direction

        :param coorarray:

         array([[ 0.278,  1.922, -0.8  ],
                [ 0.701,  1.392, -1.44 ],
                [-0.444,  1.406, -0.585]])

        :return: perturbed coordinates array
        '''

        m, n = np.shape(coorarray)

        if random:
            # gene random axis
            x, y, z = np.array(self.make_rand_vector(dims=3)) * distance
        else:
            x, y, z = vector * distance

        Rxyz = np.matrix([
            [1., 0., 0., x],
            [0., 1., 0., y],
            [0., 0., 1., z],
            [0., 0., 0., 1.],
        ])

        # add one more dimension to x, y, z
        mol_new = np.concatenate((np.transpose(coorarray), np.array([[1., 1., 1.]])), axis=0)

        # translate
        coor_translated = np.transpose((Rxyz * mol_new)[:m, :])

        return np.array(coor_translated)

    def pdb_line2list(self, line_str):
        '''
        only work for pdb file generated by vmd!
        :param line_str:
        :return:
        '''

        idx = [0, 4, 11, 17, 22, 26, 38, 46, 54, 60, 66, 76, 78]
        line_list = []
        for j in range(len(idx) - 1):
            line_list.append(line_str[idx[j]: idx[j + 1]].strip())
        return line_list
    def pdb_list2line(self, line_list):
        '''
        only work for pdb file generated by vmd!
        :param line_list:
        :return:
        '''
        format_str = '%s%7s  %-4s%5s%4s    %8s%8s%8s%6s%6s      %-5s%s'
        line_str = format_str % tuple(line_list)
        return line_str

    def get_energies_nw(self, nwout):
        fcon = open(nwout).readlines()
        DFT = [i.strip().split()[-1] for i in fcon if i.find('Total DFT energy') > 0][-1]
        OneElec = [i.strip().split()[-1] for i in fcon if i.find('One electron energy') > 0][-1]
        Coulomb = [i.strip().split()[-1] for i in fcon if i.find('Coulomb energy') > 0][-1]
        Exchange = [i.strip().split()[-1] for i in fcon if i.find('Exchange energy') > 0][-1]
        Correlation = [i.strip().split()[-1] for i in fcon if i.find('Correlation energy') > 0][-1]
        NuRepulsion = [i.strip().split()[-1] for i in fcon if i.find('Nuclear repulsion energy') > 0][-1]
        NumInt = [i.strip().split()[-1] for i in fcon if i.find('Numeric. integr. density') > 0][-1]
        Atomic = [i.strip().split()[-1] for i in fcon if i.find('Sum of atomic energies:') > 0][-1]

        results = {'DFT': DFT, 'OneElec': OneElec, 'Coulomb': Coulomb, 'Exchange': Exchange, 'Correlation': Correlation,
                   'NuRepulsion': NuRepulsion, 'NumInt': NumInt, 'Atomic': Atomic}
        return results



    def perturb_xyz(self, xyz, xyz_out, offset=0.1):
        # exist or not
        if not os.path.exists(xyz):
            logging.error('ERROR: Cannot find xyz file at %s' % xyz)
            return 1

        # read coord
        mol = [i.strip().split() for i in open(xyz).readlines()]

        # all atom lines
        mol_atoms = np.array(mol[2:])

        # all coordinate
        mol_coors = mol_atoms[:, 1:4].astype(float)

        # perturbed coordinates
        mol_coors_off = self.perturb_coor(mol_coors, offset=offset)

        # writeback coordinates
        mol_atoms[:, 1:4] = mol_coors_off
        # form new molecule
        mol_new = mol[:2] + mol_atoms.tolist()

        # write the new mol to xyz file
        with open(xyz_out, 'w') as fout:
            for line in [' '.join(i) for i in mol_new]:
                fout.write('%s\n' % line)
            fout.close()

    def perturb_pdb(self, pdb, pdb_out, offset=0.1, decimal=3, rand_iterations=0, distance=1.0, perturb=True, rotate=False, translate=False):
        '''
        only work for pdb file generated by vmd!
        :param pdb:
        :param pdb_out:
        :param offset:
        :param decimal:
        :param rand_iterations: the more rand_iterations, the far from the input position
        :return:
        '''
        #print 'perturb_pdb(pdb=%s, pdb_out=%s, offset=%s, decimal=%s, rand_iterations=%s, distance=%s, perturb=%s, rotate=%s, translate=%s)' \
        #    % (pdb, pdb_out, str(offset), str(decimal), str(rand_iterations), str(distance), str(perturb), str(rotate), str(translate))

        # read header line
        mol_header_str = [i.strip() for i in open(pdb).readlines() if i.startswith('CRYST')]

        # read all ATOM lines
        mol_atoms_str = [i.strip() for i in open(pdb).readlines() if i.startswith('ATOM')]

        # split all ATOM lines
        mol_atoms = np.array([self.pdb_line2list(i) for i in mol_atoms_str])

        # all coordinates
        mol_coors = mol_atoms[:, 5:8].astype(float)


        # init mol_coors_off
        mol_coors_off = mol_coors

        # perturb each atom randomly
        if perturb in ['True', True]:
            # perturbed coordinates
            mol_coors_off = self.perturb_coor(mol_coors_off, offset=offset)
            counter = 0
            while counter < rand_iterations:
                mol_coors_off = self.perturb_coor(mol_coors_off, offset=offset)
                counter += 1

        # rotate the molecule randomly
        if rotate in ['True', True]:
            mol_coors_off = self.rotate_coor(mol_coors_off)

        # translate the molecule randomly
        if translate in ['True', True]:
            mol_coors_off = self.translate_coor(mol_coors_off, distance=distance)


        # writeback coordinates
        mol_atoms[:, 5:8] = mol_coors_off.round(decimal)

        # new atom str
        mol_atoms_str_updated = [self.pdb_list2line(i) for i in mol_atoms]

        # combine
        if len(mol_header_str) > 0:
            mol_combine = mol_header_str + mol_atoms_str_updated + ['END']
        else:
            mol_combine = mol_atoms_str_updated + ['END']

        # write out
        with open(pdb_out, 'w') as fout:
            for line in mol_combine:
                fout.write('%s\n' % line)
            fout.close()

        return mol_combine

    def perturb_pdb_2WAT(self, pdb, pdb_out, offset=0.1, decimal=3, rand_iterations=0, mindistance=2.0, maxdistance=12.0, perturb=True, rotate=False, translate=False):
        '''
        only work for pdb file generated by vmd!

        to separate two water molecules. Only operate the 2nd water molecule

        :param pdb:
        :param pdb_out:
        :param offset:
        :param decimal:
        :param rand_iterations: the more rand_iterations, the far from the input position
        :return:
        '''
        #print 'perturb_pdb_2WAT(pdb=%s, pdb_out=%s, offset=%s, decimal=%s, rand_iterations=%s, distance=%s, perturb=%s, rotate=%s, translate=%s)' \
        #    % (pdb, pdb_out, str(offset), str(decimal), str(rand_iterations), str(distance), str(perturb), str(rotate), str(translate))

        # read header line
        mol_header_str = [i.strip() for i in open(pdb).readlines() if i.startswith('CRYST')]

        # read all ATOM lines
        mol_atoms_str = [i.strip() for i in open(pdb).readlines() if i.startswith('ATOM')]

        # split all ATOM lines
        mol_atoms = np.array([self.pdb_line2list(i) for i in mol_atoms_str])

        # all coordinates
        mol_coors = mol_atoms[-3:, 5:8].astype(float)

        # init mol_coors_off
        mol_coors_off = mol_coors

        # perturb each atom randomly
        if perturb in ['True', True]:
            # perturbed coordinates
            mol_coors_off = self.perturb_coor(mol_coors_off, offset=offset)
            counter = 0
            while counter < rand_iterations:
                mol_coors_off = self.perturb_coor(mol_coors_off, offset=offset)
                counter += 1

        # rotate the molecule randomly
        if rotate in ['True', True]:
            mol_coors_off = self.rotate_coor(mol_coors_off)

        # translate the molecule randomly
        if translate in ['True', True]:
            mol_coors_off = self.translate_coor(mol_coors_off, distance=random.uniform(mindistance, maxdistance))


        # writeback coordinates
        mol_atoms[-3:, 5:8] = mol_coors_off.round(decimal)

        # new atom str
        mol_atoms_str_updated = [self.pdb_list2line(i) for i in mol_atoms]

        # combine
        if len(mol_header_str) > 0:
            mol_combine = mol_header_str + mol_atoms_str_updated + ['END']
        else:
            mol_combine = mol_atoms_str_updated + ['END']

        # write out
        with open(pdb_out, 'w') as fout:
            for line in mol_combine:
                fout.write('%s\n' % line)
            fout.close()

        return mol_combine

    def pdb2cmat(self, pdb, outName='', outFormat='vector', saveFile=True):
        '''
        convert a pdb file to Coulomb matrix file
        :param pdb: input pdb
        :param cmat: output coulomb matrix file
        :return:
        '''

        if not (pdb.endswith('.pdb') or pdb.endswith('.PDB')):
            pdb += '.pdb'

        try:
            fcon = open(pdb).readlines()
        except:
            logging.error('Cannot open file %s' % pdb)
            return ''

        # get file name
        fname = pdb[:-4]
        if outName == '':
            outName = '%s.cmat' % fname

        # format coordinates
        coor = [
                   [self.element_charge[j[2]], float(j[4]), float(j[5]), float(j[6])]
                   for j in
                   [i.strip().split() for i in fcon if i.startswith('ATOM')]
               ]

        # convert to matrix
        coormatrix = np.matrix(coor)

        # calculate coulomb matrix
        coulomb_matrix, coulomb_vector = self.coor2cmat(coormatrix)

        # write the result
        if saveFile:
            if outFormat in ['vector', 'v', 'Vector']:
                with open(outName, 'w') as fout:
                    ## insert '\n' every 10 element
                    #idx = range(0, len(coulomb_vector), 10)
                    #output_vector = np.insert(coulomb_vector.astype('str'), idx, np.array(['\n'] * len(idx)))
                    #fout.write(','.join(output_vector))

                    fout.write(',\n'.join(coulomb_vector.astype('str')))
                    fout.write('\n')
                fout.close()
            else:
                with open(outName, 'w') as fout:
                    coulomb_str = ['\n,'.join(i.astype('str')) for i in coulomb_matrix]
                    fout.writelines(coulomb_str)
                    fout.write('\n')
                fout.close()

        return coulomb_matrix, coulomb_vector

    def cal_energy_vmd(self, psf, pdb, outName='', saveFile=True):

        # prepare
        mol_name = pdb[:-4]
        if outName == '':
            outName = '%s.en' % mol_name

        self.tempname = '%s_%d' % (mol_name, random.randint(100000, 999999))

        # call vmd and run tcl command
        VMD.evaltcl('package require namdenergy')
        VMD.evaltcl('mol new      %s' % pdb)
        VMD.evaltcl('mol addfile  %s' % psf)
        VMD.evaltcl('set al [atomselect top %s]' % self.atomselection)
        Ens = VMD.evaltcl('namdenergy -all -sel $al -tempname %s -switch %s -cutoff %s -diel %s -updatesel'
                          % (str(self.tempname), str(self.switch), str(self.cutoff), str(self.dielectronic)))
        VMD.evaltcl('mol delete all')

        # format result
        labels = ['Frame', 'Time', 'Bond', 'Angle', 'Dihed', 'Impr', 'Elec', 'VdW', 'Conf', 'Nonbond', 'Total']
        Ens = Ens.replace('{', '').replace('}', '').split()

        df = pd.DataFrame([Ens], columns=labels)

        # output
        if saveFile:
            df.to_csv(outName)

        return df

    def measure_dist_vmd(self, psf, pdb, atomindex1, atomindex2):
        # call vmd and run tcl command
        VMD.evaltcl('mol new      %s' % pdb)
        VMD.evaltcl('mol addfile  %s' % psf)
        dist_str = VMD.evaltcl('measure bond {%d %d}' % (atomindex1, atomindex2))
        VMD.evaltcl('mol delete all')

        return float(dist_str)

    def extract_moiety_vmd(self, psf, pdb, atomselection, outpdb):
        # call vmd and run tcl command
        VMD.evaltcl('mol new      %s' % pdb)
        VMD.evaltcl('mol addfile  %s' % psf)
        VMD.evaltcl('set al [atomselect top "%s"]' % atomselection)
        VMD.evaltcl('$al writepdb %s' % outpdb)
        VMD.evaltcl('mol delete all')

        return

    def opt_pdb_ob(self, pdb, outpdb, ff='UFF', method='cg', nstep=100):

        # read header line
        mol_header_str = [i.strip() for i in open(pdb).readlines() if i.startswith('CRYST')]
        # read all ATOM lines
        mol_atoms_str = [i.strip() for i in open(pdb).readlines() if i.startswith('ATOM')]
        # split all ATOM lines
        mol_atoms = np.array([self.pdb_line2list(i) for i in mol_atoms_str])
        # all coordinates
        mol_coors = mol_atoms[:, 5:8].astype(float)


        #### optimization
        # command
        cmd = '%s -ff %s -n %d -%s %s' % (self.obminimize, ff, nstep, method, pdb)
        # optimize
        optcmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        # grab output and error
        optout, opterr = optcmd.communicate()
        atomlines = [i for i in optout.splitlines() if i.startswith('ATOM')]
        mol_coors_opted = np.array([self.pdb_line2list(i) for i in atomlines])[:, 5:8].astype(float)


        # writeback coordinates
        mol_atoms[:, 5:8] = mol_coors_opted

        # new atom str
        mol_atoms_str_updated = [self.pdb_list2line(i) for i in mol_atoms]

        # combine
        if len(mol_header_str) > 0:
            mol_combine = mol_header_str + mol_atoms_str_updated + ['END']
        else:
            mol_combine = mol_atoms_str_updated + ['END']

        # write out
        with open(outpdb, 'w') as fout:
            for line in mol_combine:
                fout.write('%s\n' % line)
            fout.close()

        return mol_combine

    def opt_pdb_namd(self, psf, pdb, outpdb, nstep=100):
        template = '''
structure          %s
coordinates        %s

# Force-Field Parameters
exclude             1-4
1-4scaling          1.0
cutoff              12.
switching           on
switchdist          10.
pairlistdist        13.5


timestep            1
rigidBonds          water
useSettle           on
nonbondedFreq       1
fullElectFrequency  1
stepsPerCycle       1


temperature         100
outputname          tmpforopt
binaryoutput        no

# Input
paraTypeCharmm      on
parameters          %s

minimize %d

        ''' % (psf, pdb, self.namd2ff, int(nstep))


        # write the tmplate
        open('tmpforopt.namd', 'w').write(template)

        #### optimization
        # command
        cmd = '%s tmpforopt.namd' % self.namd2
        # optimize
        optcmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        # grab output and error
        optout, opterr = optcmd.communicate()

        try:
            shutil.move('tmpforopt.coor', outpdb)
        except:
            logging.info('Cannot find tmpforopt.coor!')
            pass

        try:
            os.remove('tmpforopt.vel')
        except:
            pass

        try:
            os.remove('tmpforopt.xsc')
        except:
            pass

        try:
            os.remove('tmpforopt.namd')
        except:
            pass

        return

    def cal_energy_NW(self, pdb, nwinput, xc='b3lyp', basis='6-31+g**', charge=0, mult=1):
        '''
        generate the input file for cal energy in NW

        '''
        # command
        cmd = '%s -ipdb %s -oxyz' % (self.obabel, pdb)

        # convert
        optcmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        # grab output and error
        optout, opterr = optcmd.communicate()

        # prepare mol_coors_str
        mol_coors = [i for i in optout.split('\n')[2:] if len(i)>0]
        mol_coors_str = '\n'.join(mol_coors)


        # nw template
        nw_template = '''start %s
title "%s"

scratch_dir .
memory 8 GB

charge %d

geometry units angstroms noautoz
 symmetry c1
%s
end

basis
 * library %s
end

dft
 decomp
 xc %s
 mult %d
end

task dft energy

        ''' % (os.path.basename(pdb), os.path.basename(pdb), charge, mol_coors_str, basis, xc, mult)


        # write out
        with open(nwinput, 'w') as fout:
            fout.writelines(nw_template)
            #for line in nw_template:
            #    fout.write('%s\n' % line)
            fout.close()

        return nw_template






    def gen_random_coor_pdbs(self, inPDB, numCopys, offset, outDir='perturb',
                             distance=1.0, perturb=True, rotate=False, translate=False, separate2WAT=False):
        # exist or not
        if not os.path.exists(inPDB):
            logging.error('ERROR: Cannot find pdb file at %s' % inPDB)
            return 1

        # prepare outDir
        if os.path.exists(outDir):
            try:
                os.removedirs(outDir)
            except:
                pass

        try:
            os.makedirs(outDir)
        except:
            pass

        mol_name = inPDB[:-4]

        counter = 0
        while counter < numCopys:
            mol_name_output = '%s/%d.pdb' % (outDir, counter)

            if separate2WAT in ['True', True]:
                self.perturb_pdb_2WAT(inPDB, mol_name_output, offset=offset, maxdistance=distance, perturb=perturb, rotate=rotate, translate=translate)
            else:
                self.perturb_pdb(inPDB, mol_name_output, offset=offset, distance=distance, perturb=perturb, rotate=rotate, translate=translate)

            counter += 1

    def convert_pdb2cmat(self, inDir, outDir='cmat', outFormat='vector'):
        # exist or not
        if not os.path.exists(inDir):
            logging.error('ERROR: Cannot find pdb directory %s' % inDir)
            return 1

        # prepare outDir
        if os.path.exists(outDir):
            try:
                os.removedirs(outDir)
            except:
                pass

        try:
            os.makedirs(outDir)
        except:
            pass

        # all pdb files
        pdb_files = sorted(os.listdir(inDir), key=lambda x: int(x.split('.')[0]))

        # output file name
        outFile = '%s.cmat' % os.path.basename(inDir)

        # convert and write
        for pdb in pdb_files:
            mol_name = pdb[:-4]
            cmat_m, cmat_v = self.pdb2cmat('%s/%s' % (inDir, pdb), outName='%s/%s.cmat' % (outDir, mol_name), outFormat=outFormat)

    def cal_energy_mm(self, inDir, psf, outDir='mmEn'):
        # exist or not
        if not os.path.exists(inDir):
            logging.error('ERROR: Cannot find pdb directory %s' % inDir)
            return 1

        # prepare outDir
        if os.path.exists(outDir):
            try:
                os.removedirs(outDir)
            except:
                pass
        try:
            os.makedirs(outDir)
        except:
            pass

        # all pdb files
        pdb_files = sorted(os.listdir(inDir), key=lambda x: int(x.split('.')[0]))

        # convert and write
        for pdb in pdb_files:
            mol_name = pdb[:-4]
            self.cal_energy_vmd(psf, '%s/%s' % (inDir, pdb), outName='%s/%s.en' % (outDir, mol_name))

    def cal_energy_dft(self, inDir, outDir='dftEn', xc='b3lyp', basis='6-31+g**', charge=0, mult=1):
        # exist or not
        if not os.path.exists(inDir):
            logging.error('ERROR: Cannot find pdb directory %s' % inDir)
            return 1

        # prepare outDir
        if os.path.exists(outDir):
            try:
                os.removedirs(outDir)
            except:
                pass
        try:
            os.makedirs(outDir)
        except:
            pass

        # all pdb files
        pdb_files = sorted(os.listdir(inDir), key=lambda x: int(x.split('.')[0]))

        # convert and write
        for pdb in pdb_files:
            mol_name = pdb[:-4]
            self.cal_energy_NW('%s/%s' % (inDir, pdb), nwinput='%s/%s.nw' % (outDir, mol_name), xc=xc, basis=basis, charge=charge, mult=mult)

    def extract_allmoiety(self, inDir, psf, atomselection, outDir='moiety'):
        # exist or not
        if not os.path.exists(inDir):
            logging.error('ERROR: Cannot find pdb directory %s' % inDir)
            return 1

        # prepare outDir
        if os.path.exists(outDir):
            try:
                os.removedirs(outDir)
            except:
                pass
        try:
            os.makedirs(outDir)
        except:
            pass

        # all pdb files
        pdb_files = sorted(os.listdir(inDir), key=lambda x: int(x.split('.')[0]))

        # convert and write
        for pdb in pdb_files:
            self.extract_moiety_vmd(psf, '%s/%s' % (inDir, pdb), atomselection, outpdb='%s/%s' % (outDir, pdb))


    def opt_allpdbs_ob(self, inDir, outDir='opted', ff='UFF', method='cg', nstep=10):
        # exist or not
        if not os.path.exists(inDir):
            logging.error('ERROR: Cannot find pdb directory %s' % inDir)
            return 1

        # prepare outDir
        if os.path.exists(outDir):
            try:
                os.removedirs(outDir)
            except:
                pass
        try:
            os.makedirs(outDir)
        except:
            pass

        # all pdb files
        pdb_files = sorted(os.listdir(inDir), key=lambda x: int(x.split('.')[0]))

        # convert and write
        for pdb in pdb_files:
            self.opt_pdb_ob('%s/%s' % (inDir, pdb), outpdb='%s/%s' % (outDir, pdb), ff=ff, method=method, nstep=nstep)

    def opt_allpdbs_namd(self, psf, inDir, outDir='opted', nstep=10):
        # exist or not
        if not os.path.exists(inDir):
            logging.error('ERROR: Cannot find pdb directory %s' % inDir)
            return 1

        # prepare outDir
        if os.path.exists(outDir):
            try:
                os.removedirs(outDir)
            except:
                pass
        try:
            os.makedirs(outDir)
        except:
            pass

        # all pdb files
        pdb_files = sorted(os.listdir(inDir), key=lambda x: int(x.split('.')[0]))

        # convert and write
        for pdb in pdb_files:
            self.opt_pdb_namd(psf, '%s/%s' % (inDir, pdb), outpdb='%s/%s' % (outDir, pdb), nstep=nstep)

    def measure_alldist_vmd(self, psf, inDir, atomindex1, atomindex2, outcsv='mmEn.csv', rewrite=False):
        # exist or not
        if not os.path.exists(inDir):
            logging.error('ERROR: Cannot find pdb directory %s' % inDir)
            return 1

        # all pdb files
        pdb_files = sorted(os.listdir(inDir), key=lambda x: int(x.split('.')[0]))

        alldist = []
        # convert and write
        for pdb in pdb_files:
            try:
                d = self.measure_dist_vmd(psf, '%s/%s' % (inDir, pdb), atomindex1, atomindex2)
            except:
                d = 0
            alldist.append(d)

        # save the result
        if os.path.exists(outcsv) and not rewrite:
            df = pd.DataFrame.from_csv(outcsv)
            df['DistOO'] = alldist
            df.to_csv(outcsv)
        else:
            df_alldist = pd.DataFrame(alldist)
            df_alldist.to_csv(outcsv)

        return alldist

    def collect_2csv(self, inDir, type='en'):
        if type in ['en', 'Energy', 'energy', 'En']:
            # get file list
            files = sorted([i for i in os.listdir(inDir) if i.endswith('.en')], key=lambda x: int(x.split('.en')[0]))
            # read all files into dataframe
            df = pd.DataFrame()
            for f in files:
                df = df.append(pd.DataFrame.from_csv('%s/%s' % (inDir, f)))

        elif type in ['cmat', 'Cmat']:
            # get file list
            files = sorted([i for i in os.listdir(inDir) if i.endswith('.cmat')], key=lambda x: int(x.split('.cmat')[0]))
            # read all files into dataframe
            df = pd.DataFrame()
            for f in files:
                fcon = open('%s/%s' % (inDir, f)).readlines()
                df_tmp = pd.DataFrame([[i.strip(',\n') for i in fcon]])
                df = df.append(df_tmp)

        elif type in ['nwout', 'NWout']:
            # get file list
            files = sorted([i for i in os.listdir(inDir) if i.endswith('.out')], key=lambda x: int(x.split('.')[0]))
            # read all files into dataframe
            df = pd.DataFrame()
            for f in files:
                try:
                    energies = self.get_energies_nw('%s/%s' % (inDir, f))
                except:
                    logging.error('\n\nFailed to get energies from %s/%s' % (inDir, f))
                    energies = {}

                df_tmp = pd.DataFrame([energies])
                df = df.append(df_tmp)

        # format
        df.index = xrange(len(df))

        # save the result to csv
        if inDir.endswith('/'):
            outName = '/'.join(inDir.split('/')[:-1])
        else:
            outName = inDir

        df.to_csv('%s.csv' % outName)

        return df


if __name__ == '__main__':
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Perform small molecule manipulation.')
    parser.add_argument("-ipdb", nargs='?', type=str, help="Input pdb file")
    parser.add_argument("-ipsf", nargs='?', type=str, help="Input psf file")
    parser.add_argument("-ocmat", nargs='?', type=str, default="", help="Output cmat file")
    parser.add_argument("-ncopy", nargs='?', type=int, default=10, help="Number of copies to generate. Default: 10")
    parser.add_argument("-offset", nargs='?', type=float, default=0.1, help="Offset value on each coordinate of the atoms in Angstrom. Default: 0.1")
    parser.add_argument("-distance", nargs='?', type=float, default=1.0, help="Distance/Max Distance for translation. Default: 1.0")
    parser.add_argument("-nstep", nargs='?', type=int, default=10, help="Number of steps for optimization Default: 10")
    parser.add_argument("-atomselection", nargs='?', type=str, default='all', help="Atomselection used by extrating moiety Default: all")
    parser.add_argument("-atomindex1", nargs='?', type=int, default=0, help="Atomsindex1 Default: 0")
    parser.add_argument("-atomindex2", nargs='?', type=int, default=3, help="Atomsindex2 Default: 3")
    parser.add_argument("-inDir", nargs='?', type=str, default='perturb', help="input dir wich contains all pdb files. Default: 'perturb'")
    parser.add_argument("-outDir", nargs='?', type=str, default='perturb', help="Output dir for coordinates perturbation. Default: 'perturb'")
    parser.add_argument("-outcsv", nargs='?', type=str, default='out.csv', help="Csv file for saveing some results. Default: 'out.csv'")

    parser.add_argument("-perturbation", nargs='?', type=str, default=True, help="Perform perturbation. Default: True")
    parser.add_argument("-rotation", nargs='?', type=str, default=False, help="Perform rotation. Default: False")
    parser.add_argument("-translation", nargs='?', type=str, default=False, help="Perform translation. Default: False")
    parser.add_argument("-sep2wat", nargs='?', type=str, default=False, help="Separate 2 water molecules")
    parser.add_argument("-extract", nargs='?', type=str, default=False, help="Extract a moiety of the input molecule")
    parser.add_argument("-convtocmat", nargs='?', help="Convert pdb files in inDir into cmat files")
    parser.add_argument("-optpdb", nargs='?', help="Optimize pdb files in inDir")
    parser.add_argument("-calDist", nargs='?', help="Calculate distance between atomindex1 and atomindex2 for pdb files in inDir with VMD.")
    parser.add_argument("-calEn", nargs='?', help="Calculate MM Energy of pdb files in inDir with VMD.")
    parser.add_argument("-calEnDFT", nargs='?', help="Generate NWchem input files for pdb files in inDir.")
    parser.add_argument("-getEnDFT", nargs='?', help="Collect DFT Energy from NWchem output files in inDir.")
    parser.add_argument("perturb", nargs='?', help="Perform perturbation on the coordinates the atoms")



    # all arguments
    args = parser.parse_args()
    sl = smallLigand()

    Pwd = os.path.abspath(os.curdir)

    print args

    if args.ipdb and args.ocmat:
        print 'Converting %s to %s ...' % (args.ipdb, args.ocmat)
        sl.pdb2cmat(args.ipdb, args.ocmat, outFormat='vector')

    if args.ipdb and args.perturb:
        print 'Performing perturbation on %s, with numCopys=%d, offset=%f, outDir=%s, ' \
              'distance=%s, perturb=%s, rotate=%s, translate=%s, separate2WAT=%s' \
              % (args.ipdb, args.ncopy, args.offset, args.outDir, str(args.distance), args.perturbation, args.rotation, args.translation, args.sep2wat)
        sl.gen_random_coor_pdbs(args.ipdb, args.ncopy, args.offset, args.outDir,
                distance=args.distance, perturb=args.perturbation, rotate=args.rotation, translate=args.translation, separate2WAT=args.sep2wat)

    if args.convtocmat and args.inDir and args.outDir:
        print 'Converting pdb files in %s to cmat files ...' % args.inDir
        sl.convert_pdb2cmat(args.inDir, outDir=args.outDir, outFormat='vector')
        sl.collect_2csv(inDir=args.outDir, type='cmat')

    if args.extract and args.atomselection and args.ipsf and args.inDir and args.outDir:
        print 'Extracting %s from files in %s' % (args.atomselection, args.inDir)
        sl.extract_allmoiety(inDir=args.inDir, psf=args.ipsf, atomselection=args.atomselection, outDir=args.outDir)

    if args.calEn and args.ipsf and args.inDir and args.outDir:
        print 'Calculating MM Energy of pdb files in %s with VMD ...' % args.inDir
        sl.cal_energy_mm('%s/%s' % (Pwd, args.inDir), '%s/%s' % (Pwd, args.ipsf), outDir=args.outDir)
        sl.collect_2csv(inDir=args.outDir, type='en')

    if args.calDist and args.ipsf and args.inDir and args.outcsv:
        print 'Measure the distance between %d and %d for pdb files in %s with VMD ...' % (args.atomindex1, args.atomindex2, args.inDir)
        sl.measure_alldist_vmd(args.ipsf, args.inDir, args.atomindex1, args.atomindex2, outcsv=args.outcsv, rewrite=False)

    if args.calEnDFT and args.inDir and args.outDir:
        print 'Generating input file for nwchem with pdb form %s ...' % args.inDir
        sl.cal_energy_dft(args.inDir, args.outDir, xc='b3lyp', basis='6-31+g**', charge=0, mult=1)

    if args.getEnDFT and args.inDir:
        print 'Collect DFT Energy from NWchem output files form %s ...' % args.inDir
        sl.collect_2csv(args.inDir, type='nwout')

    if args.optpdb and args.outDir and args.inDir and args.ipsf:
        print 'Optimize pdb files in %s ...' % args.inDir
        sl.opt_allpdbs_namd(args.ipsf, inDir=args.inDir, outDir=args.outDir, nstep=args.nstep)




    '''
    #### examples ####

    # separate the 2nd water molecule
    ~/tools/myPythonLib/PLg09/smallLigand.py -ipdb solvate.pdb -ncopy 10000 -offset 1.0 -outDir perturb -distance 20.0\
     -sep2wat True -perturbation False -rotation True -translation True perturb

    # optimize 10 steps
    ~/tools/myPythonLib/PLg09/smallLigand.py -ipsf solvate.psf -inDir perturb -outDir perturb_opted -optpdb True -nstep 10

    # convert to cmat
    ~/tools/myPythonLib/PLg09/smallLigand.py -inDir perturb_opted -outDir cmat_opted -convtocmat True

    # calculate the energies
    ~/tools/myPythonLib/PLg09/smallLigand.py -ipsf solvate.psf -inDir perturb_opted -outDir mmEn_opted -calEn True

    # extract one water
    ~/tools/myPythonLib/PLg09/smallLigand.py -ipsf solvate.psf -inDir perturb_opted -outDir wat1 -atomselection 'resid 1301' -extract True

    ~/tools/myPythonLib/PLg09/smallLigand.py -ipsf solvate.psf -inDir perturb_opted -outDir wat2 -atomselection 'resid 8425' -extract True

    # measure the oo distance
    ~/tools/myPythonLib/PLg09/smallLigand.py -ipsf solvate.psf -inDir perturb_opted -atomindex1 0 -atomindex2 3 -outcsv mmEn_opted.csv -calDist True

    # gen NWchem input
    ~/tools/myPythonLib/PLg09/smallLigand.py -inDir perturb_opted -outDir dftEn_opted -calEnDFT True

    # collect the result from DFT
    ~/tools/myPythonLib/PLg09/smallLigand.py -inDir dftEn_opted_out -getEnDFT True


    '''