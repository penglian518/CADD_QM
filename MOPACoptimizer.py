#!/home/p6n/anaconda2/bin/python
####!/usr/bin/env python
import argparse, subprocess, os

from ase.build import molecule
from ase.optimize import BFGS
from ase.io import read
from ase.calculators.mopac import MOPAC

class MOPACCalculator:
    def __init__(self):
        self.mopac = '/opt/mopac/MOPAC2016.exe'
        self.mopaclib = '/opt/mopac'
        self.task = '1SCF GRADIENTS'
        self.processors = 1
        self.charge = 0
        self.uks = False  # unrestricted Kohn-Sham, required if there is an odd number of electrons
        self.method = 'PM7'
        self.eps = 78.4

        self.fmax = 0.05
        self.steps = 200
        self.trj = 'mopac.trj' # or None
        self.rst = 'mopac.rst'  # or None
        self.log = '-'  # '-' for stdout

    def geo_opt(self, inputXYZ):
        # name
        name = inputXYZ[:-4]
        self.trj = '%s.trj' % name
        self.rst = '%s.rst' % name
        self.log = '%s.log' % name
        opted_xyz = '%s_opted.xyz' % name # requires to have trj

        # read the coordinates
        atoms = read(inputXYZ)

        # configuration
        os.environ['ASE_MOPAC_COMMAND'] = '%s PREFIX.mop 2> /dev/null' % self.mopac
        os.environ['LD_LIBRARY_PATH'] = '%s:$LD_LIBRARY_PATH' % self.mopaclib
        calc = MOPAC(label=name, uhf=self.uks, method=self.method, task='%s THREADS=%s CHARGE=%s EPS=%s' % (self.task, self.processors, self.charge, self.eps))
        #calc = MOPAC()
        #print(os.environ)
        #print(calc.parameters)

        atoms.set_calculator(calc)

        # optimization
        dyn = BFGS(atoms, trajectory=self.trj, restart=self.rst, logfile=self.log)
        dyn.run(fmax=self.fmax, steps=self.steps)

        # convert the last frame to xyz
        frm = read(self.trj, -1)
        frm.write(opted_xyz, 'xyz')


if __name__ == '__main__':
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Perform geometry optimization with MOPAC through ASE.')
    parser.add_argument("-np", nargs='?', type=int, default=1, help="Number of processors to use")
    parser.add_argument("-method", nargs='?', type=str, default="PM7", help="Functional to use")
    parser.add_argument("-steps", nargs='?', type=int, default=300, help="Max minimization steps for each rotamer")
    parser.add_argument("-charge", nargs='?', type=int, default=0, help="Charge of the molecule")
    parser.add_argument("-openshell", nargs='?', type=str, default='False', help="Open shell or not?")
    parser.add_argument("-fmax", nargs='?', type=float, default=0.01, help="Max force for convergence")
    parser.add_argument("-eps", nargs='?', type=float, default=78.4, help="dielectric constant for the solvent.")

    parser.add_argument("batch", nargs='?', help="Batch mode will require structures in one mol2 file, "
                                                 "and will output optimized structures into one pdb file.")
    parser.add_argument("mol", type=str, default="", help="Input the molecule file (with file type e.g. 'mol.xyz')")


    obabel = '/home/p6n/anaconda2/bin/obabel'
    # all arguments
    args = parser.parse_args()

    # path
    moldir = os.path.dirname(args.mol)
    molbasename = os.path.basename(args.mol)
    if len(moldir) > 0:
        os.chdir(moldir)

    molname = '.'.join(molbasename.split('.')[:-1])
    moltype = molbasename.split('.')[-1]

    # apply the arguments to the calculator
    cc = MOPACCalculator()
    cc.processors = args.np
    cc.method = args.method
    cc.fmax = args.fmax
    cc.steps = args.steps
    cc.charge = args.charge
    cc.eps = args.eps
    if args.openshell in ['False', 'false', 'F', 'f', '0', 0]:
        cc.uks = False
    elif args.openshell in ['True', 'true', 'T', 't', '1', 1]:
        cc.uks = True


    # batch mode or not
    if args.batch:
        # output pdb file
        if molname.endswith('.rot'):
            outPDB = '%s.opted.pdb' % '.'.join(molname.split('.')[:-1])
            outEN = '%s.opted.en' % '.'.join(molname.split('.')[:-1])
        else:
            outPDB = '%s.opted.pdb' % molname
            outEN = '%s.opted.en' % molname
        # clean for output
        try:
            os.remove(outPDB)
        except:
            pass
        try:
            os.remove(outEN)
        except:
            pass

        # convert structures in the mol2 file into separated xyz files
        cmd = '%s -i%s %s.%s -O %s_.xyz -m' % (obabel, moltype, molname, moltype, molname)
        cmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = cmd.communicate()

        sepXYZs = [i for i in os.listdir('.') if i.startswith('%s_' % molname)]
        sortedXYZs = sorted(sepXYZs, key=lambda x: int(x.split('_')[-1][:-4]))

        outEN_list = []
        # run the minimization
        for xyz in sortedXYZs:
            molname_sep = xyz[:-4]

            cc.rst = '%s.rst' % molname_sep
            cc.trj = '%s.trj' % molname_sep
            cc.log = '%s.log' % molname_sep
            cc.out = '%s.out' % molname_sep
            cc.mop = '%s.mop' % molname_sep
            cc.arc = '%s.arc' % molname_sep

            # delete the rst trj log files
            try:
                os.remove(cc.log)
            except:
                pass
            try:
                os.remove(cc.trj)
            except:
                pass
            try:
                os.remove(cc.rst)
            except:
                pass

            # run the minimization
            try:
                cc.geo_opt('%s.xyz' % molname_sep)
            except Exception as e:
                print('Cannot optimize structure %s, the error is:\n%s' % (xyz, e))

            # convert optimized structure to one pdb file
            try:
                # convert structures in the mol2 file into separated xyz files
                cmd = '%s -ixyz %s_opted.xyz -O %s_opted.pdb' % (obabel, molname_sep, molname_sep)
                cmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out, err = cmd.communicate()
            except:
                print('Failed to convert molecule %s_opted.xyz to pdb' % molname_sep)
            try:
                # convert structures in the mol2 file into separated xyz files
                cmd = 'cat %s_opted.pdb >> %s' % (molname_sep, outPDB)
                cmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out, err = cmd.communicate()
            except:
                print('Failed to combine molecule %s_opted.pdb to output pdb file' % molname_sep)

            # grab the energy
            try:
                en = open(cc.log).readlines()[-1].strip().split()[3]
            except:
                en = 0
                pass
            outEN_list.append('%s\n' % str(en))
            # delete the xyz and trj
            try:
                os.remove('%s.xyz' % molname_sep)
            except:
                pass
            try:
                os.remove('%s_opted.xyz' % molname_sep)
            except:
                pass
            try:
                os.remove('%s_opted.pdb' % molname_sep)
            except:
                pass
            try:
                os.remove(cc.log)
            except:
                pass
            try:
                os.remove(cc.trj)
            except:
                pass
            try:
                os.remove(cc.rst)
            except:
                pass
            try:
                os.remove(cc.mop)
            except:
                pass
#            try:
#                os.remove(cc.out)
#            except:
#                pass
            try:
                os.remove(cc.arc)
            except:
                pass

        # write the energy list
        open(outEN, 'w').writelines(outEN_list)

    else:
        # if not in xyz convert it to xyz
        if moltype not in ['xyz']:
            if moltype in ['smi']:
                cmd = '%s -i%s %s.%s -O %s.xyz --gen3D' % (obabel, moltype, molname, moltype, molname)
            else:
                # convert .pdb to .xyz
                cmd = '%s -i%s %s.%s -O %s.xyz' % (obabel, moltype, molname, moltype, molname)
            cmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out, err = cmd.communicate()

        cc.rst = '%s.rst' % molname
        cc.trj = '%s.trj' % molname
        cc.log = '%s.log' % molname
        cc.out = '%s.out' % molname
        cc.mop = '%s.mop' % molname
        cc.arc = '%s.arc' % molname

        # delete the rst trj log files
        try:
            os.remove(cc.log)
        except:
            pass
        try:
            os.remove(cc.trj)
        except:
            pass
        try:
            os.remove(cc.rst)
        except:
            pass

        # run the minimization
        try:
            cc.geo_opt('%s.xyz' % molname)
        except Exception as e:
            print('Cannot optimize structure %s.xyz' % molname)
            print('The error is:\n%s' % e)

        # delete the rst trj log files
        try:
            os.remove(cc.mop)
        except:
            pass
        try:
            os.remove(cc.trj)
        except:
            pass
        try:
            os.remove(cc.rst)
        except:
            pass
        try:
            os.remove(cc.out)
        except:
            pass
        try:
            os.remove(cc.arc)
        except:
            pass

