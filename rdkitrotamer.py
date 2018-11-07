#!/home/p6n/anaconda2/envs/rdkit/bin/python
#
# @purpose
#   to replace the obrotamer which is not good enough for sampling single bond
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Mar 15 2017
#
import subprocess, logging, argparse, os, random

from rdkit import Chem
from rdkit.Chem import AllChem

class plrdkit():
    def __init__(self):
        self.obabel = '/home/p6n/anaconda2/bin/obabel'
        return

    def readinMol(self, molFile, type='mol2', sanitize=True, removeHs=True):
        if type in ['mol2', 'Mol2']:
            m = Chem.MolFromMol2File(molFile, sanitize=sanitize, removeHs=removeHs)
        elif type in ['pdb', 'PDB']:
            m = Chem.MolFromPDBFile(molFile, sanitize=sanitize, removeHs=removeHs)
        else:
            logging.error('Cannot read %s, only mol2/pdb accepted' % molFile)
        return m


    def genConformations(self, m, n):
        #m = Chem.AddHs(m)
        ids = AllChem.EmbedMultipleConfs(m, numConfs=n)
        for id in ids:
            AllChem.UFFOptimizeMolecule(m, confId=id, maxIters=0, vdwThresh=10.0)
        # EmbedMultipleConfs returns a Boost-wrapped type which
        # cannot be pickled. Convert it to a Python list, which can.
        return m, list(ids)

    def writetoSDF(self, m, SDFfile):
        Nconfs = len(m.GetConformers())

        writer = Chem.SDWriter(SDFfile)

        for i in range(Nconfs):
            writer.write(m, i)

        writer.close()

    def SDFtoMOL2(self, SDFfile, MOL2file):
        # convert sdf to mol2
        cmd = '%s -isdf %s -O %s' % (self.obabel, SDFfile, MOL2file)
        cmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = cmd.communicate()


if __name__ == '__main__':
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Perform small molecule manipulation.')
    parser.add_argument("-ipdb", nargs='?', type=str, help="Input pdb file")
    parser.add_argument("-imol2", nargs='?', type=str, help="Input mol2 file")
    parser.add_argument("-omol2", nargs='?', type=str, help="Output mol2 file")
    parser.add_argument("-N", nargs='?', type=int, default=1000, help="Number of rotamers to generate")
    parser.add_argument("genRotamers", nargs='?', help="Perform perturbation on the coordinates the atoms")

    # all arguments
    args = parser.parse_args()
    pl = plrdkit()

    # perform conversion
    if args.ipdb and args.omol2 and args.N and args.genRotamers:
        if os.path.exists(args.ipdb):
            mol = args.ipdb
            mol_type = 'pdb'
        else:
            logging.error('Could not find input file: %s' % args.ipdb)
            exit()
    elif args.imol2 and args.omol2 and args.N and args.genRotamers:
        if os.path.exists(args.imol2):
            mol = args.imol2
            mol_type = 'mol2'
        else:
            logging.error('Could not find input file: %s' % args.imol2)
            exit()


    # read the mol file
    m = pl.readinMol(mol, type=mol_type, sanitize=True, removeHs=False)
    # gen conformations
    m, ids = pl.genConformations(m, args.N)
    # write to the tmp sdf file
    tmp_sdf = '/tmp/tmp_plrdkit_%d.sdf' %  random.randint(100000, 999999)
    pl.writetoSDF(m, tmp_sdf)
    # convert sdf to mol2
    pl.SDFtoMOL2(tmp_sdf, args.omol2)
