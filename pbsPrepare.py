#!/usr/bin/env python
#
# @purpose
#   to prepare PBS submit scripts, according to the conf dict
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 23 2016
#


def genPBS_Condo_g09(conf):

    # check and set default
    if 'PBS_queue' not in conf.keys():
        conf['PBS_queue'] = 'batch'
    if 'PBS_nodes' not in conf.keys():
        conf['PBS_nodes'] = '1'
    if 'PBS_ppn' not in conf.keys():
        conf['PBS_ppn'] = '32'
    if 'PBS_walltime' not in conf.keys():
        conf['PBS_walltime'] = '24:00:00'
    if 'PBS_qos' not in conf.keys():
        conf['PBS_qos'] = 'burst'
    if 'PBS_jobname' not in conf.keys():
        print 'ERROR: PBS_jobname is required by PBS!'
        return 1
    if 'PBS_molname' not in conf.keys():
        print 'ERROR: PBS_molname is required by PBS!'
        return 1
    if 'PBS_molname_of_previous_step' not in conf.keys():
        print 'ERROR: PBS_molname_of_previous_step is required by PBS!'
        return 1


    condo_template = '''#PBS -q %s
#PBS -l nodes=%s:ppn=%s,walltime=%s
#PBS -l qos=%s
#PBS -N %s
#PBS -W group_list=cades-user
#PBS -A bsd-burst
#PBS -j oe

mol="%s"

dir=$PBS_O_WORKDIR

module load /software/user_tools/current/modules/cades-bsd/gaussian/g09 env/cades-bsd

#userscratch=$localscratch/$USER
userscratch=/lustre/or-hydra/cades-bsd/p6n/tmp

mkdir -p $userscratch/$PBS_JOBID
export GAUSS_SCRDIR=$userscratch/$PBS_JOBID
cd $GAUSS_SCRDIR

cp $dir/%s.chk .
g09 < $dir/$mol.com > $dir/$mol.log
cp -p $mol.chk $dir/

rm -fr $userscratch/$PBS_JOBID

''' % (conf['PBS_queue'], conf['PBS_nodes'], conf['PBS_ppn'],
       conf['PBS_walltime'], conf['PBS_qos'], conf['PBS_jobname'],
       conf['PBS_molname'], conf['PBS_molname_of_previous_step']
       )

    return condo_template

def genPBS_Condo_nw(conf):

    # check and set default
    if 'PBS_queue' not in conf.keys():
        conf['PBS_queue'] = 'batch'
    if 'PBS_nodes' not in conf.keys():
        conf['PBS_nodes'] = '1'
    if 'PBS_ppn' not in conf.keys():
        conf['PBS_ppn'] = '32'
    if 'PBS_walltime' not in conf.keys():
        conf['PBS_walltime'] = '24:00:00'
    if 'PBS_qos' not in conf.keys():
        conf['PBS_qos'] = 'burst'
    if 'PBS_jobname' not in conf.keys():
        print 'ERROR: PBS_jobname is required by PBS!'
        return 1
    if 'PBS_molname' not in conf.keys():
        print 'ERROR: PBS_molname is required by PBS!'
        return 1
    if 'PBS_molname_of_previous_step' not in conf.keys():
        print 'ERROR: PBS_molname_of_previous_step is required by PBS!'
        return 1


    condo_template = '''#PBS -q %s
#PBS -l nodes=%s:ppn=%s,walltime=%s
#PBS -l qos=%s
#PBS -N %s
#PBS -W group_list=cades-user
#PBS -A bsd-burst
#PBS -j oe


mol="%s"

dir=$PBS_O_WORKDIR

module load /software/user_tools/current/modules/cades-bsd/gaussian/g09 env/cades-bsd
module load PE-gnu
module load nwchem/6.6

#userscratch=$localscratch/$USER
userscratch=/lustre/or-hydra/cades-bsd/p6n/tmp

export GAUSS_SCRDIR=$userscratch/$PBS_JOBID/$$
mkdir -p $GAUSS_SCRDIR

cd $GAUSS_SCRDIR

cp $dir/* .
mkdir tmp
mpirun -np $PBS_NP nwchem $mol.nw > $dir/$mol.out 2>&1
rm -fr tmp
cp -u ./* $dir

rm -fr $GAUSS_SCRDIR

''' % (conf['PBS_queue'], conf['PBS_nodes'], conf['PBS_ppn'],
       conf['PBS_walltime'], conf['PBS_qos'], conf['PBS_jobname'],
       conf['PBS_molname']
       )

    return condo_template

def genPBS_Condo_nw_complied(conf):

    # check and set default
    if 'PBS_queue' not in conf.keys():
        conf['PBS_queue'] = 'batch'
    if 'PBS_nodes' not in conf.keys():
        conf['PBS_nodes'] = '1'
    if 'PBS_ppn' not in conf.keys():
        conf['PBS_ppn'] = '32'
    if 'PBS_walltime' not in conf.keys():
        conf['PBS_walltime'] = '24:00:00'
    if 'PBS_qos' not in conf.keys():
        conf['PBS_qos'] = 'burst'
    if 'PBS_jobname' not in conf.keys():
        print 'ERROR: PBS_jobname is required by PBS!'
        return 1
    if 'PBS_molname' not in conf.keys():
        print 'ERROR: PBS_molname is required by PBS!'
        return 1
    if 'PBS_molname_of_previous_step' not in conf.keys():
        print 'ERROR: PBS_molname_of_previous_step is required by PBS!'
        return 1


    condo_template = '''#PBS -q %s
#PBS -l nodes=%s:ppn=%s,walltime=%s
#PBS -l qos=%s
#PBS -N %s
#PBS -W group_list=%s
#PBS -A %s
#PBS -j oe


mol="%s"

dir=$PBS_O_WORKDIR

module load /software/user_tools/current/modules/cades-bsd/gaussian/g09 env/cades-bsd
module unload PE-gnu
module unload PE-intel
module unload PE-pgi
module load PE-intel
module unload intel/16.0.1
module load intel/17.0.0

export NWCHEM_TOP='/home/p6n/tools/src/nwchem-6.6-i17/'
export nwchem=$NWCHEM_TOP/bin/LINUX64/nwchem


#userscratch=$localscratch/$USER
userscratch=/lustre/or-hydra/cades-bsd/p6n/tmp


export GAUSS_SCRDIR=$userscratch/$PBS_JOBID/$$
mkdir -p $GAUSS_SCRDIR

cd $GAUSS_SCRDIR

cp $dir/* .
mkdir tmp

for m in $mol
do
    mpirun -np $PBS_NP $nwchem $m.nw > $dir/$m.out 2>&1
done

rm -fr tmp
cp -u ./* $dir

rm -fr $GAUSS_SCRDIR

''' % (conf['PBS_queue'], conf['PBS_nodes'], conf['PBS_ppn'],
       conf['PBS_walltime'], conf['PBS_qos'], conf['PBS_jobname'],
       conf['PBS_grouplist'], conf['PBS_account'], conf['PBS_molname']
       )

    return condo_template

def genSLURM_Cascades_EMSL_nw(conf):

    # check and set default
    if 'PBS_queue' not in conf.keys():
        conf['PBS_queue'] = 'batch'
    if 'PBS_nodes' not in conf.keys():
        conf['PBS_nodes'] = '1'
    if 'PBS_ppn' not in conf.keys():
        conf['PBS_ppn'] = '32'
    if 'PBS_walltime' not in conf.keys():
        conf['PBS_walltime'] = '24:00:00'
    if 'PBS_qos' not in conf.keys():
        conf['PBS_qos'] = 'burst'
    if 'PBS_jobname' not in conf.keys():
        print 'ERROR: PBS_jobname is required by PBS!'
        return 1
    if 'PBS_molname' not in conf.keys():
        print 'ERROR: PBS_molname is required by PBS!'
        return 1
    if 'PBS_molname_of_previous_step' not in conf.keys():
        print 'ERROR: PBS_molname_of_previous_step is required by PBS!'
        return 1


    template = '''#!/bin/csh -f
#MSUB -l nodes=%s:ppn=%s,walltime=%s
#MSUB -A st49868
#MSUB -e %s.err
#MSUB -o %s.nwout
#MSUB -N %s
#MSUB -m ea
#MSUB -M lian765@emsl.pnl.gov
#MSUB -V


set mol="%s"


setenv SCRDIR "/dtemp/lian765/${SLURM_JOBID}"

set Pwd=${PWD}
#set nwchem="/home/lian765/tools/nwchem07272017/bin/nwchem-07272017-trunk"


# Standard version
setenv THE_TOPPER                "/home/scicons/cascade/apps/nwchem-6.8/"
setenv NWCHEM                    "/dtemp/scicons/bin/nwchem6.8"
setenv NWCHEM_LIBRARY_DIRECTORY     "${THE_TOPPER}/src/basis/libraries/"
setenv NWCHEM_LIBRARY_DIRECTORY_PW  "${THE_TOPPER}/src/nwpw/libraryps/"


source /etc/profile.d/modules.csh
module purge
module load intel/ips_17_u4
module load impi/5.1.2.150

setenv ARMCI_DEFAULT_SHMMAX 32768
setenv NWCHEM_BASIS_LIBRARY "/home/lian765/tools/nwchem07272017/libraries/"
setenv NWCHEM_NWPW_LIBRARY "/home/scicons/cascade/apps/nwchem-6.6/src/nwpw/libraryps/"
#this disables xeon phi offload
setenv NWC_RANKS_PER_DEVICE 0
#this disables threaded in MKL since it is better to keep it to advanced users
setenv OMP_NUM_THREADS 1
setenv MKL_NUM_THREADS 1
setenv NWC_RANKS_PER_DEVICE 0
setenv ARMCI_OPENIB_DEVICE mlx4_0
setenv OFFLOAD_INIT on_offload

setenv MPIRETURN 999

mkdir -p $SCRDIR
cd $SCRDIR
#cp $nwchem nwchem
cp $Pwd/* .
mkdir tmp

foreach m ($mol)
    echo $m
    #mpirun -n $SLURM_NPROCS ./nwchem $m.nw > $m.out
    srun --mpi=pmi2 -n $SLURM_NPROCS -K1 $NWCHEM  $Pwd/$m.nw > $Pwd/$m.out
end
setenv MPIRETURN $?

rm -fr tmp nwchem
cp -u ./* $Pwd

############################################################################
# End of the job script
############################################################################

exit $MPIRETURN
''' % (conf['PBS_nodes'], conf['PBS_ppn'],
       conf['PBS_walltime'], conf['PBS_jobname'], conf['PBS_jobname'], conf['PBS_jobname'],
       conf['PBS_molname'])



    return template

def genSLURM_Cascades_EMSL_nw_old(conf):

    # check and set default
    if 'PBS_queue' not in conf.keys():
        conf['PBS_queue'] = 'batch'
    if 'PBS_nodes' not in conf.keys():
        conf['PBS_nodes'] = '1'
    if 'PBS_ppn' not in conf.keys():
        conf['PBS_ppn'] = '32'
    if 'PBS_walltime' not in conf.keys():
        conf['PBS_walltime'] = '24:00:00'
    if 'PBS_qos' not in conf.keys():
        conf['PBS_qos'] = 'burst'
    if 'PBS_jobname' not in conf.keys():
        print 'ERROR: PBS_jobname is required by PBS!'
        return 1
    if 'PBS_molname' not in conf.keys():
        print 'ERROR: PBS_molname is required by PBS!'
        return 1
    if 'PBS_molname_of_previous_step' not in conf.keys():
        print 'ERROR: PBS_molname_of_previous_step is required by PBS!'
        return 1


    template = '''#!/bin/csh -f
#MSUB -l nodes=%s:ppn=%s,walltime=%s
#MSUB -A st49868
#MSUB -e %s.err
#MSUB -o %s.nwout
#MSUB -N %s
#MSUB -m ea
#MSUB -M lian765@emsl.pnl.gov
#MSUB -V


set mol="%s"


setenv SCRDIR "/dtemp/lian765/${SLURM_JOBID}"

set Pwd=${PWD}
#set nwchem="/home/lian765/tools/nwchem07272017/bin/nwchem-07272017-trunk"
source /msc/apps/compilers/intel/15.0.090/composer_xe_2015.0.090/bin/compilervars.csh intel64
source /msc/apps/compilers/intel/impi/5.0.1.035/intel64/bin/mpivars.csh intel64


source /etc/profile.d/modules.csh
module purge
module load nwchem/6.6
setenv ARMCI_DEFAULT_SHMMAX 32768
setenv NWCHEM_BASIS_LIBRARY "/home/lian765/tools/nwchem07272017/libraries/"
setenv NWCHEM_NWPW_LIBRARY "/home/scicons/cascade/apps/nwchem-6.6/src/nwpw/libraryps/"
setenv ARMCI_OPENIB_DEVICE mlx4_0
setenv MPIRETURN 999

mkdir -p $SCRDIR
cd $SCRDIR
#cp $nwchem nwchem
cp $Pwd/* .
mkdir tmp

foreach m ($mol)
    echo $m
    #mpirun -n $SLURM_NPROCS ./nwchem $m.nw > $m.out
    srun --mpi=pmi2 -n $SLURM_NPROCS -K1 /dtemp/scicons/bin/nwchem6.6  $Pwd/$m.nw > $Pwd/$m.out
end
setenv MPIRETURN $?

rm -fr tmp nwchem
cp -u ./* $Pwd

############################################################################
# End of the job script
############################################################################

exit $MPIRETURN
''' % (conf['PBS_nodes'], conf['PBS_ppn'],
       conf['PBS_walltime'], conf['PBS_jobname'], conf['PBS_jobname'], conf['PBS_jobname'],
       conf['PBS_molname'])



    return template









def help():
    help_doc = '''There are four steps to use this module.

    #1. import the module
    from PLg09 import pbsPrepare

    #2. set the 'conf' dictionary, follow the example bellow:
    #   NOTE: each calculation should have a conf dict!
    conf = {
        # ... ...

        # required settings
        'PBS_jobname': '',
        'PBS_molname': '',
        'PBS_molname_of_previous_step': '',


        'PBS_queue': 'batch',
        'PBS_nodes': '1',
        'PBS_ppn': '32',
        'PBS_walltime': '12:00:00',
        'PBS_qos': 'burst'

        # ... ...
    }

    #3. generate the PBS submit script for g09 on Condo
    pbs_str = pbsPrepare.genPBS_Condo_g09(conf)

    #4. write the PBS into a file
    output_pbs = '%s/%s.pbs' % (output_path, inp_name)
    with open(output_pbs, 'w') as foutpbs:
        foutpbs.write(pbs_str + '\n')
    foutpbs.close()


    About 'conf':
        'conf' is a python dictionary that collects all the settings required to generate a g09 input file.

'''
    print help_doc


if __name__ == '__main__':
    help()