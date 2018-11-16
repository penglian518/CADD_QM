#!/usr/bin/env python
#
# @purpose
#   configuration file for the whole project
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 20 2016
#
import os

#### For Step 1
step1_qm_conf = {
    # A group of calculation contains many steps
    'Group_name': 'B3LYP_631+Gd_SMD',

    # the step number of this calculation, starts from 1
    'Step': '1',

    # whether start from previous step, if 'Step' is 1, will NOT able to start from previous
    'Start_from_previous': True,
    'Geom': 'Checkpoint',
    'Guess': 'Read',

    # As 0 means False in Python, all numbers here should be input as string!
    'Charge': '0',

    # if not specified, will calculate it according to odd or even number of electrons
    # 'Spin': '1',

    # if not specified, use B3LYP and 6-31G(d)
    'functional': 'B3LYP',
    'basis_set': '6-31+G(d)',

    # put all options for Opt and Freq in braces,
    'Opt': '(MaxCyc=250)',
    'Freq': True,

    # 'Pseudo' is a switch to use pseudo potential, set it to 'Read' or True to turn on
    'Pseudo': True,
    # if turned on, should put the elements in a list.
    'Pseudo_elements': ['Hg'],
    # if turned on, default potential is SDD
    'Pseudo_potential': 'SDD',
    # if turned on, the option of Gen should put here and separated by blank. e.g. 'Gen 5D'
    'Pseudo_basis': '',

    # similart to Opt, Freq, put all options in braces for 'SCF'
    'SCF': '(MaxCyc=250, Tight)',
    'Integral': 'UltraFine',
    'NoSymm': True,

    # similart to 'Pseudo', 'SCRF' is a switch
    'SCRF': True,
    'SCRF_model': 'SMD',
    'SCRF_solvent': 'water',
    'SCRF_others': '',

    # others options could be put here
    'Others': 'pop=NBO'
}

resource_conf = {
    # number of processors to be used
    'NProc': '32',
    # memory
    'Mem': '10GB',

    # chk file name will be the same as the input file

    # path to access the coordinates in .xyz file,
    # will be read by pybel
    # the filename will be used as the name of .com .chk .log
    'path_to_input_xyz': ''
}

pbs_conf = {

    # required settings
    'PBS_jobname': '',
    'PBS_molname': '',
    'PBS_molname_of_previous_step': '',

    'PBS_queue': 'batch',
    'PBS_nodes': '1',
    'PBS_ppn': '32',
    'PBS_walltime': '12:00:00',
    'PBS_qos': 'burst'
}

project_conf = {
    'XYZ_folderpath': '..',
    'XYZ_foldername': 'structures',

    'Local_pwd': os.path.abspath('.'),
    'Local_calculation_folder_name': 'calculations',
    'Local_output_folder_name': 'output',

    'Remote_cluster_name': 'condo',
    'Remote_calculation_folder_name': '/home/p6n/workplace/Takat'
}

project_conf_sub1 = {
    'Reaction_foldername': 'reactions',
    'Reaction_dataset': 'logK_ALL.txt'
}





#### For step 1 plan B
# if step 1 is failed, use this conf to re-cal
step1_plan_B = {
    # the step number of this calculation, starts from 1
    'Substep': '2',

    # whether start from previous step, if 'Step' is 1, will NOT able to start from previous
    'Start_from_previous': True,
    'Geom': 'Checkpoint',
    'Guess': 'Harris',

    # put all options for Opt and Freq in braces,
    'Opt': '(MaxCyc=250, Cartesian, GDIIS)',
    'Freq': False,

    'NoSymm': True,
    # others options could be put here
    'Others': 'IOP(1/8=18)'
}


#### For step 2
step2_qm_conf = {
    # the step number of this calculation, starts from 1
    'Step': '2',

    # whether start from previous step, if 'Step' is 1, will NOT able to start from previous
    'Start_from_previous': True,
    'Geom': 'Checkpoint',
    'Guess': 'Read',

    # put all options for Opt and Freq in braces,
    'Opt': False,
    'Freq': True,

    'Others': 'pop=NBO'
}



#### For gas phase calculation for step 1
step1_gas = {
    'Step': '1',
    'Substep': '1',
    'Start_from_previous': False,
    'Opt': '(MaxCyc=250)',
    'Freq': True,
    'SCRF': True,
    'Others': 'pop=NBO'

}



def conf_builder_help():
    doc_str = '''
Build conf with the 'update' method of a python dictionary.

Example:
    # For step 1, build the dict from start
    conf = {}
    conf.update(step1_qm_conf)
    conf.update(resource_conf)
    conf.update(project_conf)
    conf.update(project_conf_sub)

    # For step 2, just update the keys that need to be changed
    conf.update(step2_qm_conf)
'''

    print doc_str


if __name__ == '__main__':
    conf_builder_help()
