#!/usr/bin/env python
#
# @purpose
#   to prepare input files for g09, according to the conf dict
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 20 2016
#
import os, copy
import pybel

def validate_conf(conf):
    '''
    Check the inputs and form a new conf dictionary.

    1. conf should be a dict
    2. to remove indicators/switches which show False.
        NOTE: in python False means 0, so for the int 0, use a string 0 instead.
    3. set default configurations for functional, basis set, pseudo potential, scrf
    4. set default for number of processors and memory to allocate

    :param conf: configuration dict
    :return: updated conf
    '''

    # check if the conf is a dict
    if type(conf) != dict:
        print 'WARN: a dict contain all configure for a g09 job is required, but got %s' % type(conf)
        return {}

    # clean the conf, remove options that is False or [] or ''
    newconf = {}
    for k in conf.keys():
        if conf[k] not in [False, '', []]:
            newconf[k] = conf[k]

    # for qm conf
    # check if Group name is exist
    if 'Group_name' not in newconf.keys():
        print 'WARN: no Group_name was found, using default "cal_1"'
        newconf['Group_name'] = 'cal_1'

    # check if Step is exist
    if 'Step' not in newconf.keys():
        print 'WARN: no Step was found, using default "1"'
        newconf['Step'] = '1'

    # check if Substep is exist
    if 'Substep' not in newconf.keys():
        print 'WARN: no Substep was found, using default "1"'
        newconf['Substep'] = '1'

    # check Start_from_previous
    if 'Start_from_previous' in newconf.keys():
        if newconf['Start_from_previous']:
            if int(newconf['Step']) < 2 and int(newconf['Substep']) < 2:
                print 'WARN: no previous step is found, set "Start_from_previous" to False'
                # False will be removed in the next validation run!!
                # So, hear using '0' instead. In the test function please use conf['Start_from_previous'] == True.
                newconf['Start_from_previous'] = '0'
            else:
                # check if Geom is exist
                if 'Geom' not in newconf.keys():
                    print 'WARN: no Geom was found, using default "Checkpoint".'
                    newconf['Geom'] = 'Checkpoint'

                ## sometimes read previous guess is not necessary
                ## check if Guess is exist
                #if 'Guess' not in newconf.keys():
                #    print 'WARN: no Guess was found, using default "Read".'
                #    newconf['Guess'] = 'Read'

    # check if functional is exist
    if 'functional' not in newconf.keys():
        print 'WARN: no functional was found, using the default one "B3LYP"'
        newconf['functional'] = 'B3LYP'

    # check if basis set is exist
    if 'basis_set' not in newconf.keys():
        print 'WARN: no basis set was found, using the default one "6-31G(d)"'
        newconf['basis_set'] = '6-31G(d)'

    # check if pseudo is exist
    if 'Pseudo' in newconf.keys():
        if 'Pseudo_potential' not in newconf.keys():
            newconf['Pseudo_potential'] = 'SDD'
        if 'Pseudo_basis' not in newconf.keys():
            newconf['Pseudo_basis'] = 'Gen'

    # check if scrf is exist
    if 'SCRF' in newconf.keys():
        if 'SCRF_model' not in newconf.keys():
            newconf['SCRF_model'] = 'PCM'
        if 'SCRF_solvent' not in newconf.keys():
            newconf['SCRF_solvent'] = 'water'
        if 'SCRF_Read' in newconf.keys():
            if 'SCRF_ReadConf' not in newconf.keys():
                newconf['SCRF_ReadConf'] = 'Radii=UFF'

    # for resources conf
    # check if nproc is exist
    if 'NProc' not in newconf.keys():
        print 'WARN: no NProc was found, using the default one "2"'
        newconf['NProc'] = '2'

    # check if Mem is exist
    if 'Mem' not in newconf.keys():
        print 'WARN: no Mem was found, using the default one "100MB"'
        newconf['Mem'] = '500MB'

    return newconf

def lower_and_chrop_blank_braces(string):
    return string.lower().replace(' ','').replace('\t', '').replace(',','').replace('(', '').replace(')', '')




def genFolderName(conf):
    '''
    Generate folder name according to the configure file of the calculation.
    Generally, the folder name is formed by several parts.


    type_of_cal  func  basis  [pseudo  pseudo_element]  scrf  solvent
    opt_b3lyp_6-31gd_sdd_Hg_smd_water


    1. validate the conf
    2. deal with the parts one by one
    3. combine them together to form a routine string

    :param conf: configure dict for g09 job
    :return: the folder name string
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        return ''

    foldername_list = []
    confkeys = [i.lower() for i in conf.keys()]

    # type of cal
    if 'opt' in confkeys and 'freq' in confkeys:
        foldername_list.append('opt_freq')
    elif 'opt' in confkeys:
        foldername_list.append('opt')
    elif 'freq' in confkeys:
        foldername_list.append('freq')
    else:
        foldername_list.append('cal')

    # functional
    foldername_list.append(lower_and_chrop_blank_braces(conf['functional']))

    # basis set
    foldername_list.append(lower_and_chrop_blank_braces(conf['basis_set']))

    # pseudo, if possible list the elements those using ECP
    if 'pseudo' in confkeys:
        foldername_list.append(lower_and_chrop_blank_braces(conf['Pseudo_potential']))
        try:
            ecp_ele_str = '_'.join(conf['Pseudo_elements'])
            foldername_list.append(ecp_ele_str)
        except:
            pass

    # scrf
    if 'scrf' in confkeys:
        foldername_list.append(lower_and_chrop_blank_braces(conf['SCRF_model']))
        foldername_list.append(lower_and_chrop_blank_braces(conf['SCRF_solvent']))

    # combine these options with '_'
    foldername_str = '_'.join(foldername_list)

    #print foldername_str
    return foldername_str

def genHeader(conf):
    '''
    Generate header of the g09 input file
    1. validate the conf
    2. if starts from previous, set %OldChk
    3. if path_to_input_xyz is provided, then set the name of .chk file (to the same as .xyz file).
        Otherwise, neglect the .chk file

    :param conf: configure dict for g09 job
    :return: the header string
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        return ''

    header_list = []
    # %NProc
    header_list.append('%%NProc=%s' % conf['NProc'])
    # %Mem
    header_list.append('%%Mem=%s' % conf['Mem'])

    #gen %OldChk
    if 'Start_from_previous' in conf.keys() and conf['Start_from_previous'] == True:
        input_filename = os.path.basename(conf['path_to_input_xyz'])[:-4]
        header_list.append('%%OldChk=step%s_%s_%s.chk' % (conf['preStep'], conf['preSubstep'], input_filename))

    # gen %Chk
    if 'path_to_input_xyz' in conf.keys():
        input_filename = os.path.basename(conf['path_to_input_xyz'])[:-4]
        header_list.append('%%Chk=step%s_%s_%s.chk' % (conf['Step'], conf['Substep'], input_filename))

    header_str = '\n'.join(header_list)

    #print header_str
    return header_str

def genRoutine(conf):
    '''
    Generate calculation Routine according to the configure file.
    Currently, divide the routine into the follow parts, and will be processed one by one,

    init func  basis      opt              [pseudo]      scf                  integral          nosymm   scrf
    #P B3LYP/6-31G(d) Opt=(MaxCyc=250) [Pseudo=Read] SCF=(MaxCyc=250, Tight) Integral=UltraFine NoSymm SCRF=(SMD, solvent=water)

    1. validate the conf
    2. deal with the parts one by one
    3. combine them together to form a routine string


    :param conf: configure dict for g09 job
    :return: the calculation routine string
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        return ''

    routine_list = []
    init_str = '#P'

    confkeys = [i.lower() for i in conf.keys()]

    # functional and basis set
    func_basis_str = '%s/%s' % (conf['functional'], conf['basis_set'])

    # type of cal
    cal_type_str = ''
    if 'opt' in confkeys:
        if conf['Opt'] in [True, 1]:
            cal_type_str = 'Opt'
        else:
            cal_type_str = 'Opt=%s' % conf['Opt']
    if 'freq' in confkeys:
        if conf['Freq'] in [True, 1]:
            cal_type_str += ' Freq'
        else:
            cal_type_str += ' Freq=%s' % conf['Freq']

    # pseudo
    if 'pseudo' in confkeys:
        pseudo_str = 'Pseudo=Read'
        func_basis_str = '%s/%s' % (conf['functional'], conf['Pseudo_basis'])

    # SCF
    if 'scf' in confkeys:
        scf_str = 'SCF=%s' % conf['SCF']

    # Integral
    if 'integral' in confkeys:
        intgeral_str = 'Integral=%s' % conf['Integral']

    # NoSymm
    if 'nosymm' in confkeys:
        if conf['NoSymm'] in [True, 1]:
            nosymm_str = 'NoSymm'

    # scrf
    if 'scrf' in confkeys:
        if 'scrf_read' in confkeys:
            if 'scrf_others' in confkeys:
                scrf_str = 'SCRF=(%s, solvent=%s, %s, Read)' % (conf['SCRF_model'], conf['SCRF_solvent'], conf['SCRF_others'])
            else:
                scrf_str = 'SCRF=(%s, solvent=%s, Read)' % (conf['SCRF_model'], conf['SCRF_solvent'])
        else:
            if 'scrf_others' in confkeys:
                scrf_str = 'SCRF=(%s, solvent=%s, %s)' % (conf['SCRF_model'], conf['SCRF_solvent'], conf['SCRF_others'])
            else:
                scrf_str = 'SCRF=(%s, solvent=%s)' % (conf['SCRF_model'], conf['SCRF_solvent'])

    # start from previous calculation
    if 'Start_from_previous' in conf.keys() and conf['Start_from_previous'] == True:
        geom_str = 'Geom=%s' % conf['Geom']
        if 'guess' in confkeys:
            gauss_str = 'Guess=%s' % conf['Guess']

    # others
    if 'others' in confkeys:
        others_str = conf['Others']


    # if variables exists add them to routine_list
    for var in ['init_str', 'func_basis_str', 'cal_type_str', 'pseudo_str', 'scf_str', 'intgeral_str', 'nosymm_str',
              'scrf_str', 'geom_str', 'gauss_str', 'others_str']:
        try:
            s = eval(var)
            routine_list.append(s)
        except NameError:
            pass


    # combine these options with '_'
    routine_str = ' '.join(routine_list)

    print routine_str
    return routine_str

def genChargeSpin(conf):
    '''
    Generate charge and spin strings for g09 input file.
    1. validate the conf
    2. if no charge settings, set to '0'
    3. if spin is setted by user, use it directly, otherwise, calculated it according number of electrons
    4. check the input path. NOTE: input file is necessary!
    5. read the coordinate file with pybel
    6. calculate the total number of electrons (sum of atomic numbers + charge)
    7. calculate the spin, number of single electrons + 1

    NOTE:
        The procedure for calculate spin is simple.
        Specify the conf['Spin'] value if intend to have higher multiplicity.

    :param conf:
    :return:
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        return ''

    # check charge
    if 'Charge' not in conf.keys():
        print 'WARN: No Charge is found, set default value 0'
        conf['Charge'] = '0'


    # check spin
    if 'Spin' in conf.keys():
        print 'WARN: Using specified SPIN, %s' % conf['Spin']
    else:
        print 'WARN: No specified SPIN is found, trying to calculate it.'

        # check input path
        if 'path_to_input_xyz' not in conf.keys():
            print 'WARN: No path is found to access to the xyz file, which is required to calculate the spin.'
            return ''

        # read the coordinates
        mol = pybel.readfile('xyz', conf['path_to_input_xyz']).next()

        # calculate the total number of electrons
        total_atomicnum = sum([i.atomicnum for i in mol.atoms])
        electron_num = total_atomicnum + int(conf['Charge'])

        # set spin, odd electron --> 2, even --> 1
        conf['Spin'] = str(electron_num % 2 + 1)

        print 'WARN: No specified SPIN is found, assign %s automatically' % conf['Spin']


    charge_str = '%s %s' % (conf['Charge'], conf['Spin'])

    return charge_str

def genCoordinates(conf):
    '''
    Generate coordinates strings for g09 input file.
    1. validate the conf
    2. check the input path. NOTE: input file is necessary!
    3. read the coordinate file with pybel
    4. get the coordinate part and return

    :param conf:
    :return:
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        return ''

    # check input path
    if 'path_to_input_xyz' not in conf.keys():
        print 'WARN: No path is found to access to the xyz file.'
        return ''

    # read the coordinates
    mol = pybel.readfile('xyz', conf['path_to_input_xyz']).next()
    #mol.OBMol.SetTotalCharge(conf['Charge'])
    #mol.OBMol.SetTotalSpinMultiplicity(conf['Spin'])

    # get the coordinate part
    raw_str = mol.write('gjf')
    coor_str = '\n'.join(raw_str.split('\n')[5:-2])

    return coor_str

def genECP(conf):
    '''
    Generate ECP strings in the bottom part of g09 input file.
    1. validate the conf
    2. check if Pseudo potential is called
    3. check if pseudo element list is provided
    4. generate and check the coordinates
    5. figure out ECP elements and non-ECP elements
    6. generate the final ECP string

    :param conf:
    :return:
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        return ''

    # check pseudo
    if 'Pseudo' not in conf.keys():
        print 'WARN: No pseudo-potential requirement is found.'
        return ''

    # check pseudo element
    if 'Pseudo_elements' not in conf.keys():
        print 'WARN: No element is to be applied with pseudo-potential.'
        return ''

    # generate coordinates first!
    try:
        coor_str = genCoordinates(conf)
    except:
        print 'WARN: Failed to generate coordinates for xyz file %s' % conf['path_to_input_xyz']
        return ''

    # check coor_str
    if len(coor_str) == 0:
        print 'WARN: No coordinates is found in xyz file %s' % conf['path_to_input_xyz']
        return ''

    # all elements in this input
    all_elements = set([i.split()[0] for i in coor_str.split('\n') if len(i) > 0])

    # ECP elements and non ECP ele
    # need a deepcopy of the value, not a link!
    ecp_ele = copy.deepcopy(conf['Pseudo_elements'])
    non_ecp_ele = [i for i in all_elements if i not in ecp_ele]


    ecp_ele.append(str(0))
    non_ecp_ele.append(str(0))

    ecp_ele_str = ' '.join(ecp_ele)
    non_ecp_ele_str = ' '.join(non_ecp_ele)


    # generate final ECP str
    # no ECP element, ecp_ele will only contain ['0']
    if len(ecp_ele) == 1:
        print 'WARN: No element is to be applied with pseudo-potential..'
        return ''

    # all ECP element, non_ecp_ele will contain only ['0']
    if len(non_ecp_ele) == 1:
        ecp_str = '%s\n%s\n****\n\n%s\n%s' % (ecp_ele_str, conf['Pseudo_potential'],
                                              ecp_ele_str, conf['Pseudo_potential'])
    else:
        # have both ECP and non ECP element
        ecp_str = '%s\n%s\n****\n%s\n%s\n****\n\n%s\n%s' % (non_ecp_ele_str, conf['basis_set'],
                                                            ecp_ele_str, conf['Pseudo_potential'],
                                                            ecp_ele_str, conf['Pseudo_potential'])

    return ecp_str



def genInputFile(conf):
    '''
    Generate the input file of g09.
    1. validate the conf
    2. check the input path. NOTE: input file is necessary!
    3. generate header routine charge spin coordinates sequencially
    4. check if starts from previous or not
    5. check if Pseudo potential is used

    :param conf:
    :return:  g09 input string, suggest_inp_name, conf
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        return ''

    # check input path
    if 'path_to_input_xyz' not in conf.keys():
        print 'WARN: No path is found to access to the xyz file.'
        return ''

    # get the strings for each part
    header_str = genHeader(conf)
    routine_str = genRoutine(conf)
    mol_name = os.path.basename(conf['path_to_input_xyz'])[:-4]
    charge_str = genChargeSpin(conf)
    coor_str = genCoordinates(conf)


    # start from previous or not
    if 'Start_from_previous' in conf.keys() and conf['Start_from_previous'] == True:
        # if pseudo potential is required
        if 'Pseudo' in conf.keys():
            ecp_str = genECP(conf)
            if 'SCRF_Read' in conf.keys():
                scrf_read_str = conf['SCRF_ReadConf']
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            '',
                            ecp_str,
                            '',
                            scrf_read_str,
                            ''
                            '']
            else:
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            '',
                            ecp_str,
                            '',
                            '']

        else:
            if 'SCRF_Read' in conf.keys():
                scrf_read_str = conf['SCRF_ReadConf']
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            '',
                            scrf_read_str,
                            '',
                            '']
            else:
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            '',
                            '']
    else:
        # if pseudo potential is required
        if 'Pseudo' in conf.keys():
            ecp_str = genECP(conf)
            if 'SCRF_Read' in conf.keys():
                scrf_read_str = conf['SCRF_ReadConf']
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            coor_str,
                            '',
                            ecp_str,
                            '',
                            scrf_read_str,
                            '',
                            '']
            else:
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            coor_str,
                            '',
                            ecp_str,
                            '',
                            '']

        else:
            if 'SCRF_Read' in conf.keys():
                scrf_read_str = conf['SCRF_ReadConf']
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            coor_str,
                            '',
                            scrf_read_str,
                            '',
                            '']
            else:
                # form the final input string
                inp_list = [header_str,
                            '',
                            routine_str,
                            '',
                            mol_name,
                            '',
                            charge_str,
                            coor_str,
                            '',
                            '']
    # combine
    inp_str = '\n'.join(inp_list)

    # suggest input file name, without suffix
    suggest_inp_name = 'step%s_%s_%s' % (conf['Step'], conf['Substep'], mol_name)
    return inp_str, suggest_inp_name



def help():
    help_doc = '''There are four steps to use this module.

    #1. import the module
    from PLg09 import g09prepare

    #2. set the 'conf' dictionary, follow the example bellow:
    #   NOTE: each calculation should have a conf dict!
    conf = {
        # ... ...

        'functional': 'B3LYP',
        'basis_set': '6-31G(d)'

        # ... ...
    }

    #3. generate the g09 input file according to 'conf'
    g09inputfile = g09prepare.genInputFile(conf)

    #4. write the string into a file
    with open('output_file', 'w') as fout:
        fout.write(g09inputfile + '\n')
    fout.close


    About 'conf':
        'conf' is a python dictionary that collects all the settings required to generate a g09 input file.


    A Example to form conf:

    qm_conf = {
        # A group of calculation contains many steps
        'Group_name': 'cal_1',

        # the step number of this calculation, starts from 1
        'Step': '2',

        # whether start from previous step, if 'Step' is 1, will NOT able to start from previous
        'Start_from_previous': True,
        'Geom': 'Checkpoint',
        'Guess': 'Read',

        # As 0 means False in Python, all numbers here should be input as string!
        'Charge': '0',

        # if not specified, will calculate it according to odd or even number of electrons
        #'Spin': '1',

        # if not specified, use B3LYP and 6-31G(d)
        'functional': 'B3LYP',
        'basis_set': '6-31G(d)',

        # put all options for Opt and Freq in braces,
        'Opt': '(MaxCyc=250)',
        'Freq': False,

        # 'Pseudo' is a switch to use pseudo potential, set it to 'Read' or True to turn on
        'Pseudo': False,
        # if turned on, should put the elements in a list.
        'Pseudo_elements': ['Hg'],
        # if turned on, default potential is SDD
        'Pseudo_potential': 'SDD',
        # if turned on, the option of Gen should put here and separated by blank. e.g. 'Gen 5D'
        'Pseudo_basis': 'Gen',

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
        'Others': ''
    }

    resource_conf = {
        # number of processors to be used
        'NProc': '32',
        # memory
        'Mem': '7000MB',

        # chk file name will be the same as the input file

        # path to access the coordinates in .xyz file,
        # will be read by pybel
        # the filename will be used as the name of .com .chk .log
        'path_to_input_xyz': ''
    }


    project_conf = {
        'XYZ_folderpath': '../',
        'XYZ_foldername': 'structures'
    }


    # combine different confs
    conf = {}
    conf.update(qm_conf)
    conf.update(resource_conf)
    conf.update(project_conf)

'''
    print help_doc




if __name__ == '__main__':
    help()