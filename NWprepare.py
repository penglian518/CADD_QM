#!/usr/bin/env python
#
# @purpose
#   to prepare input files for NWChem, according to the conf dict
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 20 2016
#
import os, copy, logging, re
import pybel
from g09prepare import validate_conf, lower_and_chrop_blank_braces, genFolderName

logging.basicConfig(level=logging.INFO)

def suggest_jobname_title(conf):
    # check input path
    if 'path_to_input_xyz' not in conf.keys():
        logging.error('ERROR: No path is found to access to the xyz file.')
        return ''

    if 'NEB' in conf.keys() and 'NEB_path_name' in conf.keys():
        mol_name = conf['NEB_path_name']
    elif 'String' in conf.keys() and 'String_path_name' in conf.keys():
        mol_name = conf['String_path_name']
    elif 'PythonStringEn' in conf.keys() and 'Path_name' in conf.keys():
        mol_name = 'StringEN_%s' % conf['Path_name']
    elif 'CalStringEn' in conf.keys() and 'Path_name' in conf.keys():
        mol_name = '%s_%s' % (conf['Path_name'], str(conf['ith_Bead']))
    else:
        mol_name = os.path.basename(conf['path_to_input_xyz'])[:-4]

    title = 'step%s_%s_%s' % (conf['Step'], conf['Substep'], mol_name)

    #gen %OldChk
    if 'Start_from_previous' in conf.keys() and conf['Start_from_previous'] == True:
        title = 'step%s_%s_%s' % (conf['preStep'], conf['preSubstep'], mol_name)


    #job_name = '%s_%s' % (conf['Group_name'], mol_name)
    job_name = '%s' % mol_name


    return job_name, title

def convBasis_G09toNW(g09basis):
    # to lower and remove blank space
    basis = g09basis.lower().replace(' ', '')

    # gaussian type
    if basis.startswith('6-31') or basis.startswith('3-21'):
        nwbasis = basis.replace('(d)', '*').replace('(p)', '*').replace('(d,p)','**')
    elif basis == 'sdd':
        nwbasis = 'stuttgart_rsc_1997_ecp'
    elif basis == 'sdd_g09':
        nwbasis = 'sdd_from_gaussian'
    else:
        nwbasis = basis
        logging.info('INFO: unkonwn type of basis sets, "%s". ' % nwbasis)

    return nwbasis



def genHeader(conf):
    '''
    Generate header of the NWChem input file
    1. validate the conf
    2. if starts from previous, set restart
    3. if path_to_input_xyz is provided, then set the name of .chk file (to the same as .xyz file).
        Otherwise, neglect the .chk file

    :param conf: configure dict for g09 job
    :return: the header string
    '''

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    header_list = ['echo']

    ## for start new calculation
    # jobname and title
    Jname, title = suggest_jobname_title(conf)
    if Jname == '':
        logging.info('Failed to gen title for the job. Check "path_to_input_xyz"')
        return ''

    # for restart
    if 'Start_from_previous' in conf.keys() and conf['Start_from_previous'] == True:
        # start
        header_list.append('restart %s' % Jname)
    else:
        # start
        header_list.append('start %s' % Jname)

    # title
    header_list.append('title %s' % title)
    # scratch
    header_list.append('scratch_dir ./tmp')
    # memory
    header_list.append('memory %s %s' % (''.join(re.findall('[0-9]+', conf['Mem'])), ''.join(re.findall('[a-zA-Z]+', conf['Mem']))))

    header_str = '\n'.join(header_list)

    #print header_str
    return header_str

def getChargeSpin(conf):
    '''
    Generate charge and spin strings for nwchem input file.
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
        logging.info('WARN: No Charge is found, set default value 0')
        conf['Charge'] = '0'


    # check spin
    if 'Spin' in conf.keys():
        logging.info('WARN: Using specified SPIN, %s' % conf['Spin'])
    else:
        logging.info('WARN: No specified SPIN is found, trying to calculate it.')

        # check input path
        if 'path_to_input_xyz' not in conf.keys():
            logging.error('ERROR: No path is found to access to the xyz file, which is required to calculate the spin.')
            return ''

        # read the coordinates
        mol = pybel.readfile('xyz', conf['path_to_input_xyz']).next()

        # calculate the total number of electrons
        total_atomicnum = sum([i.atomicnum for i in mol.atoms])
        electron_num = total_atomicnum + int(conf['Charge'])

        # set spin, odd electron --> 2, even --> 1
        conf['Spin'] = str(electron_num % 2 + 1)

        logging.info('WARN: No specified SPIN is found, assign %s automatically' % conf['Spin'])

    charge_str = conf['Charge']
    spin_str = conf['Spin']

    return charge_str, spin_str

def genCharge(conf):
    charge, spin = getChargeSpin(conf)
    charge_str = 'charge %s' % str(charge)
    return charge_str


def genCoordinates(conf, path_to_xyz='path_to_input_xyz'):
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
    if path_to_xyz not in conf.keys():
        logging.info('WARN: No path is found to access to the xyz file.')
        return ''

    # read the coordinates
    mol = pybel.readfile('xyz', conf[path_to_xyz]).next()
    #mol.OBMol.SetTotalCharge(conf['Charge'])
    #mol.OBMol.SetTotalSpinMultiplicity(conf['Spin'])

    # get the coordinate part
    raw_str = mol.write('gjf')
    coor_str = '\n'.join(raw_str.split('\n')[5:-2])

    return coor_str

def genGeometry(conf):
    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    # generate NEB geometry or NOT
    if 'neb' in confkeys:
        if 'neb_path_to_startgoem' in confkeys:
            # start geom
            geo_list = ['geometry nocenter noautoz noautosym']
            coor_str = genCoordinates(conf, path_to_xyz='NEB_path_to_startgoem')
            geo_list.append(coor_str)
            geo_list.append('end\n')

        if 'neb_path_to_midgoem' in confkeys:
            # mid geom
            geo_list += ['geometry midgeom nocenter noautoz noautosym']
            coor_str = genCoordinates(conf, path_to_xyz='NEB_path_to_midgoem')
            geo_list.append(coor_str)
            geo_list.append('end\n')

        if 'neb_path_to_endgoem' in confkeys:
            # end geom
            geo_list += ['geometry endgeom nocenter noautoz noautosym']
            coor_str = genCoordinates(conf, path_to_xyz='NEB_path_to_endgoem')
            geo_list.append(coor_str)
            geo_list.append('end')
    elif 'string' in confkeys:
        if 'string_path_to_startgoem' in confkeys:
            # start geom
            geo_list = ['geometry nocenter noautoz noautosym']
            coor_str = genCoordinates(conf, path_to_xyz='String_path_to_startgoem')
            geo_list.append(coor_str)
            geo_list.append('end\n')

        if 'string_path_to_midgoem' in confkeys:
            # mid geom
            geo_list += ['geometry midgeom nocenter noautoz noautosym']
            coor_str = genCoordinates(conf, path_to_xyz='String_path_to_midgoem')
            geo_list.append(coor_str)
            geo_list.append('end\n')

        if 'string_path_to_endgoem' in confkeys:
            # end geom
            geo_list += ['geometry endgeom nocenter noautoz noautosym']
            coor_str = genCoordinates(conf, path_to_xyz='String_path_to_endgoem')
            geo_list.append(coor_str)
            geo_list.append('end')
    elif 'geometry_load' in confkeys:
        # NoSymm
        if 'nosymm' in confkeys:
            if conf['NoSymm'] in [True, 1]:
                geo_list = ['geometry units angstroms noautoz noautosym']
        else:
            geo_list = ['geometry units angstroms noautoz']

        # load command
        load_str = '    load frame %d %s' % (conf['ith_Bead'], conf['Geometry_load'])
        geo_list.append(load_str.replace('PATH_NAME', conf['Path_name']))

        # end
        geo_list.append('end')
    else:
        # NoSymm
        if 'nosymm' in confkeys:
            if conf['NoSymm'] in [True, 1]:
                geo_list = ['geometry units angstroms noautoz noautosym']
        else:
            geo_list = ['geometry units angstroms noautoz']

        # coordinates
        coor_str = genCoordinates(conf)
        geo_list.append(coor_str)

        # end
        geo_list.append('end')

    geo_str = '\n'.join(geo_list)

    return geo_str




def genBasis(conf):

    basispart_list = ['basis']
    ecppart_list = ['ecp']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    # check input path
    if 'path_to_input_xyz' not in conf.keys():
        logging.info('WARN: No path is found to access to the xyz file.')
        return ''

    # read the coordinates
    mol = pybel.readfile('xyz', conf['path_to_input_xyz']).next()
    ele_Table = pybel.ob.OBElementTable()

    # get element labels
    atomicNumbers = sorted(list(set([i.atomicnum for i in mol.atoms])))
    all_elements = [ele_Table.GetSymbol(i) for i in atomicNumbers]

    # pseudo
    if 'pseudo' in confkeys:
        for ele in all_elements:
            if ele in conf['Pseudo_elements']:
                if 'NOAuto_Basis_set' not in conf.keys():
                    basispart_list.append('    %s library %s' % (ele, convBasis_G09toNW(conf['Pseudo_potential'])))
                # ecp part
                ecppart_list.append('    %s library %s' % (ele, convBasis_G09toNW(conf['Pseudo_potential'])))
            else:
                if 'NOAuto_Basis_set' not in conf.keys():
                    basispart_list.append('    %s library %s' % (ele, convBasis_G09toNW(conf['basis_set'])))
    else:
        for ele in all_elements:
            if 'NOAuto_Basis_set' not in conf.keys():
                basispart_list.append('    %s library %s' % (ele, convBasis_G09toNW(conf['basis_set'])))



    # additional keywords for basis set
    if 'basis_set_others_nw' in confkeys:
        basispart_list.append('    %s' % conf['basis_set_others_nw'])

    # ending
    basispart_list.append('end')

    # ecp part ending
    if len(ecppart_list) > 1:
        ecppart_list.append('end')

        basispart_list.append('')
        basispart_list += ecppart_list

    # convert to str
    basispart_str = '\n'.join(basispart_list)

    return basispart_str

def genDISP(conf):
    disppart_list = ['disp']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    # scf g09
    try:
        opt_str = conf['Others']
    except:
        opt_str = ''

    opt_list = opt_str.lower().split(' ')

    flag_empiricaldispersion = False
    for s in opt_list:
        if s.startswith('empiricaldispersion'):
            flag_empiricaldispersion = True

            disp = s.split('=')[1]
            if disp in ['gd3bj', 'd3bj']:
                disppart_list.append(' vdw 4')
            elif disp in ['d3']:
                disppart_list.append(' vdw 3')

    # if no empirical dispersion option, turn off this function in NWChem
    if not flag_empiricaldispersion:
        disppart_list.append('    off')

    if not flag_empiricaldispersion:
        disppart_list = []

    disppart_str = ' '.join(disppart_list)

    return disppart_str


def genDFT(conf):

    dftpart_list = ['dft']
    dftpart_list.append('    decomp')

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    # scf g09
    try:
        scf_str = conf['SCF']
    except:
        scf_str = ''

    scf_list = scf_str.lower().replace('(', '').replace(')', '').replace(' ', '').split(',')

    for s in scf_list:
        if s.startswith('maxcyc'):
            dftpart_list.append('    iterations %s' % s.split('=')[1])
        elif s.startswith('tight'):
            dftpart_list.append('    convergence density 1e-8')
            dftpart_list.append('    convergence energy 1e-6')
            dftpart_list.append('    tolerances tight')

    # grid
    if 'integral' in confkeys:
        grid_str = conf['Integral'].lower()

        if grid_str in ['ultrafine']:
            dftpart_list.append('    grid xfine')
        elif grid_str in ['', 'fine']:
            dftpart_list.append('    grid fine')

    # XC
    xc_str = conf['functional'].lower()

    if xc_str in ['b3lyp', 'pbe0']:
        dftpart_list.append('    xc %s' % xc_str)
    elif xc_str in ['m06l', 'm06-l']:
        dftpart_list.append('    xc m06-l')
    elif xc_str in ['m062x', 'm06-2x']:
        dftpart_list.append('    xc m06-2x')
    elif xc_str in ['m052x', 'm05-2x']:
        dftpart_list.append('    xc m05-2x')
    elif xc_str in ['b3pw91']:
        dftpart_list.append('    xc hfexch 0.20 slater 0.80 becke88 nonlocal 0.72 perdew91 0.81 pw91lda 1.00')
    elif xc_str in ['blyp']:
        dftpart_list.append('    xc becke88 lyp')
    else:
        logging.error('ERROR: Cannot find a equivalent XC in NWChem. Please define yourself')
        return ''

    # disp
    disp_str = genDISP(conf)
    if len(disp_str) > 0:
        dftpart_list.append('    %s' % disp_str)

    # mult
    charge, mult = getChargeSpin(conf)
    dftpart_list.append('    mult %s' % mult)

    # additional keyword for DFT
    if 'dft_others_nw' in confkeys:
        dftpart_list.append('    %s' % conf['DFT_others_nw'])

    # ending
    dftpart_list.append('end')
    dftpart_str = '\n'.join(dftpart_list)

    return dftpart_str

def genRelativistic(conf):
    relativisticpart_list = ['relativistic']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    if 'zora' in confkeys:
        relativisticpart_list.append('    zora on')
    if 'zora_others' in confkeys:
        relativisticpart_list.append('    %s' % conf['Zora_others'])
    if 'modelpotential' in confkeys:
        relativisticpart_list.append('    modelpotential %s' % str(conf['Modelpotential']))

    # ending
    relativisticpart_list.append('end')
    relativisticpart_str = '\n'.join(relativisticpart_list)

    return relativisticpart_str

def genSO(conf):
    sopart_list = ['so']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    if 'so_str' in confkeys:
        sopart_list.append(conf['SO_str'])
    # ending
    sopart_list.append('end')
    sopart_str = '\n'.join(sopart_list)

    return sopart_str


def genDriver(conf):
    driverpart_list = ['driver']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    ## title
    #Jname, title = suggest_jobname_title(conf)
    #driverpart_list.append('    xyz %s' % title)

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    # scf g09
    try:
        opt_str = conf['Opt']
    except:
        opt_str = ''

    opt_list = opt_str.lower().replace('(', '').replace(')', '').replace(' ', '').split(',')

    for s in opt_list:
        if s.startswith('maxcyc'):
            driverpart_list.append('    maxiter %s' % s.split('=')[1])
        elif s.startswith('tight'):
            driverpart_list.append('    GMAX 0.000015\n    GRMS 0.00001\n    XMAX 0.00006\n    XRMS 0.00004')
        elif s.startswith('nwstep'):
            if 'saddle' in confkeys:
                driverpart_list.append('    sadstp %s' % s.split('=')[1])
            else:
                driverpart_list.append('    trust %s' % s.split('=')[1])

    # initial hessian
    if 'guess' in confkeys:
        if conf['Guess'] in ['Harris', 'harris']:
            driverpart_list.append('    inhess 1')
        elif conf['Guess'] in ['Read', 'read']:
            driverpart_list.append('    inhess 0')

    # initial hessian
    if 'nwinhess' in confkeys:
        driverpart_list.append('    inhess %s' % str(conf['NWInhess']))
        driverpart_list.append('    maxiter 250')

    # additional keyword for Opt
    if 'opt_others_nw' in confkeys:
        driverpart_list.append('    %s' % conf['Opt_others_nw'])

    # ending
    driverpart_list.append('end')
    driverpart_str = '\n'.join(driverpart_list)

    return driverpart_str

def genFreq(conf):
    freqpart_list = ['freq']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]


    if 'nwreusehess' in confkeys:
        freqpart_list.append('    reuse')

    if 'nwfreqanimate' in confkeys:
        freqpart_list.append('    animate')

    # ending
    freqpart_list.append('end')
    freqpart_str = '\n'.join(freqpart_list)

    return freqpart_str

def genCosmo(conf):
    cosmopart_list = ['cosmo']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    if 'scrf' in confkeys:
        if 'scrf_solvent' in confkeys:
            cosmopart_list.append('    solvent %s' % conf['SCRF_solvent'])
        if 'scrf_model' in confkeys:
            if conf['SCRF_model'] in ['SMD', 'smd']:
                cosmopart_list.append('    do_cosmo_smd .true.')
        if 'scrf_others_nw' in confkeys:
            cosmopart_list.append('    %s' % conf['SCRF_others_nw'])
        if 'scrf_parameters_nw' in confkeys:
            cosmopart_list.append('    parameters %s' % conf['SCRF_parameters_nw'])

        # ending
        cosmopart_list.append('end')
        cosmopart_str = '\n'.join(cosmopart_list)
    else:
        cosmopart_str = ''

    return cosmopart_str

def genNEB(conf):
    nebpart_list = ['neb']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    if 'neb_nbeads' in confkeys:
        nebpart_list.append('    nbeads %s' % str(conf['NEB_nbeads']))
    if 'neb_kbeads' in confkeys:
        nebpart_list.append('    kbeads %s' % str(conf['NEB_kbeads']))
    if 'neb_maxiter' in confkeys:
        nebpart_list.append('    maxiter %s' % str(conf['NEB_maxiter']))
    if 'neb_stepsize' in confkeys:
        nebpart_list.append('    stepsize %s' % str(conf['NEB_stepsize']))
    if 'neb_convergence' in confkeys and conf['NEB_convergence'] not in ['default', 'Default']:
        nebpart_list.append('    %s' % str(conf['NEB_convergence']))
    if 'neb_xyzpath' in confkeys:
        nebpart_list.append('    xyz_path ./%s' % str(os.path.basename(conf['NEB_xyzpath'])))
    if 'neb_print_shift' in confkeys:
        nebpart_list.append('    print_shift %s' % str(conf['NEB_print_shift']))


    # ending
    nebpart_list.append('end')
    nebpart_str = '\n'.join(nebpart_list)

    return nebpart_str

def genString(conf):
    '''string method'''
    stringpart_list = ['string']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    if 'string_nbeads' in confkeys:
        stringpart_list.append('    nbeads %s' % str(conf['String_nbeads']))
    if 'string_maxiter' in confkeys:
        stringpart_list.append('    maxiter %s' % str(conf['String_maxiter']))
    if 'string_stepsize' in confkeys:
        stringpart_list.append('    stepsize %s' % str(conf['String_stepsize']))
    if 'string_interpol' in confkeys:
        stringpart_list.append('    interpol %s' % str(conf['String_interpol']))
    if 'string_freeze1' in confkeys:
        stringpart_list.append('    freeze1 %s' % str(conf['String_freeze1']))
    if 'string_freezen' in confkeys:
        stringpart_list.append('    freezen %s' % str(conf['String_freezen']))
    if 'string_tol' in confkeys:
        stringpart_list.append('    %s' % str(conf['String_tol']))
    if 'string_xyzpath' in confkeys:
        stringpart_list.append('    xyzpath %s' % str(conf['String_xyzpath']))
    if 'string_print_shift' in confkeys:
        stringpart_list.append('    print_shift %s' % str(conf['String_print_shift']))


    # ending
    stringpart_list.append('end')
    stringpart_str = '\n'.join(stringpart_list)

    return stringpart_str

def genPythonStringEn(conf):
    '''Python script'''
    pythonpart_list = ['python']

    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    if 'python_script' in confkeys and 'path_name' in confkeys and 'nbeads' in confkeys and 'theory' in confkeys:
        script = conf['Python_script']
        script = script.replace('PATH_NAME', conf['Path_name']).replace('NBEADS', str(1+conf['NBeads'])).replace('THEORY', conf['Theory'])
        pythonpart_list.append(script)
    else:
        print conf

    # ending
    pythonpart_list.append('end')
    pythonpart_str = '\n'.join(pythonpart_list)

    return pythonpart_str



def genTask(conf):
    # validate the conf
    conf = validate_conf(conf)
    if len(conf) < 1:
        logging.error('Failed to validate conf!')
        return ''

    # all conf keys
    confkeys = [i.lower() for i in conf.keys()]

    taskpart_list = []

    if 'nwcalhess' in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft hessian')
        else:
            taskpart_list.append('task scf hessian')

    if 'energy' in confkeys and 'theory' in confkeys:
        taskpart_list.append('task %s energy' % conf['Theory'])

    if 'opt' in confkeys and 'saddle' not in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft optimize')
        else:
            taskpart_list.append('task scf optimize')

    if 'saddle' in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft saddle')
        else:
            taskpart_list.append('task scf saddle')


        #taskpart_list.append('task scf hessian')

    if 'freq' in confkeys:
        if 'scrf' in confkeys:
            #taskpart_list.append('task dft hessian')
            taskpart_list.append('task dft frequencies')
        else:
            #taskpart_list.append('task scf hessian')
            taskpart_list.append('task scf frequencies')

    if 'nwfindsaddle' in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft saddle')
        else:
            taskpart_list.append('task scf saddle')

    if 'neb' in confkeys and 'neb_ignore' in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft neb ignore')
        else:
            taskpart_list.append('task scf neb ignore')

    if 'neb' in confkeys and 'neb_ignore' not in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft neb')
        else:
            taskpart_list.append('task scf neb')

    if 'string' in confkeys and 'string_ignore' in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft string ignore')
        else:
            taskpart_list.append('task scf string ignore')

    if 'string' in confkeys and 'string_ignore' not in confkeys:
        if 'scrf' in confkeys:
            taskpart_list.append('task dft string')
        else:
            taskpart_list.append('task scf string')

    if 'pythonstringen' in confkeys:
        taskpart_list.append('task python')


    # ending
    taskpart_str = '\n'.join(taskpart_list)

    return taskpart_str



def genInputFile(conf):
    '''
    Generate input file for NWChem.
    1. validate the conf
    2. check the input path. NOTE: input file is necessary!
    3. generate header routine charge spin coordinates sequencially
    4. check if starts from previous or not
    5. check if Pseudo potential is used

    :param conf:
    :return:  NWChem input string, suggest_inp_name, conf
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
    charge_str = genCharge(conf)
    geometry_str = genGeometry(conf)
    basis_str = genBasis(conf)
    dft_str = genDFT(conf)
    relativistic_str = genRelativistic(conf)
    so_str = genSO(conf)
    driver_str = genDriver(conf)
    freq_str = genFreq(conf)
    neb_str = genNEB(conf)
    string_str = genString(conf)
    pythonstringen_str = genPythonStringEn(conf)
    cosmo_str = genCosmo(conf)
    task_str = genTask(conf)

    mol_name = os.path.basename(conf['path_to_input_xyz'])[:-4]

    # start from previous or not
    if 'Start_from_previous' in conf.keys() and conf['Start_from_previous'] == True:
        inp_list = [header_str, '']

        if 'Opt' in conf.keys():
            inp_list += [driver_str, '']
        if 'NWFindSaddle' in conf.keys():
            inp_list += [driver_str, '']
        if 'Freq' in conf.keys():
            inp_list += [freq_str, '']
        if 'Relativistic' in conf.keys():
            inp_list += [relativistic_str, '']
        if 'SO' in conf.keys():
            inp_list += [so_str, '']
        if 'NEB' in conf.keys():
            inp_list += [neb_str, '']
            mol_name = conf['NEB_path_name']
        if 'String' in conf.keys():
            inp_list += [string_str, '']
            mol_name = conf['String_path_name']
        if 'Python' in conf.keys():
            inp_list += [pythonstringen_str, '']
            mol_name = conf['Path_name']
        if 'CalStringEn' in conf.keys():
            mol_name = conf['Path_name']

        inp_list.append(task_str)

    else:
        inp_list = [header_str, '']
        inp_list += [charge_str, '']
        if 'NoGeometry' not in conf.keys():
            inp_list += [geometry_str, '']
        inp_list += [basis_str, '']
        if 'SO' in conf.keys():
            inp_list += [so_str, '']
        inp_list += [dft_str, '']

        if 'SCRF' in conf.keys():
            inp_list += [cosmo_str, '']
        if 'Opt' in conf.keys():
            inp_list += [driver_str, '']
        if 'Freq' in conf.keys():
            inp_list += [freq_str, '']
        if 'Relativistic' in conf.keys():
            inp_list += [relativistic_str, '']
        if 'NEB' in conf.keys():
            inp_list += [neb_str, '']
            mol_name = conf['NEB_path_name']
        if 'String' in conf.keys():
            inp_list += [string_str, '']
            mol_name = conf['String_path_name']
        if 'PythonStringEn' in conf.keys():
            inp_list += [pythonstringen_str, '']
            mol_name = conf['Path_name']
        if 'CalStringEn' in conf.keys():
            mol_name = conf['Path_name']

        inp_list.append(task_str)

    # combine
    inp_str = '\n'.join(inp_list)

    # suggest input file name, without suffix
    if 'CalStringEn' in conf.keys():
        suggest_inp_name = '%s_%s' % (mol_name, conf['ith_Bead'])
    else:
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