#!/usr/bin/env python
#
# @purpose
#   to start a new calculation environment.
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 20 2016
#

import subprocess, json, os, time, math, inspect, shutil, logging
from PLg09 import constants, g09prepare, pbsPrepare, g09checkResults
import pandas as pd
from chempy import Reaction, Substance

from CalAndPlot import CalAndPlot
from PrepareAndJobControl import PrepareAndJobControl


def Cal_for_picked_mols(conf, picked_mols=(), supercycle=5):
    logging.info('Starts function "Cal_for_picked_mols"...')

    PrepareJC = PrepareAndJobControl()
    conf = g09prepare.validate_conf(conf)
    logging.basicConfig(level=logging.INFO)

    if len(picked_mols) == 0:
        logging.error('ERROR: No molecules are specified.')
        return

    ## Init for step 1
    #PrepareJC.InitNewCal(conf, deleteDir=True, picked_mols_list=picked_mols)

    ## Check and recal for step 1
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=picked_mols)
    #PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)

    ## cruise for step 1
    try:
        PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True,
                                  plan_B=step1_plan_B, picked_mols_list=picked_mols, sleep_time_min=5)
    except:
        logging.error('\n\nERROR: Something wrong in Cruise for step%s_%s' % (str(conf['Step']), str(conf['Substep'])))
        exit
        #pass

    # <!-- some check function here, to make sure Step 1 is finished -->
    try:
        failed_mols = PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, picked_mols_list=picked_mols, returnMols=True)
    except:
        failed_mols = []
        logging.error('ERROR: Failed to get "failed_mols" for step%s_%s' % (str(conf['Step']), str(conf['Substep'])))

    ## check cycle for step 1
    counter = 1
    while counter <= supercycle:
        if len(failed_mols) > 0:
            logging.info('\nThis is the %d cycle for step%s:' % (counter, str(conf['Step'])))
            logging.info('There are %d failed molecules found:\n%s\n' % (len(failed_mols), ' '.join(failed_mols)))

            # pertubate the xyz and then recalculate
            for mol in failed_mols:
                PrepareJC.perturb_xyz(conf, mol, offset_factor=0.1)

            ## cruise for step 1
            try:
                PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True,
                                              plan_B=step1_plan_B, picked_mols_list=failed_mols, sleep_time_min=5)
            except Exception as e:
                logging.error('\n\nERROR: Something wrong in Cruise for step%s_%s' % (str(conf['Step']), str(conf['Substep'])))
                logging.error('\n\nERROR Message:\n%s\n\n' % e)
                pass

            # <!-- some check function here, to make sure Step 1 is finished -->
            try:
                failed_mols = PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, picked_mols_list=failed_mols, returnMols=True)
            except:
                failed_mols = []
                logging.error('ERROR: Failed to get "failed_mols" for step%s_%s' % (str(conf['Step']), str(conf['Substep'])))

        else:
            logging.warning('No failed mols were found, continue to next step.')
            break

        counter += 1



    #### For step 2 ####
    conf.update(step2_qm_conf)

    ## Init for step 2
    #PrepareJC.InitNewCal(conf, picked_mols_list=picked_mols)
    #PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)

    ## Check for step 2
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=picked_mols)


    ## cruise for step 2
    PrepareJC.Cruise_for_one_step(conf, max_cycles=1, initNew=True, deleteDir=False, gen_inp_for_planB=False,
                                  plan_B=step1_plan_B, picked_mols_list=picked_mols, sleep_time_min=5)


    # <!-- some check function here, to make sure Step 2 is finished -->

    return

def Check_imaginary_freq(conf, picked_mols=()):
    # <!-- some check function here, to make sure there is no imaginary frequency found. If not, collect the molecules
    # and throw them into Cal_for_picked_mol -->
    logging.info('Starts function "Check_imaginary_freq(conf, picked_mols=(%s))"' % ', '.join(picked_mols))

    PrepareJC = PrepareAndJobControl()

    folder_name, folder_path, conf_path, xyz_folder, submit_all_file = PrepareJC.gen_fundamental_vars(conf)
    # find the command
    g09pyGauss = '%s/g09pyGauss.py' % os.path.dirname(inspect.getfile(constants))

    if picked_mols in ['reaction']:
        # read reactions
        reactions, reactants, products = PrepareJC.read_reactions(conf)
        all_strucutres = [i[:-4] for i in os.listdir(xyz_folder) if i.endswith('.xyz')]
        picked_mols = all_strucutres

    mols_with_imaginary_freq = []
    for cpd in picked_mols:
        # plot
        folder = '%s/%s' % (folder_path, cpd)

        plot_cmds = '%s freq %s' % (g09pyGauss, folder)
        freq_str = subprocess.check_output(plot_cmds, shell=True)

        if len([i for i in freq_str.split('\n') if i.startswith('WARN')]) > 0:
            mols_with_imaginary_freq.append(cpd)

            logging.warning('Found mol with imaginary freq %s ' % cpd)
            logging.info(freq_str)

    return mols_with_imaginary_freq





def Cal_for_whole_reaction_dataset(conf):
    logging.info('Starts function "Cal_for_whole_reaction_dataset"...')

    PrepareJC = PrepareAndJobControl()
    Analyze = CalAndPlot()


    ## Init for step 1
    #PrepareJC.InitNewCal(conf, picked_mols_list=['reaction'])

    ## Check and recal for step 1
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=['reaction'])
    #PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)


    ## cruise for step 1
    PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True,
                                  plan_B=step1_plan_B, picked_mols_list=['reaction'], sleep_time_min=5)

    # <!-- some check function here, to make sure Step 1 is finished -->

    #### For step 2 ####
    conf.update(step2_qm_conf)

    ## Init for step 2
    #PrepareJC.InitNewCal(conf, picked_mols_list=['reaction'])
    #PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)

    ## Check for step 2
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=['reaction'])


    ## cruise for step 2
    PrepareJC.Cruise_for_one_step(conf, max_cycles=1, initNew=True, gen_inp_for_planB=False, plan_B=step1_plan_B,
                                  picked_mols_list=['all'], sleep_time_min=5)


    # <!-- some check function here, to make sure Step 2 is finished -->

    # <!-- some check function here, to make sure there is no imaginary frequency found. If not, collect the molecules
    # and throw them into Cal_for_picked_mol -->

    #### For plotting ####
    # Calculate the logK
    Analyze.Cal_logK_for_a_set_of_reactions(conf)

    # plot the figures
    Analyze.Plot_results_for_a_set_of_reactions(conf, overwrite=True)

    # get xyz
    Analyze.Get_xyz_for_a_set_of_reactions(conf)


    return


def Gas_phase_correction(conf, Correction=True):
    PrepareJC = PrepareAndJobControl()
    if Correction:
        #conf.update(step1_gas)

        picked_mols = ['OH_l_g_r_-', 'Hg_l_g_r_+2', 'Cl_l_g_r_-']
        PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True,
                            plan_B=step1_plan_B, picked_mols_list=picked_mols)

def Cal_and_plot_results_for_a_set_of_reactions(conf, DataSets, PlotDataSets=['logK_ALL.txt']):
    for txt in DataSets:
        conf.update({'Reaction_dataset': txt})
        Analyze = CalAndPlot()
        #### For calculating the logK ####
        Analyze.Cal_logK_for_a_set_of_reactions(conf)

        if txt in PlotDataSets:
            #### For plotting the results ####
            Analyze.Plot_results_for_a_set_of_reactions(conf, overwrite=True)
            #### For getting xyz ####
            Analyze.Get_xyz_for_a_set_of_reactions(conf)

def Cal_and_plot_results_for_a_set_of_comformations(conf):
    Analyze = CalAndPlot()
    #### For calculating the logK ####
    #Analyze.Cal_logK_for_a_set_of_reactions(conf)

    #### For plotting the results ####
    Analyze.Plot_results_for_a_set_of_conformations(conf, overwrite=True)
    #### For getting xyz ####
    Analyze.Get_xyz_for_a_set_of_conformations(conf)
    #### For Bolzmann weighting ####
    Analyze.Bolzmann_weighting(conf, eps=0.03, minSamples=1)


def main(conf, picked_mols, supercycle=5):
    #supercycle = 5
    Cal_for_picked_mols(conf, picked_mols=picked_mols, supercycle=supercycle)

    mols_with_imaginary_freq = Check_imaginary_freq(conf, picked_mols=picked_mols)

    try:
        logging.INFO('Molecules with imaginary frequences: %s' % ','.join(mols_with_imaginary_freq))
    except:
        pass

    counter = 1
    while len(mols_with_imaginary_freq) > 0 and counter <= supercycle:
        PrepareJC = PrepareAndJobControl()

        # pertubate the xyz and then recalculate
        for mol in mols_with_imaginary_freq:
            PrepareJC.perturb_xyz(conf, mol, offset_factor=0.1)

        Cal_for_picked_mols(conf, picked_mols=mols_with_imaginary_freq, supercycle=supercycle)

        mols_with_imaginary_freq = Check_imaginary_freq(conf, picked_mols=mols_with_imaginary_freq)

        counter += 1


def Cal_conformers_for_Bolzmann(mols):
    '''

    :param mols: could be unfinished conformers or new molecules for Bolzmann weighting or mix!
                 the prefix '_l_Xe_r_' is required for unfinished conformers.
    :return:
    '''
    # find all molecules
    molecules = []
    for m in mols:
        if m.startswith('_l_Xe_r_'):
            # for restart molecules
            molecules.append('_'.join(m.replace('_l_Xe_r_', '').split('_')[1:]))
        else:
            molecules.append(m)
    molecules = list(set(molecules))

    # build a dict wich contains molecules and all the conformers
    molecules_conformers = {}
    for m in molecules:
        if m in mols:
            # for initiate molecule
            molecules_conformers[m] = ['all']
        else:
            molecules_conformers[m] = [a for a in mols if '_'.join(a.replace('_l_Xe_r_', '').split('_')[1:]) == m]

    # run the calculation
    for m in sorted(molecules_conformers.keys()):
        conf.update({'XYZ_foldername': 'structures/%s/xyz' % m})
        # cal all molecules in the "XYZ_foldername"
        picked_mols = molecules_conformers[m]

        logging.info('\n%s: %s\n' % (m, json.dumps(picked_mols)))

        main(conf, list(set(picked_mols)), supercycle=10)

        '''
        if len(picked_mols) == 1:
            if picked_mols[0] in ['all']:
                # plot
                Cal_and_plot_results_for_a_set_of_comformations(conf)
        '''

def seperated_steps_old():
    #### For step 1 ####
    # combine different confs
    conf = {}
    conf.update(step1_qm_conf)
    conf.update(resource_conf)
    conf.update(project_conf)
    conf.update(project_conf_sub1)

    PrepareJC = PrepareAndJobControl()
    Analyze = CalAndPlot()

    # pick molecules to deal with
    # picked_mols = ['HgCl+', 'HgCl2', 'HgCl3-', 'Ti_HgCl3-', 'HgCl4-2', 'P_HgCl4-2', 'H2O', 'Hg_l_H2O_r_+2',
    # 'Ti_Hg_l_H2O_r_3+2', 'Ti_Hg_l_H2O_r_2Cl+']


    # picked_mols = ['Hg_l_H2O_r_2+2', 'Hg_l_NH3_r_+2', 'Hg_l_NH3_r_2+2', 'Hg_l_NH3_r_3+2', 'Ti_Hg_l_NH3_r_3+2',
    #               'Hg_l_NH3_r_4+2', 'P_Hg_l_NH3_r_4+2', 'Hg_l_H2O_r__l_NH3_r_+2',
    #               'Ti_Hg_l_H2O_r__l_NH3_r_2+2', 'Y_Hg_l_H2O_r_3+2']

    # picked_mols = ['CH3NH2', 'CH3NH3+', 'Hg_l_CH3NH2_r_Cl+', 'Hg_l_CH3NH2_r_+2', 'Hg_l_CH3NH2_r_2+2', 'Hg_l_CH3NH2_r_3+2', 'Hg_l_CH3NH2_r_4+2', ]
    # picked_mols = ['P_Hg_l_CH3COO_r_+', 'P_Hg_l_CH3COO_r_2', 'Hg_l_CH3COO_r__l_H2O_r_+']

    # picked_mols = picked_mols + ['Hg_l_NH3_r_4+2', 'Ti_Hg_l_H2O_r_Cl2', 'Ti_Hg_l_NH3_r_3+2', 'Hg_l_NH3_r_3+2']
    # picked_mols = ['Ti_Hg_l_H2O_r_Cl2', 'HgCl4-2', 'HgCl3-',
    #               'Ti_Hg_l_H2O_r_2_l_NH3_r_+2', 'Hg_l_NH3_r_3+2', 'Ti_Hg_l_NH3_r_3+2', 'P_Hg_l_NH3_r_4+2', 'Hg_l_NH3_r_4+2', 'Ti_Hg_l_H2O_r__l_NH3_r_2+2',
    #               'Y_Hg_l_OH_r__l_H2O_r_2+', 'CH3COO-']

    picked_mols = ['Hg_l_CH3COO_r__l_H2O_r_+', 'P_Hg_l_CH3COO_r__l_H2O_r_+']

    ## Init for step 1
    # PrepareJC.InitNewCal(conf, picked_mols_list=['reaction'])
    # PrepareJC.InitNewCal(conf, deleteDir=True, picked_mols_list=picked_mols)

    ## Check and recal for step 1
    # PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=['reaction'])
    # PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=picked_mols)
    # PrepareJC.Rsync_local_to_remote(conf)
    # PrepareJC.Submit_jobs(conf)

    ## Complement for step 1
    # PrepareJC.InitNewCal(conf, picked_mols_list=['42_a.-1'])
    # PrepareJC.Rsync_local_to_remote(conf)
    # PrepareJC.Submit_jobs(conf)

    ## cruise for step 1
    # PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=['reaction'], sleep_time_min=5)
    # PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=picked_mols, sleep_time_min=5)


    #### For step 2 ####
    conf.update(step2_qm_conf)

    ## Init for step 2
    # PrepareJC.InitNewCal(conf, picked_mols_list=['reaction'])
    # PrepareJC.InitNewCal(conf, picked_mols_list=picked_mols)
    # PrepareJC.Rsync_local_to_remote(conf)
    # PrepareJC.Submit_jobs(conf)

    ## Check for step 2
    # PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=['reaction'])
    PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=picked_mols)

    ## Complement for step 2
    # PrepareJC.InitNewCal(conf, picked_mols_list=['42_a.-1'])
    # PrepareJC.Rsync_local_to_remote(conf)
    # PrepareJC.Submit_jobs(conf)

    ## cruise for step 2
    # PrepareJC.Cruise_for_one_step(conf, max_cycles=1, initNew=True, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=['all'], sleep_time_min=5)
    # PrepareJC.Cruise_for_one_step(conf, max_cycles=1, initNew=True, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=picked_mols, sleep_time_min=5)




    #### For gas phase correction of OH- ####

    Correction = 0
    if Correction:
        conf.update(step1_gas)

        picked_mols = ['OH_l_g_r_-', 'Hg_l_g_r_+2', 'Cl_l_g_r_-']

        # PrepareJC.InitNewCal(conf, deleteDir=True, picked_mols_list=picked_mols)

        ## Check and recal for step 1
        # PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=picked_mols)
        # PrepareJC.Rsync_local_to_remote(conf)
        # PrepareJC.Submit_jobs(conf)

        ## Complement for step 1
        # PrepareJC.InitNewCal(conf, picked_mols_list=['42_a.-1'])
        # PrepareJC.Rsync_local_to_remote(conf)
        # PrepareJC.Submit_jobs(conf)

        ## cruise for step 1
        # PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, gen_inp_for_planB=True, plan_B=step1_plan_B,
        #                    picked_mols_list=['Hg_l_OH_r_3-'], sleep_time_min=5)
        PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True,
                                      plan_B=step1_plan_B, picked_mols_list=picked_mols, sleep_time_min=5)

    #### For plotting ####
    # Calculate the logK
    Analyze.Cal_logK_for_a_set_of_reactions(conf)

    # plot the figures
    Analyze.Plot_results_for_a_set_of_reactions(conf, overwrite=True)

    # get xyz
    Analyze.Get_xyz_for_a_set_of_reactions(conf)


if __name__ == '__main__':
    # Each calculation should have an independent conf dict.
    # One can generate or modify the conf dict using scripts.
    # In this part, only a template is required.
    logging.basicConfig(level=logging.INFO)

    #### ----|||| important part ||||---- ####
    ##### controls which project to work on

    #from Project_conf_M062X import *
    # from Project_conf_M062X_Bondi import *
    from Project_conf_M062X_SAS_Alpha0485 import *

    #### ----|||| important part ||||---- ####

    reaction_sets = ['pKa_cys.txt']
    #reaction_sets = ['Gsolv_bench.txt']


    for r in reaction_sets:
        conf = {}
        conf.update(step1_qm_conf)
        conf.update(resource_conf)
        conf.update(project_conf)
        conf.update(project_conf_sub1)

        #### main cycle ####
        conf.update({'Reaction_dataset': r})


        ######## for Bolzmann weighting, Restart!! #########
        '''
        molecules = [
    "_l_Xe_r_9_Hg_l_SCH2CH2COO_r_2-2",
    "_l_Xe_r_5_Hg_l_SCH2CH2COO_r_2-2",
    "_l_Xe_r_4_Hg_l_SCH2CH2COO_r_2-2"
]


        Cal_conformers_for_Bolzmann(molecules)

        '''
        ######## for normal calculation #########
        # cal all molecules in "Reaction_dataset"

        molecules = ['reaction']
        molecules = [
            "La5021"
            ]


        main(conf, list(set(molecules)), supercycle=10)


    #### For gas phase correction of OH- ####
    #conf.update(step1_gas)
    #Gas_phase_correction(conf)


    #### Cal and Plot ####
    txt_list = reaction_sets

    conf = {}
    conf.update(step1_qm_conf)
    conf.update(resource_conf)
    conf.update(project_conf)
    conf.update(project_conf_sub1)
    conf.update(step2_qm_conf)


    ######## for normal calculation #########
    Cal_and_plot_results_for_a_set_of_reactions(conf, txt_list, PlotDataSets=txt_list)


    '''
    ######### For Bolzmann weighting #########
    # molecules = ['Hg_l_CH3CHSCOO_r_2-2', 'Hg_l_CH3CHSCOO_r_', 'CH3CHSCOO-2', 'Hg_l_CH3CHSCOOHV1_r_+', 'CH3CHSCOOH-']

    #molecules = ['Hg_l_CH3NH2_r_2+2', 'Hg_l_CH3NH2_r_3+2', 'Hg_l_CH3NH2_r_4+2']

    #molecules = ['Hg_l_SCH2CH2COOH_r_+', 'Hg_l_SCH2CH2COO_r_', 'Hg_l_SCH2CH2COO_r_2-2', 'Hg_l_SCHCH2COO2_r_-',
    #             'Hg_l_SCHCH2COOH2_r_+', 'SCHCH2COO2-3', 'SCHCH2COOH2-', 'SHCHCH2COO2-2', 'SHCHCH2COOH2']

    #molecules = ['Hg_l_SCHCH2COO2_r_2-4']
    molecules = ['CH3CHSCOO-2', 'Hg_l_CH3CHSCOO_r_', 'Hg_l_CH3CHSCOO_r_2-2']

    for CSmol in molecules:
        conf.update({'XYZ_foldername': 'structures/%s/xyz' % CSmol})

        Cal_and_plot_results_for_a_set_of_comformations(conf)

    '''
