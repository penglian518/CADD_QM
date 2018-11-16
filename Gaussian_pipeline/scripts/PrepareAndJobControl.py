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
import numpy as np
#from chempy import Reaction, Substance

class PrepareAndJobControl:
    def __init__(self):
        logging.basicConfig(level=logging.INFO)

    def Rsync_log(self, conf, mol):
        Local_folder = '%s/%s/%s/%s/' % (conf['Local_pwd'], conf['XYZ_folderpath'], conf['Local_calculation_folder_name'], conf['Group_name'])
        Remote_folder = '%s/%s' % (conf['Remote_calculation_folder_name'], conf['Group_name'])
        Cmd = 'rsync -av --include="*/" --include="*.log" --include="*.com" --exclude="*" %s:%s/%s %s' % \
              (conf['Remote_cluster_name'], Remote_folder, mol, Local_folder)
    
        logging.info('\nRsync log file from %s:' % conf['Remote_cluster_name'])
        logging.info(Cmd + '\n')
    
        os.system(Cmd)
        # It's weird that the subprocess does not do the filtering.
        #Output = subprocess.check_output(Cmd.split())
        #print Output
        #return Output

    def Rsync_local_to_remote(self, conf):
        Local_folder = '%s/%s/%s/%s' % (conf['Local_pwd'], conf['XYZ_folderpath'], conf['Local_calculation_folder_name'], conf['Group_name'])
        Remote_folder = '%s/' % (conf['Remote_calculation_folder_name'])
        Cmd = 'rsync -av %s %s:%s' % \
              (Local_folder, conf['Remote_cluster_name'], Remote_folder)

        logging.info('\nRsync the local files to %s:' % conf['Remote_cluster_name'])
        logging.info(Cmd + '\n')

        os.system(Cmd)

    def Submit_jobs(self, conf):
        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        Remote_folder = '%s/%s/' % (conf['Remote_calculation_folder_name'], conf['Group_name'])

        Cmd = 'ssh %s "%s/%s"' % \
              (conf['Remote_cluster_name'], Remote_folder, os.path.basename(submit_all_file))

        # run the command and capture the output
        Output = subprocess.check_output(Cmd.split())

        # format the output
        Output_list = []
        i = 0
        while i < len(Output.strip().split('\n')):
            if i % 2:
                Output_list.append((Output.strip().split('\n')[i - 1], Output.strip().split('\n')[i]))
            i += 1

        # show the results
        logging.info('\nRun the submit script %s:' % os.path.basename(submit_all_file))
        logging.info(Cmd + '\n')
        logging.info(Output_list)

        #os.system(Cmd)

        # in the format of [(step1_1_10_a.-1, 34030.or-condo-pbs01),(),(),...]
        return Output_list

    def Check_current_jobs(self, conf):
        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        Cmd = 'ssh %s /opt/torque/bin/qstat -u p6n' % conf['Remote_cluster_name']

        # run the command and capture the output
        Output = subprocess.check_output(Cmd.split())

        Output_list = []
        if len(Output) > 0:
            # put the queue info into a dataframe
            df = pd.DataFrame([i.split() for i in Output.split('\n')[5:-1]])
            # filter the completed jobs, and get the current job IDs
            Output_list = list(df[df.iloc[:, 9] != 'C'].iloc[:, 0].values)

        return Output_list

    def Delete_remote_mols(self, conf, mols):
        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        Remote_folder = '%s/%s/' % (conf['Remote_calculation_folder_name'], conf['Group_name'])

        mols_tobe_deleted = ['%s/%s' % (Remote_folder, m) for m in mols]

        Cmd = 'ssh %s rm -fr %s' % (conf['Remote_cluster_name'], ' '.join(mols_tobe_deleted))

        # run the command and capture the output
        Output = subprocess.check_output(Cmd.split())

        # show the results
        if len(Output) > 0:
            logging.info(Output)

        logging.info('\nDeleted remote molecules %s:' % ' '.join(mols_tobe_deleted))
        return



    def gen_fundamental_vars(self, conf):
        '''
        This function services as a supplement of conf dict

        :param conf:
        :return:
        '''
        # folder_name = 'calculations'
        folder_name = conf['Local_calculation_folder_name']
        # folder to save all calculations in this group, relative path
        folder_path = '../%s/%s' % (folder_name, conf['Group_name'])
        # folder to save every 'conf's
        conf_path = '%s/conf' % folder_path

        # make the calculation folder
        try:
            os.makedirs(conf_path)
        except:
            pass

        # input molecules, in .xyz format
        xyz_folder = '%s/%s/%s' % (conf['Local_pwd'], conf['XYZ_folderpath'], conf['XYZ_foldername'])

        # script to submit all jobs
        submit_all_file = '%s/step%s_submitall.sh' % (folder_path, conf['Step'])

        return folder_name, folder_path, conf_path, xyz_folder, submit_all_file

    def gen_submit_all_script(self, mol_names, step_num):
        '''

        :param mol_names: should be a string of mol names separated by ' ' . e.g. 'mol1 mol2 mol3 mol4 ...'
        :param step_num:
        :return:
        '''

        # prepare the submit scripts
        submit_all_str = '''#!/bin/bash

        molecules="%s"

        for mol in $molecules
        do
            cd $mol
            qsub step%s_${mol}.pbs
            cd ..
        done
        ''' % (mol_names, str(step_num))

        return submit_all_str

    def gen_submit_all_script_v2(self, mol_step_list):
        '''

        :param mol_step_list: e.g. '[(mol1, 1_1), (mol2, 1_2), (mol3, 2_1) ...]'

        :return:
        '''

        # prepare the submit scripts
        submit_all_list = ['#!/bin/bash\n\nDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"\ncd $DIR\n']
        for mol in mol_step_list:
            submit_all_list.append('cd %s; echo step%s_%s; qsub step%s_%s.pbs; cd ..' % (mol[0], mol[1], mol[0], mol[1], mol[0]))
        submit_all_list.append('exit\n')

        submit_all_str = '\n'.join(submit_all_list)

        return submit_all_str

    def gen_com_pbs(self, conf, xyz, xyz_folder, folder_path, conf_path, CalReaction=False, CheckPseudo=True, deleteDir=False):
        '''

        :param conf:
        :param xyz: molecule name
        :param xyz_folder: abs path to the folder that contains xyz files
        :param folder_path: relative path to the folder that contains this group of calculations
        :param conf_path:  path to the folder that contains the back conf files
        :return:
        '''
        # setting output for each molecule
        output_path = '%s/%s' % (folder_path, xyz)
        if deleteDir:
            try:
                os.removedirs(output_path)
                logging.info('Will delete dir: %s' % output_path)
            except:
                logging.info('Failed to delete dir : %s' % output_path)
                logging.info('Will delete the files in dir : %s' % output_path)
                try:
                    for i in os.listdir(output_path):
                        try:
                            os.remove('%s/%s' % (output_path, i))
                            logging.info('Deleting: %s/%s' % (output_path, i))
                        except:
                            pass
                except:
                    pass
        try:
            os.mkdir(output_path)
        except:
            pass

        # determine conf['preStep'] and conf['preSubstep'], according to .log files in the folder!
        if int(conf['Step']) == 1:
            conf['preStep'] = conf['Step']
            conf['preSubstep'] = str(len([i for i in os.listdir(output_path) if i[-4:] == '.log' and i.startswith('step%s_' % conf['Step'])]))
        elif conf['Step'] > 1:
            pre_substep_num = len([i for i in os.listdir(output_path) if i[-4:] == '.log' and i.startswith('step%s_' % conf['Step'])])
            if pre_substep_num == 0:
                conf['preStep'] = str(int(conf['Step']) - 1)
                conf['preSubstep'] = str(len(
                    [i for i in os.listdir(output_path) if i[-4:] == '.log' and i.startswith('step%s_' % conf['preStep'])]))
            else:
                conf['preStep'] = conf['Step']
                conf['preSubstep'] = str(pre_substep_num)


        ## G09 input
        # input coordinates file
        xyz_file = '%s/%s.xyz' % (xyz_folder, xyz)
        logging.info('\nDealing with %s' % xyz_file)
        conf['path_to_input_xyz'] = xyz_file

        # molecule charge
        # read charge from the note line in xyz file
        mol_charge = int(open(xyz_file).readlines()[1].split()[0])
        conf['Charge'] = mol_charge

        # use Pseudo potential or not
        if CheckPseudo:
            with open(conf['path_to_input_xyz']) as fin:
                fcon = fin.readlines()
            fin.close()

            # get the elements
            elements_in_xyz = list(set([i.strip().split()[0] for i in fcon[2:] if len(i) > 1]))
            #print elements_in_xyz

            # common element within conf
            element_common = [i for i in elements_in_xyz if i in conf['Pseudo_elements']]
            #print element_common


            if len(element_common) == 0:
                conf['Pseudo'] = False
            else:
                conf['Pseudo'] = True


        # generate the input string
        g09inp_str, inp_name = g09prepare.genInputFile(conf)



        # write the inp file
        output_file = '%s/%s.com' % (output_path, inp_name)
        with open(output_file, 'w') as fout:
            fout.write(g09inp_str + '\n')
        fout.close()

        ## PBS submit
        conf['PBS_jobname'] = inp_name
        conf['PBS_molname'] = inp_name
        conf['PBS_molname_of_previous_step'] = 'step%s_%s_%s' % (conf['preStep'], conf['preSubstep'], xyz)

        pbs_str = pbsPrepare.genPBS_Condo_g09(conf)

        # write the pbs file
        output_pbs = '%s/%s.pbs' % (output_path, inp_name)
        with open(output_pbs, 'w') as foutpbs:
            foutpbs.write(pbs_str + '\n')
        foutpbs.close()

        # save the conf for debug purpose
        conf_output_file = '%s/%s.conf' % (conf_path, inp_name)
        json.dump(conf, open(conf_output_file, 'w'), sort_keys=True, indent=4)

    def perturb_xyz(self, conf, xyz, offset_factor=0.1):
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        xyz_file = '%s/%s.xyz' % (xyz_folder, xyz)

        # exist or not
        if not os.path.exists(xyz_folder):
            logging.error('ERROR: Cannot find xyz file at %s' % xyz_file)
            return 1

        # backup the original copy
        N_copys_of_xyz = len([i for i in os.listdir(xyz_folder) if i.startswith('%s.xyz' % xyz)])
        shutil.copyfile(xyz_file, '%s.%d' % (xyz_file, N_copys_of_xyz))

        # read coord
        mol = [i.strip().split() for i in open(xyz_file).readlines()]

        mol_atoms = np.array(mol[2:])
        mol_coors = mol_atoms[:, 1:4].astype(float)

        # gen offset values for mol coor.
        # The offset_factor should better not larger than 0.5
        offset_matrix = np.random.random(mol_coors.shape) * offset_factor
        # combine
        mol_coors_off = mol_coors + offset_matrix
        # writeback coordinates
        mol_atoms[:, 1:4] = mol_coors_off
        # form new molecule
        mol_new = mol[:2] + mol_atoms.tolist()

        # write the new mol to xyz file
        with open(xyz_file, 'w') as fout:
            for line in [' '.join(i) for i in mol_new]:
                fout.write('%s\n' % line)
            fout.close()




    def read_reactions(self, conf):

        reactions_file = '%s/%s/%s' % (conf['XYZ_folderpath'], conf['Reaction_foldername'], conf['Reaction_dataset'])

        # if the reaction dataset file exist
        if not os.path.exists(reactions_file):
            logging.info('\n%s not exist!\n' % reactions_file)
            return []

        # open the dataset
        with open(reactions_file) as fin:
            fcon = fin.readlines()
        fin.close()

        # read reactions
        fcon = [l for l in fcon if not l.startswith('#')]
        reactions = [Reaction.from_string(i.strip()) for i in fcon if len(i) > 1]

        # reactants
        all_reac = [i.reac.keys() for i in reactions]
        reactants = []
        for r in all_reac:
            reactants += r

        # products
        all_prod = [i.prod.keys() for i in reactions]
        products = []
        for p in all_prod:
            products += p

        # drop duplicates
        reactants = list(set(reactants))
        products = list(set(products))


        # change '(' to '__l__', and ')' to '__r__'
        reactants_bashsafe = [r.replace('(', '_l_').replace(')', '_r_') for r in reactants]
        products_bashsafe = [p.replace('(', '_l_').replace(')', '_r_') for p in products]

        return reactions, reactants_bashsafe, products_bashsafe



    def InitNewCal(self, conf, picked_mols_list=('all'), deleteDir=False):
        '''
        Prepare input and PBS script for mols in the picked_mols_list.

        It used to start new calculations. Either for step 1 or the following steps

        :param picked_mols_list: only mol name is required, no .xyz is required. e.g. ['mol1', 'mol3']
        '''
        logging.info('Starts function "InitNewCal"...')

        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        if conf['Pseudo']:
            CheckPseudo = True
        else:
            CheckPseudo = False

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        # pick molecules to be calculated in this list, without suffix '.xyz'
        if len(picked_mols_list) == 1 and picked_mols_list[0] == 'all':
            CalReaction = False
            xyz_file_names = [i[:-4] for i in os.listdir(xyz_folder) if i.endswith('.xyz')]

        elif len(picked_mols_list) == 1 and picked_mols_list[0] in ['reaction', 'reactions', 'Reaction', 'Reactions']:
            CalReaction = True

            # read reactions
            reactions, reactants, products = self.read_reactions(conf)
            all_strucutres = [i[:-4] for i in os.listdir(xyz_folder) if i.endswith('.xyz')]

            # check if structure file exist
            xyz_file_names = []
            for cpd in list(set(reactants + products)):
                if cpd in all_strucutres:
                    xyz_file_names.append(cpd)
                else:
                    logging.error('Cannot find the strucutre file for %s' % cpd)
                    return 1
        else:
            #CalReaction = False
            CalReaction = True
            xyz_file_names = picked_mols_list

        # generate g09 input files
        mol_step_list = []
        for xyz in xyz_file_names:
            if xyz.endswith('_l_g_r_'):
                conf.update({'SCRF': False})
            else:
                conf.update({'SCRF': True})
            self.gen_com_pbs(conf, xyz, xyz_folder, folder_path, conf_path, CalReaction, CheckPseudo, deleteDir)
            mol_step_list.append((xyz, '%s_%s' % (conf['Step'], conf['Substep'])))

        # prepare the submit scripts
        #submit_all_str = self.gen_submit_all_script(' '.join(xyz_file_names), '%s_%s' % (conf['Step'], conf['Substep']))
        submit_all_str = self.gen_submit_all_script_v2(mol_step_list)

        # write the submitall script
        with open(submit_all_file, 'w') as fout_submit:
            fout_submit.write(submit_all_str + '\n')
        fout_submit.close()

        # make it executable
        os.system('chmod +x %s' % submit_all_file)


        # delete remote mols
        if deleteDir:
            self.Delete_remote_mols(conf, xyz_file_names)

        return

    def Check_and_recal(self, conf, gen_inp_for_planB=True, plan_B={}, picked_mols_list=('all'), returnMols=False):
        '''
        The purpose is to keep the calculations in that step finish, if not create a new substep with planB conf.

        Better turn off gen_inp_for_planB, if there is no planB configurations.


        The logical is list bellow:
        Starts from substep 1, check log files one by one, till find the 'first' one which is normally finished.
        If none of them finished normally, use plan_B update conf then generate the input, submit scripts for next substep of calculation

        :param conf:
        :return:
        '''
        logging.info('Starts function "Check_and_recal" ...')

        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        if conf['Pseudo']:
            CheckPseudo = True
        else:
            CheckPseudo = False


        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        # pick molecules to be calculated in this list, without suffix '.xyz'
        if len(picked_mols_list) == 1 and picked_mols_list[0] == 'all':
            CalReaction = False
            xyz_file_names = [i[:-4] for i in os.listdir(xyz_folder) if i.endswith('.xyz')]
        elif len(picked_mols_list) == 1 and picked_mols_list[0] in ['reaction', 'reactions', 'Reaction', 'Reactions']:
            CalReaction = True

            # read reactions
            reactions, reactants, products = self.read_reactions(conf)
            all_strucutres = [i[:-4] for i in os.listdir(xyz_folder) if i.endswith('.xyz')]

            # check if structure file exist
            xyz_file_names = []
            for cpd in list(set(reactants + products)):
                if cpd in all_strucutres:
                    xyz_file_names.append(cpd)
                else:
                    logging.error('Cannot find the strucutre file for %s' % cpd)
                    return 1
        else:
            #CalReaction = False
            CalReaction = True
            xyz_file_names = picked_mols_list

        # get the current log files
        for mol in xyz_file_names:
            self.Rsync_log(conf, mol)


        # check g09 log file one by one
        failed_mols = []
        mol_step_list = []
        for mol in xyz_file_names:
            ## check G09 output
            # path for each mol
            mol_path = '%s/%s' % (folder_path, mol)

            if not os.path.exists(mol_path):
                logging.warning('\nWARN: %s is not exist!' % mol_path)
                continue

            # all log files for this mol, for conf['Step']
            log_files = [i for i in os.listdir(mol_path) if i[-4:] == '.log' and i.startswith('step%s_' % conf['Step'])]
            if len(log_files) == 0:
                logging.warning('WARN: no log files were found for %s' % mol)
                continue

            #### core logical bellow
            # check log one by one
            i_log = 1
            step_num_to_use = False
            while i_log <= len(log_files):
                log_file = '%s/step%s_%s_%s.log' % (mol_path, conf['Step'], i_log, mol)
                return_from_checkfail = g09checkResults.checkfail(log_file, silence=True)

                # if log finished normally then jump to next mol, else check next log file
                if return_from_checkfail == 'pass':
                    step_num_to_use = False
                    break
                else:
                    i_log += 1
                    step_num_to_use = i_log
                    continue

            # step_num_to_use not False means all log files of the mol not finish normally.
            # step1 failed
            if step_num_to_use:
                # give possible reason for fail
                logging.info('\n' + ' '.join(return_from_checkfail))
                #print mol, step_num_to_use

                # save the failed mol names
                failed_mols.append(mol)

                # generate input for plan B
                if gen_inp_for_planB:
                    # prepare for step1 plan B
                    conf.update(plan_B)
                    conf['Substep'] = step_num_to_use

                    # generate the input and pbs scripts for plan B
                    self.gen_com_pbs(conf, mol, xyz_folder, folder_path, conf_path, CalReaction)
                    mol_step_list.append((mol, '%s_%s' % (conf['Step'], conf['Substep'])))

        if gen_inp_for_planB:
            # if there are some mols assign to take plan B then prepare submit all scripts for them
            if len(failed_mols) > 0:
                # prepare the submit all scripts
                #submit_all_str = self.gen_submit_all_script(' '.join(failed_mols), conf['Step'])
                submit_all_str = self.gen_submit_all_script_v2(mol_step_list)
            else:
                # all finished, nothing to submit
                submit_all_str = '#!/bin/bash\necho "Nothing to submit.\n"'

            # write the submit all script
            with open(submit_all_file, 'w') as fout_submit:
                fout_submit.write(submit_all_str + '\n')
            fout_submit.close()

            # make it executable
            os.system('chmod +x %s' % submit_all_file)

        # prepare for return code
        if len(failed_mols) > 0:
            All_mols_are_converged = False
            logging.warning('The following mols are not converged:\n%s' % ', '.join(failed_mols))
        else:
            All_mols_are_converged = True

        if returnMols:
            return failed_mols
        else:
            return All_mols_are_converged


    def JobControl_upload_submit_and_waite_for_done(self, conf, sleep_time_min=20):
        logging.info('Starts function "JobControl_upload_submit_and_waite_for_done(conf, sleep_time_min=%d)"...' % sleep_time_min)
        # upload to server
        self.Rsync_local_to_remote(conf)

        # submit the jobs
        jobs_list = self.Submit_jobs(conf)
        submited_jobIDs = [i[1] for i in jobs_list]

        # check current jobs
        time.sleep(10)
        current_all_jobIDs = self.Check_current_jobs(conf)
        current_running_jobIDs = [i for i in submited_jobIDs if i in current_all_jobIDs]

        # waite until all the jobs finish
        while len(current_running_jobIDs) > 0:
            logging.info('%d jobs are left at %s.' % (len(current_running_jobIDs), time.ctime()))
            time.sleep(sleep_time_min * 60)

            current_all_jobIDs = self.Check_current_jobs(conf)
            current_running_jobIDs = [i for i in submited_jobIDs if i in current_all_jobIDs]

        logging.info('All calculations are finished at %s' % time.ctime())


    def Cruise_for_one_step(self, conf, max_cycles=10, initNew=True, deleteDir=False, gen_inp_for_planB=True,
                            plan_B={}, picked_mols_list=('all'), sleep_time_min=20):
        logging.info('Starts function "Cruise_for_one_step"...')

        #### init the job
        if initNew:
            self.InitNewCal(conf, picked_mols_list, deleteDir=deleteDir)

        # upload files, submit jobs, and waite for done
        self.JobControl_upload_submit_and_waite_for_done(conf, sleep_time_min=sleep_time_min)

        #### check for step 1
        All_converged = self.Check_and_recal(conf, gen_inp_for_planB, plan_B, picked_mols_list)

        # if not all converged, upload, submit, waite for finish, and check again
        current_cycle = 1
        while All_converged != True and current_cycle <= max_cycles:
            logging.info('\nCurrent cycle: %d, not all mols are converged, starting for recheck at  %s\n' % (current_cycle, time.ctime()))

            # upload files, submit jobs, and waite for done
            self.JobControl_upload_submit_and_waite_for_done(conf)

            # check again
            All_converged = self.Check_and_recal(conf, gen_inp_for_planB, plan_B, picked_mols_list)

            # current cycle number
            current_cycle += 1

        if All_converged == True and current_cycle <= max_cycles:
            logging.info('Step %s is finished at %s.' % (conf['Step'], time.ctime()))
        else:
            logging.info('Max cycles are reached for step %s.' % conf['Step'])



if __name__ == '__main__':
    # Each calculation should have an independent conf dict.
    # One can generate or modify the conf dict using scripts.
    # In this part, only a template is required.

    # controls which project to work on
    #from Project_conf_M062X import *
    # from Project_conf_M062X_Bondi import *
    from Project_conf_M062X_SAS_Alpha0485 import *

    #### For step 1 ####
    # combine different confs
    conf = {}
    conf.update(step1_qm_conf)
    conf.update(resource_conf)
    conf.update(project_conf)
    conf.update(project_conf_sub1)

    PrepareJC = PrepareAndJobControl()
    #Analyze = CalAndPlot()

    # pick molecules to deal with
    #picked_mols = ['Hg_l_CH3COO_r__l_H2O_r_+', 'P_Hg_l_CH3COO_r__l_H2O_r_+']

    #conf.update({'Reaction_dataset': 'Gsovl_MNSolv_anion.txt'})
    #conf.update({'Reaction_dataset': 'Gsovl_MNSolv_cation.txt'})
    #conf.update({'Reaction_dataset': 'Gsovl_MNSolv_mol.txt'})

    #picked_mols = ['reaction']
    picked_mols = ['Cys_H0', 'Cys_H1_C', 'Cys_H1_N', 'Cys_H1_R', 'Cys_H2_C', 'Cys_H2_N', 'Cys_H2_R', 'Cys_H3']

    ## Init for step 1
    #PrepareJC.InitNewCal(conf, picked_mols_list=['reaction'])
    PrepareJC.InitNewCal(conf, deleteDir=True, picked_mols_list=picked_mols)

    ## Check and recal for step 1
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=['reaction'])
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=picked_mols)
    PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)

    ## Complement for step 1
    #PrepareJC.InitNewCal(conf, picked_mols_list=['42_a.-1'])
    #PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)

    ## cruise for step 1
    #PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=['reaction'])
    #PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=picked_mols)


    #### For step 2 ####
    conf.update(step2_qm_conf)

    ## Init for step 2
    #PrepareJC.InitNewCal(conf, picked_mols_list=['reaction'])
    #PrepareJC.InitNewCal(conf, picked_mols_list=picked_mols)
    #PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)

    ## Check for step 2
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=['reaction'])
    #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=picked_mols)

    ## Complement for step 2
    #PrepareJC.InitNewCal(conf, picked_mols_list=['42_a.-1'])
    #PrepareJC.Rsync_local_to_remote(conf)
    #PrepareJC.Submit_jobs(conf)

    ## cruise for step 2
    #PrepareJC.Cruise_for_one_step(conf, max_cycles=1, initNew=True, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=['reaction'])
    #PrepareJC.Cruise_for_one_step(conf, max_cycles=1, initNew=True, gen_inp_for_planB=False, plan_B=step1_plan_B, picked_mols_list=picked_mols)



    #### For gas phase correction of OH- ####
    Correction = 0
    if Correction:
        conf.update(step1_gas)

        #picked_mols = ['OH_l_g_r_-', 'Hg_l_g_r_+2', 'Cl_l_g_r_-']
        picked_mols = ['citric_l_g_r_', 'nacitric_l_g_r_',  'sa_l_g_r_', 'urea_l_g_r_' ]

        PrepareJC.InitNewCal(conf, deleteDir=True, picked_mols_list=picked_mols)

        ## Check and recal for step 1
        #PrepareJC.Check_and_recal(conf, gen_inp_for_planB=True, plan_B=step1_plan_B, picked_mols_list=picked_mols)
        PrepareJC.Rsync_local_to_remote(conf)
        PrepareJC.Submit_jobs(conf)

        ## Complement for step 1
        #PrepareJC.InitNewCal(conf, picked_mols_list=['42_a.-1'])
        #PrepareJC.Rsync_local_to_remote(conf)
        #PrepareJC.Submit_jobs(conf)

        ## cruise for step 1
        #PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, gen_inp_for_planB=True, plan_B=step1_plan_B,
        #                    picked_mols_list=['Hg_l_OH_r_3-'])
        #PrepareJC.Cruise_for_one_step(conf, max_cycles=5, initNew=True, deleteDir=True, gen_inp_for_planB=True,
        #                    plan_B=step1_plan_B, picked_mols_list=picked_mols)
