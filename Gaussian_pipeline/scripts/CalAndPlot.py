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
from copy import deepcopy
import pybel
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from sklearn import linear_model

from CSearchRandom import CSearchRand


class CalAndPlot:
    def __init__(self):
        logging.basicConfig(level=logging.INFO)

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

    def read_reactions(self, conf):

        reactions_file = '%s/%s/%s' % (conf['XYZ_folderpath'], conf['Reaction_foldername'], conf['Reaction_dataset'])
        try:
            logKorB = conf['Reaction_dataset'].split('_')[0]
        except:
            logKorB = 'logK'

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

        return logKorB, reactions, reactants_bashsafe, products_bashsafe

    def read_conformations(self, conf):
        # input molecules, in .xyz format
        xyz_folder = '%s/%s/%s' % (conf['Local_pwd'], conf['XYZ_folderpath'], conf['XYZ_foldername'])
        cpds = sorted(os.listdir(xyz_folder), key=lambda x: int(x.split('_')[4]))
        compands_bashsafe = [i[:-4] for i in cpds if i.endswith('.xyz')]

        return compands_bashsafe



    def get_energies_no_corrections(self, conf, mol, folder_path):
        # path for the mol
        mol_path = '%s/%s' % (folder_path, mol)
        if not os.path.exists(mol_path):
            logging.warning('\nWARN: %s is not exist!' % mol_path)

        # all log files for this mol, for conf['Step']
        log_files = [i for i in os.listdir(mol_path) if i[-4:] == '.log' and i.startswith('step%s_' % conf['Step'])]
        if len(log_files) == 0:
            logging.warning('WARN: no log files were found for %s.' % mol)

        # check the last log
        i_log = len(log_files)
        log_file = '%s/step%s_%s_%s.log' % (mol_path, conf['Step'], i_log, mol)
        try:
            return_from_checkfail = g09checkResults.checkNormalTermination(log_file)
        except:
            logging.warn('WARN: failed to check if the calculation terminate normally for file %s' % log_file)
            return_from_checkfail =False

        if return_from_checkfail != True:
            logging.error('ERROR: the calculation %s not finished normally.\n' % log_file)
            return []

        # get the energies
        #print 'Getting final energies from %s ...' % log_file
        Energies = g09checkResults.finalE(open(log_file).readlines())

        return Energies


    def get_energies(self, conf, mol, folder_path, correctProton=True, correctOH=True, correctCl=False, correctHg=False):
        # proton correction
        if mol in ['H+'] and correctProton:
            return [0, 0, 0, 0, 0, 0, 0, 0, -270.30/constants.h2kcal, 0]

        # OH- correction
        if mol in ['OH-'] and correctOH:
            # get therm correction to Gibbs Free Energy in aq phase
            Energies = self.get_energies_no_corrections(conf, mol, folder_path)
            dGthermCorr = Energies[4]

            # get electronic energy at gas phase (Last SCF Done energy) for OH-. (OH_l_g_r_-)
            conf_tmp = deepcopy(conf)
            conf_tmp.update({'Step': '1'})
            Energies_gas = self.get_energies_no_corrections(conf_tmp, 'OH_l_g_r_-', folder_path)
            Escf = Energies_gas[0]

            # correct the result with experimental value
            dGsolv_corr = dGthermCorr + Escf + -104.7/constants.h2kcal

            return [0, 0, 0, 0, 0, 0, 0, 0, dGsolv_corr, 0]

        # Cl- correction
        if mol in ['Cl-'] and correctCl:
            # get therm correction to Gibbs Free Energy in aq phase
            Energies = self.get_energies_no_corrections(conf, mol, folder_path)
            dGthermCorr = Energies[4]

            # get electronic energy at gas phase (Last SCF Done energy) for Cl-. (Cl_l_g_r_-)
            conf_tmp = deepcopy(conf)
            conf_tmp.update({'Step': '1'})
            Energies_gas = self.get_energies_no_corrections(conf_tmp, 'Cl_l_g_r_-', folder_path)
            Escf = Energies_gas[0]

            # correct the result with experimental value
            dGsolv_corr = dGthermCorr + Escf + -74.5/constants.h2kcal

            return [0, 0, 0, 0, 0, 0, 0, 0, dGsolv_corr, 0]

        # Hg+2 correction
        if mol in ['Hg+2'] and correctHg:
            # get therm correction to Gibbs Free Energy in aq phase
            Energies = self.get_energies_no_corrections(conf, mol, folder_path)
            dGthermCorr = Energies[4]

            # get electronic energy at gas phase (Last SCF Done energy) for Hg+2. (Hg_l_g_r_+2)
            conf_tmp = deepcopy(conf)
            conf_tmp.update({'Step': '1'})
            Energies_gas = self.get_energies_no_corrections(conf_tmp, 'Hg_l_g_r_+2', folder_path)
            Escf = Energies_gas[0]

            # correct the result with experimental value
            dGsolv_corr = dGthermCorr + Escf + -446.2/constants.h2kcal

            return [0, 0, 0, 0, 0, 0, 0, 0, dGsolv_corr, 0]


        # no correction. If conformation search using the lowest energy

        ##### to use the lowest structure from bolzman weighting ####
        use_lowest_bolzman_conformation = True
        if use_lowest_bolzman_conformation:
            # check if there is conformation searched structures
            bolz_dir = '../%s/%s/bolzmann/%s' % (conf['Local_output_folder_name'], conf['Group_name'], mol)
            pre_dir = '../%s/%s/' % (conf['Local_output_folder_name'], conf['Group_name'])
            try:
                bolz_conformations = sorted(os.listdir('%s/xyz' % bolz_dir),
                                            key=lambda x: int(x.split('_l_Xe_r_')[1].split('_')[0]))
            except:
                bolz_conformations = []

            # if there is, used the smallest energy one
            if len(bolz_conformations) > 0:
                # for xyz the label 0 means the one with smallest energy,
                # but, for figures have to lookup the csv file to find the 'old frame' number
                bolz_mol = bolz_conformations[0][:-4]
                df_bolz = pd.DataFrame.from_csv('%s/%s.Bolzmann.csv' % (bolz_dir, mol))
                bolz_mol_newFrame = int(bolz_mol.split('_l_Xe_r_')[1].split('_')[0])
                bolz_mol_oldFrame = list(df_bolz[df_bolz.Frame == int(bolz_mol_newFrame)]['Frame_old'])[0]
                bolz_mol_old = bolz_mol.replace('_l_Xe_r_%d_' % bolz_mol_newFrame, '_l_Xe_r_%d_' % bolz_mol_oldFrame)

                # update the xyz
                pre_xyz = '%s/xyz/%s.xyz' % (pre_dir, mol)
                cur_xyz = '%s/xyz/%s.xyz' % (bolz_dir, bolz_mol)
                shutil.copy(cur_xyz, pre_xyz)
                logging.info('\nget_energies: Copying %s to %s' % (cur_xyz, pre_xyz))

                # update the figures
                bolz_figs = [i for i in os.listdir('%s/fig' % pre_dir) if i.startswith('%s_' % bolz_mol_old)]
                for fig in bolz_figs:
                    pre_fig = '%s/fig/%s' % (pre_dir, fig)
                    cur_fig = '%s/fig/%s' % (pre_dir, fig.replace(bolz_mol_old, mol))

                    try:
                        shutil.copy(pre_fig, cur_fig)
                        logging.info('\nget_energies: Copying %s to %s' % (pre_fig, cur_fig))
                    except:
                        logging.warn('\nWARN: Failed to copy figure %s to %s\n' % (pre_fig, cur_fig))
                        print '\nWARN: Failed to copy figure %s to %s\n' % (pre_fig, cur_fig)
                        pass

                mol = bolz_mol_old

        Energies = self.get_energies_no_corrections(conf, mol, folder_path)

        return Energies




    def logK_formula(self, Gaq_M, Gaq_L, Gaq_ML):
        '''
        M + L <-> ML   K

        lnK = -(deltaG)/RT

        :param Gaq_M:  in a.u. Value from 'Sum of electronic and thermal Free Energies'
        :param Gaq_L: in a.u. Value from 'Sum of electronic and thermal Free Energies'
        :param Gaq_ML: in a.u.

        :return in kcal/mol.
        '''
        return -1*((Gaq_ML - Gaq_M - Gaq_L)*constants.h2kcal)/(2.303 * constants.RT)


    def get_num_of_molecules_from_reaction(self, reaction, wat):
        # num of H2O molecules
        try:
            num_reac_h2o = reaction.reac[wat]
        except KeyError:
            num_reac_h2o = 0

        try:
            num_prod_h2o = reaction.prod[wat]
        except KeyError:
            num_prod_h2o = 0

        return {'Num_reac_mol': num_reac_h2o, 'Num_prod_mol': num_prod_h2o}

    def logK_formula_reaction(self, reactants, products, reaction, waterCorrection=True):
        '''
        M + L <-> ML   K

        lnK = -(deltaG)/RT

        :param reactants:  a list of dG in a.u. (Value from 'Sum of electronic and thermal Free Energies') for each reactant
        :param products: a list of dG in a.u. (Value from 'Sum of electronic and thermal Free Energies') for each product
        :param reaction: the reaction obj from chempy for this reaction
        :param waterCorrection: switch on the correction on water

        :return logK in kcal/mol.
        '''

        logK = -1*((sum(products) - sum(reactants)) * constants.h2kcal)/(2.303 * constants.RT)

        if waterCorrection:
            # correction for water molecule and water cluster
            WAT_dict = {'H2O': 1, 'WAT2': 2, 'WAT3': 3, 'WAT4': 4, 'WAT5': 5, 'WAT6': 6}

            dGcorr_WAT = 0.
            for w in WAT_dict.keys():
                Num_WAT = self.get_num_of_molecules_from_reaction(reaction, w)
                dGcorr_WAT += (Num_WAT['Num_prod_mol'] - Num_WAT['Num_reac_mol']) * (constants.RT * math.log(55.34/WAT_dict[w]))

            return logK + dGcorr_WAT/(2.303 * constants.RT)
        else:
            return logK


    def pKa_formula(self, Gaq_A, Gaq_AH, Gaq_H=-270.30):
        '''

        AH <-> A + H K

        pKa = (deltaG)/RT

        :param Gaq_A:  in a.u. Value from 'Sum of electronic and thermal Free Energies'
        :param Gaq_AH: in a.u. Value from 'Sum of electronic and thermal Free Energies'
        :param Gaq_H: in kcal/mol.
        '''
        return ((Gaq_A - Gaq_AH)*constants.h2kcal + Gaq_H)/(2.303 * constants.RT)

    def pKa_formula_reaction(self, reactants, products, reaction, waterCorrection=True):
        '''
        AH <-> A + H   Ka

        pKa = (deltaG)/RT

        :param reactants:  a list of dG in a.u. (Value from 'Sum of electronic and thermal Free Energies') for each reactant
        :param products: a list of dG in a.u. (Value from 'Sum of electronic and thermal Free Energies') for each product

        :return pKa in kcal/mol.
        '''
        pKa = ((sum(products) - sum(reactants)) * constants.h2kcal)/(2.303 * constants.RT)

        if waterCorrection:
            # correction for water molecule and water cluster
            WAT_dict = {'H2O': 1, 'WAT2': 2, 'WAT3': 3, 'WAT4': 4, 'WAT5': 5, 'WAT6': 6}

            dGcorr_WAT = 0.
            for w in WAT_dict.keys():
                Num_WAT = self.get_num_of_molecules_from_reaction(reaction, w)
                dGcorr_WAT += (Num_WAT['Num_prod_mol'] - Num_WAT['Num_reac_mol']) * (
                constants.RT * math.log(55.34 / WAT_dict[w]))

            return pKa + dGcorr_WAT/(2.303 * constants.RT)
        else:
            return pKa


    def cal_logK_for_a_reaction(self, conf, reaction, logKorB='logK', silence=True, calpKa=False, save=True):
        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        # get reactants and products
        reactants = reaction.reac.keys()
        products = reaction.prod.keys()

        # bashsafe
        #reactants_bsafe = [i.replace('(', '_l_').replace(')', '_r_') for i in reactants]
        #products_bsafe = [i.replace('(', '_l_').replace(')', '_r_') for i in products]
        reactants_bsafe = {i: i.replace('(', '_l_').replace(')', '_r_') for i in reactants}
        products_bsafe = {i: i.replace('(', '_l_').replace(')', '_r_') for i in products}


        # get energies for both reactants and products
        #En_reactants = [self.get_energies(conf, i, folder_path) for i in reactants_bsafe]
        #En_products = [self.get_energies(conf, i, folder_path) for i in products_bsafe]

        if conf['Reaction_dataset'] in ['Gsolv_bench.txt']:
            En_reactants = {i: self.get_energies(conf, reactants_bsafe[i], folder_path, correctProton=False, correctOH=False, correctCl=False, correctHg=False) for i in reactants}
            En_products = {i: self.get_energies(conf, products_bsafe[i], folder_path, correctProton=False, correctOH=False, correctCl=False, correctHg=False) for i in products}
        else:
            En_reactants = {i: self.get_energies(conf, reactants_bsafe[i], folder_path, correctProton=True, correctOH=True, correctCl=False, correctHg=False) for i in reactants}
            En_products = {i: self.get_energies(conf, products_bsafe[i], folder_path, correctProton=True, correctOH=True, correctCl=False, correctHg=False) for i in products}




        # please refer g09checkResults.finalE for sequence of energies
        # For DFT it is (finalE, cor_z, cor_u, cor_h, cor_g, z, u, h, g, len(steps))
        # For MP2 it is (finalE, cor_z, cor_u, cor_h, cor_g, z, u, h, g, mp2E_tot, mp2E_cor, len(steps))
        if sum([len(En_reactants[i]) for i in reactants]) == 10 * len(reactants) \
                and sum([len(En_products[i]) for i in products]) == 10 * len(products):
            #logK = self.logK_formula_reaction([en[8] for en in En_reactants], [en[8] for en in En_products], reaction)

            # gen the deltaG energy list for both reactants and products
            dG_reactants = []
            for r in reactants:
                dG_reactants += [En_reactants[r][8]] * reaction.reac[r]

            dG_products = []
            for p in products:
                dG_products += [En_products[p][8]] * reaction.prod[p]

            # calculate the logK/pKa
            if calpKa:
                pKa = self.pKa_formula_reaction(dG_reactants, dG_products, reaction, waterCorrection=True)
            else:
                logK = self.logK_formula_reaction(dG_reactants, dG_products, reaction, waterCorrection=True)


        else:
            logging.error('ERROR: Check the energy items for reactants and products')
            # show details for reactants
            counter = 0
            while counter < len(reactants):
                logging.info('Reactant: %s, %s, %s' % (reactants[counter], reactants_bsafe[reactants[counter]], json.dumps(En_reactants[reactants[counter]])))
                counter += 1
            # show details for products
            counter = 0
            while counter < len(products):
                logging.info('Products: %s, %s, %s' % (products[counter], products_bsafe[products[counter]], json.dumps(En_products[products[counter]])))
                counter += 1

            # error value
            if calpKa:
                pKa = 0
            else:
                logK = 0


        # grab and format the results
        output_series = pd.Series()

        if reaction.param == None:
            if calpKa:
                result_str = 'pKa = %.3f, Exp. = %s, Diff = %s' % (pKa, str(reaction.param), str(reaction.param))
                output_series['Constant'] = 'pKa'
                output_series['Calculated'] = pKa
                output_series['Experimental'] = str(reaction.param)
                output_series['Difference'] = str(reaction.param)
            else:
                result_str = 'logK = %.3f, Exp. = %s, Diff = %s' % (logK, str(reaction.param), str(reaction.param))
                output_series['Constant'] = logKorB
                output_series['Calculated'] = logK
                output_series['Experimental'] = str(reaction.param)
                output_series['Difference'] = str(reaction.param)
        else:
            if calpKa:
                result_str = 'pKa = %.3f, Exp. = %.3f, Diff. = %.3f' % (
                pKa, float(reaction.param), float(pKa - reaction.param))
                output_series['Constant'] = 'pKa'
                output_series['Calculated'] = pKa
                output_series['Experimental'] = float(reaction.param)
                output_series['Difference'] = float(pKa - reaction.param)
            else:
                result_str = 'logK = %.3f, Exp. = %.3f, Diff. = %.3f' % (
                logK, float(reaction.param), float(logK - reaction.param))
                output_series['Constant'] = logKorB
                output_series['Calculated'] = logK
                output_series['Experimental'] = float(reaction.param)
                output_series['Difference'] = float(logK - reaction.param)


        output_str = '#%s\nReaction: %s; %s\n' % ('-' * 50, reaction.string(), result_str)
        output_series['Reaction'] = reaction.string()

        # show details for reactants
        counter = 0
        reactant_str = ''
        reactants_dict = {}
        reactants_bsafe_dict = {}
        while counter < len(reactants):
            reactant_str += 'Reactant: %s, %s\n' % (reactants[counter], str(En_reactants[reactants[counter]][8]))
            reactants_dict[reactants[counter]] = En_reactants[reactants[counter]][8]
            reactants_bsafe_dict[reactants_bsafe[reactants[counter]]] = En_reactants[reactants[counter]][8]
            counter += 1
        output_series['Reactants'] = json.dumps(reactants_dict)
        output_series['ReactantsBSafe'] = json.dumps(reactants_bsafe_dict)

        # show details for products
        counter = 0
        product_str = ''
        products_dict = {}
        products_bsafe_dict = {}
        while counter < len(products):
            product_str += 'Products: %s, %s\n' % (products[counter], str(En_products[products[counter]][8]))
            products_dict[products[counter]] = En_products[products[counter]][8]
            products_bsafe_dict[products_bsafe[products[counter]]] = En_products[products[counter]][8]
            counter += 1
        output_series['Products'] = json.dumps(products_dict)
        output_series['ProductsBSafe'] = json.dumps(products_bsafe_dict)

        output_str = '%s\n%s%s' % (output_str, reactant_str, product_str)

        # add the energy difference between products and reactant
        output_series['deltaG'] = (sum(dG_products) - sum(dG_reactants))*constants.h2kcal

        # show the results
        if not silence:
            logging.info(output_str)

        # save the results
        if save:
            output_dir = '../%s' % conf['Local_output_folder_name']
            output_file = '%s/%s' % (output_dir, conf['Reaction_dataset'])
            try:
                os.mkdir(output_dir)
            except:
                pass

            # save txt
            with open(output_file, 'a') as fout:
                fout.write(output_str)
                fout.close()

        # return the results
        return output_series

    def Cal_logK_for_a_set_of_reactions(self, conf):
        logKorB, reactions, reactants_bashsafe, products_bashsafe = self.read_reactions(conf)

        df = pd.DataFrame()
        # cal the logK
        for reaction in reactions:
            try:
                dfi = self.cal_logK_for_a_reaction(conf, reaction, logKorB=logKorB, silence=False, calpKa=False)
                df = df.append(dfi, ignore_index=True)
                # self.cal_logK_for_a_reaction(conf, reaction, logKorB=logKorB, silence=False, calpKa=True)
            except:
                logging.error('ERROR: Failed to calculate logK for the reaction %s' % reaction.string())

        # save the results
        output_dir = '../%s/%s/csv' % (conf['Local_output_folder_name'], conf['Group_name'])
        # output_dir = '../%s' % conf['Local_output_folder_name']
        output_file = '%s/%s.csv' % (output_dir, conf['Reaction_dataset'][:-4])
        try:
            os.makedirs(output_dir)
        except:
            pass
        df.to_csv(output_file, encoding='utf8')
        return


    def Plot_linear_regression_dG_explogK(self, conf):
        # datafiles
        output_dir = '../%s/%s/csv' % (conf['Local_output_folder_name'], conf['Group_name'])
        csv_file = '%s/%s.csv' % (output_dir, conf['Reaction_dataset'][:-4])
        fig_file = '%s/../fig/linear_dG_explogK_%s.png' % (output_dir, conf['Reaction_dataset'][:-4])

        # read the csv file
        df = pd.DataFrame.from_csv(csv_file)
        if 'deltaG' not in df.columns:
            self.Cal_logK_for_a_set_of_reactions(conf)
            df = pd.DataFrame.from_csv(csv_file)

        # prepare data
        X = df['deltaG'].values.reshape(len(df), 1)
        y = df['Experimental'].values

        # linear regression
        # Create linear regression object
        regr = linear_model.LinearRegression()
        # Train the model using the training sets
        regr.fit(X, y)
        # mean squared error
        mse = np.mean((regr.predict(X) - y) ** 2)

        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.scatter(X, y, color='black')
        ax.plot(X, regr.predict(X), color='blue', linewidth=2)
        #ax.set_xlabel('deltaG (kcal/mol)')
        #ax.set_ylabel('pKa (experiment)')
        ax.set_xlabel('Calculation')
        ax.set_ylabel('Experiment')
        ax.set_title('Linear Regression (deltaG - pKa)')

        text_str = 'y = %.4fx + %.4f\nR**2 = %.4f\nMSE = %.2f' % (regr.coef_[0], regr.intercept_, regr.score(X, y), mse)
        if regr.coef_[0] < 0:
            text_x = 0.8
        else:
            text_x = 0.2

        ax.text(text_x, 0.8, text_str, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        # save fig
        fig.savefig(fig_file)

        return

    def Plot_linear_regression_callogK_explogK(self, conf):
        # datafiles
        output_dir = '../%s/%s/csv' % (conf['Local_output_folder_name'], conf['Group_name'])
        csv_file = '%s/%s.csv' % (output_dir, conf['Reaction_dataset'][:-4])
        fig_file = '%s/../fig/linear_callogK_explogK_%s.png' % (output_dir, conf['Reaction_dataset'][:-4])

        # read the csv file
        df = pd.DataFrame.from_csv(csv_file)

        # prepare data
        X = df['Calculated'].values.reshape(len(df), 1)
        y = df['Experimental'].values

        # linear regression
        # Create linear regression object
        regr = linear_model.LinearRegression()
        # Train the model using the training sets
        regr.fit(X, y)
        # mean squared error
        mse = np.mean((regr.predict(X) - y) ** 2)

        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        max_X = max(X)
        max_y = max(y)
        max_Xory = max([max_X, max_y])

        ax.scatter(X, y, color='black')
        #ax.plot(X, regr.predict(X), color='red', linewidth=2)
        ax.plot(np.linspace(0, max_Xory, 20), np.linspace(0, max_Xory, 20), color='red', linewidth=2)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Calculation')
        ax.set_ylabel('Experiment')
        ax.set_title('Linear Regression (logK/pKa - logK/pKa)')

        text_str = 'y = %.4fx + %.4f\nR**2 = %.4f\nMSE = %.2f' % (regr.coef_[0], regr.intercept_, regr.score(X, y), mse)
        if regr.coef_[0] < 0:
            text_x = 0.8
        else:
            text_x = 0.2

        #ax.text(text_x, 0.8, text_str, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        # save fig
        fig.savefig(fig_file)

        return

    def Plot_linear_regression_caldG_expdG(self, conf):
        # datafiles
        output_dir = '../%s/%s/csv' % (conf['Local_output_folder_name'], conf['Group_name'])
        csv_file = '%s/%s.csv' % (output_dir, conf['Reaction_dataset'][:-4])
        fig_file = '%s/../fig/linear_caldG_expdG_%s.png' % (output_dir, conf['Reaction_dataset'][:-4])

        # read the csv file
        df = pd.DataFrame.from_csv(csv_file)

        # prepare data
        X = df['deltaG'].values.reshape(len(df), 1)
        y = df['Experimental'].values

        # linear regression
        # Create linear regression object
        regr = linear_model.LinearRegression()
        # Train the model using the training sets
        regr.fit(X, y)
        # mean squared error
        mse = np.mean((regr.predict(X) - y) ** 2) ### it's not mean signed error !!!

        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        max_X = max(X)
        max_y = max(y)
        max_Xory = max([max_X, max_y])
        min_X = min(X)
        min_y = min(y)
        min_Xory = min([min_X, min_y])

        ax.scatter(X, y, color='black')
        ax.plot(X, regr.predict(X), color='blue', linewidth=2)
        ax.plot(np.linspace(min_Xory, max_Xory, 20), np.linspace(min_Xory, max_Xory, 20), color='red', linewidth=2)
        #ax.set_xlim(left=0)
        #ax.set_ylim(bottom=0)
        ax.set_xlabel('deltaG')
        ax.set_ylabel('Experiment')
        ax.set_title('Linear Regression (deltaG)')

        text_str = 'y = %.4fx + %.4f\nR**2 = %.4f\nFittingError = %.2f' % (regr.coef_[0], regr.intercept_, regr.score(X, y), np.sqrt(mse))
        if regr.coef_[0] < 0:
            text_x = 0.8
        else:
            text_x = 0.2

        ax.text(text_x, 0.8, text_str, color='blue', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        # save fig
        fig.savefig(fig_file)

        return

    def Plot_linear_regression_charge_explogK(self, conf):
        # datafiles
        output_dir = '../%s/%s/csv' % (conf['Local_output_folder_name'], conf['Group_name'])
        csv_file = '%s/%s.csv' % (output_dir, conf['Reaction_dataset'][:-4])
        fig_file = '%s/../fig/linear_charge_explogK_%s.png' % (output_dir, conf['Reaction_dataset'][:-4])

        # read the csv file
        df = pd.DataFrame.from_csv(csv_file)

        # get all reactants for each reaction
        reactants = [eval(i).keys() for i in df.Reactants]

        # get the 'Laxx-' reactant
        reactants_list = []
        for i in reactants:
            for r in i:
                if r.startswith('L'):
                    reactants_list.append(r)

        # get charge
        def get_charge(reactant, type):
            '''
            type should in the bellow
            [u'No_', u'Atom', u'Mulliken', u'APT', u'Natural', u'Core', u'Valence', u'Rydberg', u'Total']
            '''

            dfi = pd.DataFrame.from_csv('%s/../charge/%s.csv' % (output_dir, reactant))
            charge = min(list(dfi[dfi.Atom == 'S'][type]))
            return charge

        natural = np.array([get_charge(r, 'Natural') for r in reactants_list])


        # prepare data
        X = natural.reshape(len(natural), 1)
        y = df['Experimental'].values

        # linear regression
        # Create linear regression object
        regr = linear_model.LinearRegression()
        # Train the model using the training sets
        regr.fit(X, y)
        # mean squared error
        mse = np.mean((regr.predict(X) - y) ** 2)

        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.scatter(X, y, color='black')
        ax.plot(X, regr.predict(X), color='red', linewidth=2)
        ax.set_xlabel('Natural Charge')
        ax.set_ylabel('pKa (experiment)')
        ax.set_title('Linear Regression (Natural charge - pKa)')

        text_str = 'y = %.4fx + %.4f\nR**2 = %.4f\nMSE = %.2f' % (regr.coef_[0], regr.intercept_, regr.score(X, y), mse)
        if regr.coef_[0] < 0:
            text_x = 0.8
        else:
            text_x = 0.2

        ax.text(text_x, 0.8, text_str, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        # save fig
        fig.savefig(fig_file)

        return



    def plot_result_for_a_compound(self, conf, cpd, overwrite=True):
        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        # find the command
        g09pyGauss = '%s/g09pyGauss.py' % os.path.dirname(inspect.getfile(constants))

        # plot
        logging.info('%s\nPlot for %s' % ('-' * 50, cpd))
        folder = '%s/%s' % (folder_path, cpd)

        # if cpd is samplinged, use the bolzmann sampling version.
        bolz_dir = '../%s/%s/bolzmann/%s' % (conf['Local_output_folder_name'], conf['Group_name'], cpd)
        if os.path.exists(bolz_dir):
            self.get_energies(conf, cpd, folder_path)
            return

        plot_cmds = '%s all %s %s' % (g09pyGauss, str(overwrite), folder)
        plotp = subprocess.Popen(plot_cmds, shell=True).wait()


        # move ALL figures to ../output/figures/
        #output_dir = '../%s_figures' % conf['Local_output_folder_name']
        output_dir = '../%s/%s/fig' % (conf['Local_output_folder_name'], conf['Group_name'])
        try:
            os.makedirs(output_dir)
        except:
            pass

        # move
        for f in [i for i in os.listdir(folder) if i.endswith('.png')]:
            dst_file = '%s/%s' % (output_dir, f)
            if os.path.exists(dst_file):
                if overwrite:
                    os.remove(dst_file)
                else:
                    return

            shutil.move('%s/%s' % (folder, f), dst_file)

        return

    def Plot_results_for_a_set_of_reactions(self, conf, overwrite=True):
        logKorB, reactions, reactants_bashsafe, products_bashsafe = self.read_reactions(conf)

        # plot the results
        for cpd in reactants_bashsafe + products_bashsafe:
            try:
                self.plot_result_for_a_compound(conf, cpd, overwrite=overwrite)
            except:
                logging.error('ERROR: Failed to plot result for the compound %s' % cpd)

        return


    def Plot_results_for_a_set_of_conformations(self, conf, overwrite=True):
        compounds_bashsafe = self.read_conformations(conf)

        # plot the results
        for cpd in compounds_bashsafe:
            try:
                self.plot_result_for_a_compound(conf, cpd, overwrite=overwrite)
            except:
                logging.error('ERROR: Failed to plot result for the compound %s' % cpd)

        return



    def get_xyz_for_a_compound(self, conf, cpd, logfile_step=2, overwrite=True):
        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        # input/output
        folder = '%s/%s' % (folder_path, cpd)
        output_dir = '../%s/%s/xyz' % (conf['Local_output_folder_name'], conf['Group_name'])
        fout = '%s/%s.xyz' % (output_dir, cpd)

        try:
            os.makedirs(output_dir)
        except:
            pass


        #### for bolzmann sampling structures
        bolz_dir = '../%s/%s/bolzmann/%s' % (conf['Local_output_folder_name'], conf['Group_name'], cpd)
        pre_dir = '../%s/%s/' % (conf['Local_output_folder_name'], conf['Group_name'])
        try:
            bolz_conformations = sorted(os.listdir('%s/xyz' % bolz_dir),
                                        key=lambda x: int(x.split('_l_Xe_r_')[1].split('_')[0]))
        except:
            bolz_conformations = []

        # check Bolzmann first!
        # if there is, used the smallest energy one
        if len(bolz_conformations) > 0:
            bolz_mol = bolz_conformations[0][:-4]

            # update the xyz
            pre_xyz = '%s/xyz/%s.xyz' % (pre_dir, cpd)
            cur_xyz = '%s/xyz/%s.xyz' % (bolz_dir, bolz_mol)
            shutil.copy(cur_xyz, pre_xyz)

            logging.info('Using bolzimann structure.\nCopying %s to %s' % (cur_xyz, pre_xyz))

            return
        else:
            # read log files
            try:
                logfiles = [i for i in os.listdir(folder) if i.endswith('.log') and i.startswith('step%s_' % str(logfile_step))]
            except:
                logging.error('Cannot find log files for %s for Step %s' % (folder, str(logfile_step)))
                return

            mol = pybel.readfile('g09', '%s/%s' % (folder, logfiles[-1])).next()

            # write xyz
            mol.write('xyz', fout, overwrite=overwrite)

        return

    def get_charge_for_a_compound(self, conf, cpd, logfile_step=2, overwrite=True):
        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        # input/output
        folder = '%s/%s' % (folder_path, cpd)
        output_dir = '../%s/%s/charge' % (conf['Local_output_folder_name'], conf['Group_name'])
        fout = '%s/%s.csv' % (output_dir, cpd)


        #### for bolzmann sampling structures
        bolz_dir = '../%s/%s/bolzmann/%s' % (conf['Local_output_folder_name'], conf['Group_name'], cpd)
        pre_dir = '../%s/%s/' % (conf['Local_output_folder_name'], conf['Group_name'])
        try:
            bolz_conformations = os.listdir('%s/xyz' % bolz_dir)
        except:
            bolz_conformations = []

        # check Bolzmann first!
        # if there is, used the one with lowest energy
        if len(bolz_conformations) > 0:
            bolz_mol = bolz_conformations[0][:-4]
            folder = '%s/%s' % (folder_path, bolz_mol)

            # update the logfiles
            try:
                logfiles = [i for i in os.listdir('%s' % folder)
                            if i.endswith('.log') and i.startswith('step%s_' % str(logfile_step))]
            except:
                logging.error('Cannot find log files for %s for Step %s' % (folder, str(logfile_step)))
                return
        else:
            #### for normal structures

            # read log files
            try:
                logfiles = [i for i in os.listdir(folder) if i.endswith('.log') and i.startswith('step%s_' % str(logfile_step))]
            except:
                logging.error('Cannot find log files for %s for Step %s' % (folder, str(logfile_step)))
                return

        try:
            os.makedirs(output_dir)
        except:
            pass


        fcon = open('%s/%s' % (folder, logfiles[-1])).readlines()
        mul_charge = g09checkResults.getMullikenCharge(fcon)
        apt_charge = g09checkResults.getAPTCharge(fcon)
        natural_charge = g09checkResults.getNaturalPop(fcon)

        df_mul = pd.DataFrame(mul_charge, columns=['No_', 'Atom', 'Mulliken'])
        df_apt = pd.DataFrame(apt_charge, columns=['No_', 'Atom', 'APT'])
        df_natural = pd.DataFrame(natural_charge, columns=['Atom', 'No_', 'Natural', 'Core', 'Valence', 'Rydberg', 'Total'])

        df = pd.concat([df_mul, df_apt['APT'], df_natural[range(2, 7)]], axis=1)

        # save the result
        df.to_csv(fout)

        return

    def Get_xyz_for_a_set_of_reactions(self, conf):
        logKorB, reactions, reactants_bashsafe, products_bashsafe = self.read_reactions(conf)

        # get the results
        for cpd in reactants_bashsafe + products_bashsafe:
            # Deal with the OH(g)-, Cl(g)-, Hg(g)+2 correction
            mol = ''
            if cpd in ['OH-']:
                mol = 'OH_l_g_r_-'
            elif cpd in ['Cl-']:
                mol = 'Cl_l_g_r_-'
            elif cpd in ['Hg+2']:
                mol = 'Hg_l_g_r_+2'

            # get xyz for OH(g)-, Cl(g)-, Hg(g)+2 correction
            if mol != '':
                try:
                    self.get_xyz_for_a_compound(conf, mol, logfile_step=1, overwrite=True)
                except:
                    logging.error('ERROR: Failed to get xyz for %s' % mol)

            # get xyz & charge pop for cpd
            try:
                self.get_xyz_for_a_compound(conf, cpd, logfile_step=2, overwrite=True)
                self.get_charge_for_a_compound(conf, cpd, logfile_step=2)
            except:
                logging.error('ERROR: Failed to get xyz for %s' % cpd)

    def Get_xyz_for_a_set_of_conformations(self, conf):
        compounds_bashsafe = self.read_conformations(conf)

        # get the results
        for cpd in compounds_bashsafe:
            # Deal with the OH(g)-, Cl(g)-, Hg(g)+2 correction
            mol = ''
            if cpd in ['OH-']:
                mol = 'OH_l_g_r_-'
            elif cpd in ['Cl-']:
                mol = 'Cl_l_g_r_-'
            elif cpd in ['Hg+2']:
                mol = 'Hg_l_g_r_+2'

            # get xyz for OH(g)-, Cl(g)-, Hg(g)+2 correction
            if mol != '':
                try:
                    self.get_xyz_for_a_compound(conf, mol, logfile_step=1, overwrite=True)
                except:
                    logging.error('ERROR: Failed to get xyz for %s' % mol)

            # get xyz for cpd
            try:
                self.get_xyz_for_a_compound(conf, cpd, logfile_step=2, overwrite=True)
                self.get_charge_for_a_compound(conf, cpd, logfile_step=2)
            except:
                logging.error('ERROR: Failed to get xyz for %s' % cpd)






    def calBolzmannWeightedEN(self, inCsv, outCsv):
        '''

        :param inCsv: Energy unit is hatree!
        :return:  Energy unit is hatree!
        '''


        # read inCsv
        df = pd.DataFrame.from_csv(inCsv)

        # cal deltaEn
        df['deltaEn_kcalmol-1'] = (df['En'] - min(df['En'])) * constants.h2kcal

        # cal ratio of each conformation to the lowest one
        df['Ratio2Lowest'] = np.exp(-(df['deltaEn_kcalmol-1']) / (constants.kb_kcalmolK * constants.T))

        # cal energy contribution of each conformation to the Bolzmann weighted energy
        df['En_fraction'] = df['En'] * (df['Ratio2Lowest'] / df['Ratio2Lowest'].sum())

        # Bolzmann weighted energy
        En_weighted = df['En_fraction'].sum()

        # save the csv
        df.to_csv(outCsv)

        return En_weighted

    def Bolzmann_weighting(self, conf, eps=0.03, minSamples=1):
        cs = CSearchRand()

        # validate the conf first
        conf = g09prepare.validate_conf(conf)

        # prepare vars
        folder_name, folder_path, conf_path, xyz_folder, submit_all_file = self.gen_fundamental_vars(conf)

        # compounds from conf['XYZ_foldername']
        compounds_bashsafe = self.read_conformations(conf)

        # determine mol name
        mol = xyz_folder.split('/')[-2]

        # input/output
        xyz_dir = '../%s/%s/xyz' % (conf['Local_output_folder_name'], conf['Group_name'])
        bolz_dir = '../%s/%s/bolzmann/%s' % (conf['Local_output_folder_name'], conf['Group_name'], mol)

        try:
            os.makedirs(bolz_dir)
        except:
            pass

        # variables
        pdb_en = '%s/%s.opted.en' % (bolz_dir, mol)
        pdb_combine = '%s/%s.opted.pdb' % (bolz_dir, mol)
        pdb_rms = '%s/%s.opted.rms' % (bolz_dir, mol)

        clusterFig = '%s/%s.cluster.png' % (bolz_dir, mol)
        clusterCsv = '%s/%s.cluster.csv' % (bolz_dir, mol)

        uniqueCsv = '%s/%s.unique.csv' % (bolz_dir, mol)
        uniquePDB = '%s/%s.unique.pdb' % (bolz_dir, mol)

        bolzCsv = '%s/%s.Bolzmann.csv' % (bolz_dir, mol)

        # get free energy for all conformations
        ##Ens = [self.get_energies(conf, i, folder_path)[8] for i in self.read_conformations(conf)]
        #Ens = [self.get_energies(conf, i, folder_path)[8] for i in compounds_bashsafe]
        #Frames = [int(i.split('_')[4]) for i in compounds_bashsafe]
        compounds_bashsafe_updated = []
        Ens = []
        Frames = []

        for c in compounds_bashsafe:
            try:
                i_En = self.get_energies(conf, c, folder_path)[8]
            except:
                i_En = 0.0

            i_Frames = int(c.split('_')[4])

            print c, i_En, i_Frames

            compounds_bashsafe_updated.append(c)
            Ens.append(i_En)
            Frames.append(i_Frames)



        # combine Ens with compound names
        cpd_Ens = zip(compounds_bashsafe_updated, Ens, Frames)
        cpd_Ens_unsorted = np.array(zip(compounds_bashsafe_updated, Ens, Frames))
        cpd_Ens = np.array(sorted(cpd_Ens, key=lambda x: x[1]))

        # output en
        df_tmp = pd.DataFrame(cpd_Ens[:, [1, 2]])
        #open(pdb_en, 'w').writelines("\n".join(cpd_Ens[:, 1]))
        df_tmp.to_csv(pdb_en, header=False, index=False)

        # convert xyz files to one pdb file
        fout_pdb = open(pdb_combine, 'w')
        # structures will be picked from xxx.opted.pdb by its sequence.
        # so, don't sort the write sequence of the structures here.
        for cpd in cpd_Ens_unsorted[:, 0]:
            xyz_file = '%s/%s.xyz' % (xyz_dir, cpd)
            m = pybel.readfile('xyz', xyz_file).next()
            fout_pdb.writelines(m.write('pdb'))
        fout_pdb.close()

        # gene rms
        cs.calcRMS(pdb_combine, pdb_rms, debug=False)

        # plot cluster
        cs.plotCluster(pdb_en, pdb_rms, clusterFig, clusterCsv, eps=eps, minSamples=minSamples, debug=False)

        # pick structures
        print clusterCsv
        print pdb_combine
        print uniqueCsv
        print uniquePDB
        cs.uniqueClusters(clusterCsv, pdb_combine, uniqueCsv, uniquePDB, debug=False)

        # cal Bolzmann weighted energy
        En_weighed = self.calBolzmannWeightedEN(uniqueCsv, bolzCsv)
        logging.info('The Bolzmann weighted energy for %s is %s' % (mol, str(En_weighed)))

        # split pdb
        cs.splitePDBtoXYZ(uniquePDB, '%s/xyz/' % bolz_dir, cs.outPrefix, mol, debug=False)

        return En_weighed




if __name__ == '__main__':
    # Each calculation should have an independent conf dict.
    # One can generate or modify the conf dict using scripts.
    # In this part, only a template is required.

    # controls which project to work on
    from Project_conf_M062X_SAS_Alpha0485 import *

    project_list = [
        'Project_conf_M062X_SAS_Alpha0485',
        'Project_conf_M062X_Bondi',
        'Project_conf_M062X',

    ]



    for m in project_list:
        exec('from %s import *' % m)

        #### For conf build ####
        # combine different confs
        conf = {}
        conf.update(step1_qm_conf)
        conf.update(resource_conf)
        conf.update(project_conf)
        conf.update(project_conf_sub1)

        ######### For calculate the pKa #########
        conf.update(step2_qm_conf)

        # reaction dataset to deal with. (Will update that from Project_conf.py)
        txt_list = ['Gsovl_MNSolv_anion.txt']
        txt_list = ['Gsovl_MNSolv_cation.txt']
        #txt_list = ['Gsovl_MNSolv_mol.txt']


        for txt in txt_list:
            conf.update({'Reaction_dataset': txt})
            Analyze = CalAndPlot()
            #### For calculating the logK ####
            Analyze.Cal_logK_for_a_set_of_reactions(conf)

            #### For linear regression ####
            #Analyze.Plot_linear_regression_dG_explogK(conf)
            #Analyze.Plot_linear_regression_charge_explogK(conf)
            #Analyze.Plot_linear_regression_callogK_explogK(conf)
            Analyze.Plot_linear_regression_caldG_expdG(conf)


            if txt in [
                       #'logK_Haworth_w0.txt',
                       #'Gsovl_MNSolv_anion.txt',
                       ]:
                logging.info('Ploting %s for %s' % (m, txt))
                #### For plotting the results ####
                Analyze.Plot_results_for_a_set_of_reactions(conf, overwrite=True)
                #### For getting xyz ####

            #Analyze.Get_xyz_for_a_set_of_reactions(conf)
        '''

        ######### For Conformation sampling & Bolzmann weighting #########
        #molecules = ['Hg_l_SCH2CH2COOH_r_+', 'Hg_l_SCH2CH2COO_r_', 'Hg_l_SCH2CH2COO_r_2-2', 'Hg_l_SCHCH2COO2_r_-',
        #         'Hg_l_SCHCH2COOH2_r_+', 'SCHCH2COO2-3', 'SCHCH2COOH2-', 'SHCHCH2COO2-2', 'SHCHCH2COOH2']
        #molecules += ['Hg_l_SCHCH2COO2_r_2-4']

        #molecules = ['Hg_l_CH3NH2_r_2+2', 'Hg_l_CH3NH2_r_3+2', 'Hg_l_CH3NH2_r_4+2']

        molecules = ['Hg_l_SCH2CH2COOH_r_+', 'Hg_l_SCH2CH2COO_r_', 'Hg_l_SCH2CH2COO_r_2-2', 'Hg_l_SCHCH2COO2_r_-',
                 'Hg_l_SCHCH2COOH2_r_+', 'SCHCH2COO2-3', 'SCHCH2COOH2-', 'SHCHCH2COO2-2', 'SHCHCH2COOH2']
        molecules += ['Hg_l_SCHCH2COO2_r_2-4']

        molecules = [
                        "CH3CHSCOO-2",
                        "CH3CHSCOOH-",
                        "SCHCH2COO2-3",
                        "SCHCH2COOH2-",
                        "SHCHCH2COO2-2",
                        "SHCHCH2COOH2",
                        "Hg_l_CH3CHSCOOHV1_r_+",
                        "Hg_l_CH3CHSCOO_r_",
                        "Hg_l_CH3CHSCOO_r_2-2",
                        "Hg_l_CH3NH2_r_2+2",
                        "Hg_l_CH3NH2_r_4+2",
                        "Hg_l_CH3NH2_r_3+2",
                        "Hg_l_SCH2CH2COOH_r_+",
                        "Hg_l_SCH2CH2COO_r_",
                        "Hg_l_SCH2CH2COO_r_2-2",
                        "Hg_l_SCHCH2COO2_r_-",
                        "Hg_l_SCHCH2COO2_r_2-4",
                        "Hg_l_SCHCH2COOH2_r_+"
        ]


        Analyze = CalAndPlot()
        for cpd in molecules:
            conf.update({'XYZ_foldername': 'structures/%s/xyz' % cpd})
            #### For plotting the results ####
            Analyze.Plot_results_for_a_set_of_conformations(conf, overwrite=True)
            #### For getting xyz ####
            Analyze.Get_xyz_for_a_set_of_conformations(conf)
            #### For Bolzmann weighting ####
            Analyze.Bolzmann_weighting(conf, eps=0.06, minSamples=1)


        #molecules = ['SCH2CH2COO-2', 'SCH2CH2COOH-', 'SHCH2CH2COOH', 'SHCH2CH2COO-']
        #Analyze = CalAndPlot()
        #for cpd in molecules:
        #    Analyze.plot_result_for_a_compound(conf, cpd)

        '''
        '''
        ######### For Bolzmann weighting #########
        molecules = ['Hg_l_SCH2CH2COOH_r_+', 'Hg_l_SCH2CH2COO_r_', 'Hg_l_SCH2CH2COO_r_2-2', 'Hg_l_SCHCH2COO2_r_-',
                 'Hg_l_SCHCH2COOH2_r_+', 'SCHCH2COO2-3', 'SCHCH2COOH2-', 'SHCHCH2COO2-2', 'SHCHCH2COOH2']

        #molecules = ['Hg_l_CH3NH2_r_2+2', 'Hg_l_CH3NH2_r_3+2', 'Hg_l_CH3NH2_r_4+2']

        molecules += ['Hg_l_SCH2CH2COO_r_']
        for cpd in molecules:
            conf.update({'XYZ_foldername': 'structures/%s/xyz' % cpd})
            Analyze = CalAndPlot()
            Analyze.Bolzmann_weighting(conf, eps=0.06, minSamples=1)
        '''

