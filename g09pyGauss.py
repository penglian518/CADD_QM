#!/home/p6n/anaconda2/envs/pg_env/bin/python
#
# @Note
# This script requires a separate virtual environment, pg_env, to host pygauss
#
# Have to used numpy=1.9.2 instead of the default 1.9.3, which will give the error _global xxx module!!!
#
# @purpose
#   to analysis and plot from g09 calculations
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Nov 30 2016
#


import os, argparse
import pygauss as pg


def get_mol(path_to_folder):
    '''
    generate mol file and file existance dict for pygauss analysis.

    :param path_to_folder:
    :return:
    '''

    folder_obj = pg.Folder(path_to_folder)

    # collect files in the folder
    files_in_folder = os.listdir(path_to_folder)
    init_files = [i for i in files_in_folder if i.startswith('step1_1_') and i.endswith('.com')]
    opt_files = [i for i in files_in_folder if i.startswith('step1_') and i.endswith('.log')]
    freq_files = [i for i in files_in_folder if i.startswith('step2_') and i.endswith('.log')]

    # file existance. Assume all file exist.
    file_existance = {
        'folder': path_to_folder,
        'folder_name': os.path.basename(path_to_folder),
        'init_files': 1,
        'opt_files': 1,
        'freq_files': 1,
        'nbo_files': 1
        }

    # no init file, return
    if len(init_files) == 0:
        print 'WARN: cannot find init file for gaussian calculation from %s.' % path_to_folder
        mol = pg.molecule.Molecule(folder_obj=folder_obj)
        file_existance['init_files'] = 0
    else:
        file_existance['init_fname'] = init_files[0]

    # no opt file
    if len(opt_files) == 0:
        print 'WARN: cannot find opt file for gaussian calculation from %s.' % path_to_folder
        print 'Will only deal init file %s.' % init_files[0]
        mol = pg.molecule.Molecule(folder_obj=folder_obj, init_fname=init_files[0])
        file_existance['opt_files'] = 0
    else:
        # this is a list of file names
        file_existance['opt_fname'] = opt_files

    # freq files, nbo file
    if len(freq_files) == 0:
        mol = pg.molecule.Molecule(folder_obj=folder_obj, init_fname=init_files[0], opt_fname=opt_files)
        file_existance['freq_files'] = 0
        file_existance['nbo_files'] = 0
    else:
        mol = pg.molecule.Molecule(folder_obj=folder_obj, init_fname=init_files[0], opt_fname=opt_files, freq_fname=freq_files[-1], nbo_fname=freq_files[-1])
        file_existance['freq_fname'] = freq_files[-1]
        file_existance['nbo_fname'] = freq_files[-1]


    return mol, file_existance


def plot_init_structure(path_to_folder, represent='ball_stick', rotations=[[0.0, 0.0, 0.0]], overwrite=True):
    # get mol and file existance
    mol, file_existance = get_mol(path_to_folder)

    if file_existance['init_files'] > 0:
        # seems there is some bugs in tuning the fig size of this show_initial function
        image_obj = mol.show_initial(represent=represent, rotations=rotations, width=640, height=480, zoom=1.0)
        fig_out = '%s/%s_init.%s' % (path_to_folder, file_existance['folder_name'], image_obj.format)

        if os.path.exists(fig_out) and not overwrite:
            #print 'Fig %s is exist, will skip to overwrite' % fig_out
            pass
        else:
            # write the fig
            with open(fig_out, 'w') as fout:
                fout.write(image_obj.data)
                fout.close()
    else:
        print 'There is no init file found from %s.' % path_to_folder

    return

def plot_opt_structure(path_to_folder, represent='ball_stick', rotations=[[0.0, 0.0, 0.0]], overwrite=True):
    # get mol and file existance
    mol, file_existance = get_mol(path_to_folder)

    if file_existance['opt_files'] > 0:
        image_obj = mol.show_optimisation(represent=represent, rotations=rotations, width=640, height=480, zoom=1.0)
        #fig_out = '%s/opt_%s.%s' % (path_to_folder, file_existance['opt_fname'][-1][:-4], image_obj.format)
        fig_out = '%s/%s_opted.%s' % (path_to_folder, file_existance['folder_name'], image_obj.format)

        if os.path.exists(fig_out) and not overwrite:
            #print 'Fig %s is exist, will skip to overwrite' % fig_out
            pass
        else:
            # write the fig
            with open(fig_out, 'w') as fout:
                fout.write(image_obj.data)
                fout.close()
    else:
        print 'There is no opt file found from %s.' % path_to_folder

    return

def plot_nbo_charge(path_to_folder, represent='ball_stick', rotations=[[0.0, 0.0, 0.0]], overwrite=True):
    # get mol and file existance
    mol, file_existance = get_mol(path_to_folder)

    if file_existance['freq_files'] > 0:
        image_obj = mol.show_nbo_charges(represent=represent, rotations=rotations, width=640, height=480, zoom=1.0)
        fig_out = '%s/%s_nboCharge.%s' % (path_to_folder, file_existance['folder_name'], image_obj.format)

        if os.path.exists(fig_out) and not overwrite:
            #print 'Fig %s is exist, will skip to overwrite' % fig_out
            pass
        else:
            # write the fig
            with open(fig_out, 'w') as fout:
                fout.write(image_obj.data)
                fout.close()
    else:
        print 'There is no freq file found from %s.' % path_to_folder

    return

def plot_freq(path_to_folder, overwrite=True):
    # get mol and file existance
    mol, file_existance = get_mol(path_to_folder)

    if file_existance['freq_files'] > 0:
        # get freq
        try:
            fr = mol.get_freq_analysis()
        except:
            fr = ''
            print 'Failed to read freq from the file %s' % path_to_folder
            return fr

        fr_imaginary = fr[fr.ix[:, 0] <= 0]
        if len(fr_imaginary) > 0:
            print 'WARN: %d imaginary frequency found from %s !!!' % (len(fr_imaginary), path_to_folder)

        # plot freq
        ax = mol.plot_freq_analysis()
        image_obj = ax.get_figure()
        fig_out = '%s/%s_freq.png' % (path_to_folder, file_existance['folder_name'])

        if os.path.exists(fig_out) and not overwrite:
            #print 'Fig %s is exist, will skip to overwrite' % fig_out
            pass
        else:
            # write the fig
            image_obj.savefig(fig_out)
    else:
        print 'There is no freq file found from %s.' % path_to_folder

    return fr

def plot_opt_energy(path_to_folder, overwrite=True):
    # get mol and file existance
    mol, file_existance = get_mol(path_to_folder)

    if file_existance['opt_files'] > 0:
        # plot freq
        ax = mol.plot_opt_energy(units='hartree')
        image_obj = ax.get_figure()
        fig_out = '%s/%s_optEN.png' % (path_to_folder, file_existance['folder_name'])

        if os.path.exists(fig_out) and not overwrite:
            #print 'Fig %s is exist, will skip to overwrite' % fig_out
            pass
        else:
            # write the fig
            image_obj.savefig(fig_out)
    else:
        print 'There is no opt file found from %s.' % path_to_folder

    return

def plot_homo_lumo(path_to_folder, overwrite=True):
    # get mol and file existance
    mol, file_existance = get_mol(path_to_folder)

    if file_existance['freq_files'] > 0:
        # plot freq
        ax = mol.plot_dos(per_energy=1, lbound=-20, ubound=10, legend_size=12)
        image_obj = ax.get_figure()
        fig_out = '%s/%s_gap.png' % (path_to_folder, file_existance['folder_name'])

        if os.path.exists(fig_out) and not overwrite:
            #print 'Fig %s is exist, will skip to overwrite' % fig_out
            pass
        else:
            # write the fig
            image_obj.savefig(fig_out)
    else:
        print 'There is no freq file found from %s.' % path_to_folder

    return

def plot_all(path_to_folder, represent='ball_stick', rotations=[[0.0, 0.0, 0.0]], overwrite=True):

    try:
        plot_init_structure(path_to_folder, represent=represent, rotations=rotations, overwrite=overwrite)
    except:
        print 'Failed to plot init structure'

    try:
        plot_opt_structure(path_to_folder, represent=represent, rotations=rotations, overwrite=overwrite)
    except Exception as e:
        print 'Failed to plot opt structure'
        print '{}'.format(e)

    try:
        plot_nbo_charge(path_to_folder, represent=represent, rotations=rotations, overwrite=overwrite)
    except:
        print 'Failed to plot nbo charge'

    try:
        plot_freq(path_to_folder, overwrite=overwrite)
    except:
        print 'Failed to plot freq'

    try:
        plot_opt_energy(path_to_folder, overwrite=overwrite)
    except:
        print 'Failed to plot opt energies'

    try:
        plot_homo_lumo(path_to_folder, overwrite=overwrite)
    except:
        print 'Failed to plot homo lumo gap'

    return


def main():

    return


if __name__ == '__main__':

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Plot the result from Gaussian calculations with pyGauss.')
    parser.add_argument("type", nargs='?', type=str, default="all", help="Type for the plot. Default: plot all")
    parser.add_argument("overwrite", nargs='?', type=str, default="False", help="Overwrite the fig or not. Default: False")
    parser.add_argument("folder", type=str, default="", help="Input the path to the calculation folder")

    # all arguments
    args = parser.parse_args()

    # check the folder
    if not os.path.exists(args.folder):
        print 'Folder %s NOT exist!' % args.folder
        exit()

    # perform the code
    if args.type in ['all', 'a', 'All']:
        if args.overwrite in ['False', False, 0, '0']:
            # plot all
            plot_all(path_to_folder=args.folder, represent='ball_stick', rotations=[[0.0, 0.0, 0.0]], overwrite=False)
        else:
            # plot all
            plot_all(path_to_folder=args.folder, represent='ball_stick', rotations=[[0.0, 0.0, 0.0]], overwrite=True)
    elif args.type in ['freq', 'Freq']:
        freq = plot_freq(path_to_folder=args.folder, overwrite=args.overwrite)
        print freq

