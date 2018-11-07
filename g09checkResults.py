#!/usr/bin/env python
#
# @purpose
#   to check the output for g09
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 28 2016
#


from sys import argv, exit
from glob import glob
from os import path, listdir
from math import ceil
#import numpy as np

## To do:
## build in MP2 functionality.

# Defining colors
class color:
   prp = '\033[95m'
   cyn = '\033[96m'
   dkcyn = '\033[36m'
   blu = '\033[94m'
   grn = '\033[92m'
   ylo = '\033[93m'
   red = '\033[91m'
   blk = '\033[99m'
   bld = '\033[1m'
   uln = '\033[4m'
   end = '\033[0m'


## Defining functions
def checkfail(file, silence=False):
    o = open(file, 'r').readlines()
    olen = len(o)
    fail = 'n'

    reason = ''
    warn = ''
    n = olen-10

    if n < 0:
        warn = color.prp + file + ' log file is empty.' + color.end
        if not silence:
            print warn
        fail = 'y'
        reason = 'empt'
    else:
        while n < olen-1:
            if 'Convergence failure -- run terminated.' in o[n]:
                warn = color.prp + file + ' geometry optimization failed to converge.' + color.end
                if not silence:
                    print warn
                fail = 'y'
                reason = 'conv'
                break
            elif 'The combination of multiplicity' in o[n]:
                warn =  color.prp + file + ' failed because of incorrect charge and multiplicty.' + color.end
                if not silence:
                    print warn
                fail = 'y'
                reason = 'uerr'
                break
            elif 'Error termination request processed by link 9999.' in o[n]:
                warn = color.prp + file + ' reached max wallclock time.'+ color.end
                if not silence:
                    print warn
                fail = 'y'
                reason = 'time'
                break
            elif 'Error termination via Lnk1e' in o[n]:
                warn =  color.prp + file + ' failed. Check output file.'+ color.end
                if not silence:
                    print warn
                fail = 'y'
                reason = 'unk'
                break
            else:
                n += 1

    if fail != 'y':
        return 'pass'
    else:
        if silence:
            return 'fail', reason, warn
        else:
            return 'fail', reason

    #o.close()

def checkNormalTermination(file):
    try:
        fcon = open(file).readlines()
    except IOError:
        print 'File %s is not exist.' % file
        return False
    return fcon[-1].startswith(' Normal termination of Gaussian 09')

def getSCFEnergies(file):
    try:
        fcon = open(file).readlines()
    except IOError:
        print 'File %s is not exist.' % file
        return False

    scf_lines = [i.strip() for i in fcon if i.startswith(' SCF Done:')]
    SCF_Ens = [float(i.split('=')[1].split()[0]) for i in scf_lines]

    return SCF_Ens

def getSCFEnergies_from_fcon(fcon):
    scf_lines = [i.strip() for i in fcon if i.startswith(' SCF Done:')]
    SCF_Ens = [float(i.split('=')[1].split()[0]) for i in scf_lines]

    return SCF_Ens


def getChargeAndSpin(fcon):
    charge_lines = [i for i in fcon if i.startswith(' Charge =')]
    charge = charge_lines[0].split()[2]
    multiplicity = charge_lines[0].split()[5]

    return (int(charge), int(multiplicity))

def getNAtom(fcon):
    natom_lines = [i for i in fcon if i.startswith(' NAtoms= ')]
    natom = natom_lines[0].split()[1]

    return int(natom)

def getMullikenCharge(fcon):
    natom = getNAtom(fcon)

    # get index for the flag line
    counter = 0
    idx = 0
    for i in fcon:
        if i.startswith(' Mulliken charges:'):
            idx = counter
        counter += 1

    rawtext = fcon[idx+2 : idx+2+natom]

    # [['1', 'S', '-0.901663'], ..., ['No.', 'Atom', 'Charge']]
    charges = [i.strip().split() for i in rawtext]

    return charges

def getAPTCharge(fcon):
    natom = getNAtom(fcon)

    # get index for the flag line
    counter = 0
    idx = 0
    for i in fcon:
        if i.startswith(' APT charges:'):
            idx = counter
        counter += 1

    rawtext = fcon[idx+2 : idx+2+natom]

    # [['1', 'S', '-0.901663'], ..., ['No.', 'Atom', 'Charge']]
    charges = [i.strip().split() for i in rawtext]

    return charges

def getNaturalPop(fcon):
    natom = getNAtom(fcon)

    # get index for the flag line
    counter = 0
    idx = 0
    for i in fcon:
        if i.startswith(' Summary of Natural Population Analysis:'):
            idx = counter
        counter += 1

    rawtext = fcon[idx+6 : idx+6+natom]

    # [['Atom', 'No', 'Natural Charge', 'Core', 'Valence', 'Rydberg', 'Total']]
    charges = [i.strip().split() for i in rawtext]

    return charges





# Fetching file stats
def getstats(file):
    optcompl = 'n'
    completed = 'n'
    proc = '0'
    nodes = '1'
    date = 'Running'
    step = '1'
    for i, line in enumerate(file):
        if 'Gaussian 09:' in line:
            vers = line.split()[2].split('-')[1]
        elif '%chk=' in line.lower():
            chk = line.split('=')[1].strip()
        try:
            if '%nprocshared' in line.lower() or '%nproc' in line.lower():
                proc = line.split('=')[1]
        except:
                proc = '0'
        try:
            if '%nproclinda' in line.lower():
                nodes = line.split('=')[1]
        except:
            nodes = '0'
        if '%mem' in line.lower().strip():
            mem = line.split('=')[1]
            x = 0
            keywds = []
            while x < 10:
                if '--------' in file[i + x]:
                    y = 1
                    while y < 10-x:
                        if '--------' in file[i + x + y]:
                            break
                        else:
                            keywds.append(file[i + x + y].strip())
                            y += 1
                    break
                else:
                    x += 1
            keywds = ''.join(keywds)
            j = 0
            for k in keywds.split():
                if k.split('=')[0].lower() == 'opt':
                    job = 'opt'
                    j = 1
                elif k.split('=')[0].lower() == 'freq':
                    job = 'freq'
                    j = 1
            if j == 0:
                job = 'sp'
        elif 'Charge =' in line:
            Q = line.split()[2]
            M = line.split()[5]
        elif 'Framework group' in line:
            stoi = line.split()[2].split('[')[-1].strip(']')
            ptgp = line.split()[2].split('[')[0]
        try:
            if 'Stoichiometry' in line:
                stoi = line.split()[1]
            elif 'Full point group' in line:
                ptgp = line.split()[3]
        except:
            continue
        if 'NAtoms=' in line:
            natoms = int(line.split()[1])
        elif 'basis functions,' in line:
            nbasis = line.split()[0]
        elif 'Step number' in line:
            step = line.split()[2]
        elif 'Stationary point found' in line:
            optcompl = 'y'
        elif 'Normal termination' in line:
            completed = 'y'
            date = ' '.join(line.split()[-5:])
        else:
            continue
    return vers, chk, proc, nodes, mem, keywds, job, natoms, Q, M, nbasis, stoi, ptgp, optcompl, completed, date, step

def getfreq(file):
    freqs = []
    ifreqs = []
    count = len(file)-1
    maxcount = file[count]
    while count > 0:
        if file[count].find('Frequencies --') != -1:
            for n in file[count].split()[2:]:
                if float(n) < 0:
                    ifreqs.append(float(n))
                    freqs.append(float(n))
                else:
                    freqs.append(float(n))
        count -= 1
#    print freqs
    return freqs, ifreqs

def get_lowvib(freqs, cutoff):
    low_freqs = []
    for n in sorted(freqs):
        if n < cutoff:
                low_freqs.append(n)
    mat = []
    r = 0
    while r < ceil(len(low_freqs) / 6.0):
        c = 0
        row = []
        while c < 6:
            try:
                low_freqs[c + 6 * r]
            except:
                break
            row.append(low_freqs[c + 6 * r])
            c += 1
        mat.append(row)
        r += 1
    return len(low_freqs), mat

# Fetch convergence
def getconv(file):
    conv_field = []
    conv_value = []
    conv_thresh = []
    conv_step = []
    count = 0
    conv = 'n'
    maxcount = len(file)-1
    while count <= maxcount:
        field = []
        value = []
        thresh = []
        step = []
        if file[count].find('Maximum Force') != -1:
            conv = 'y'
            j = 0
            max_j = 3
            while j <= max_j:
                field.append(str(file[count + j].split()[0] + ' ' + file[count +j].split()[1]))
                value.append(file[count +j].split()[2])
                thresh.append(file[count +j].split()[3])
                if file[count +j].split()[4] == 'NO':
                    step.append(color.red + '[NO]' + color.end)
                else:
                    step.append(color.blu + '[YES]' + color.end)
                j += 1
            conv_value.append(value)
            conv_thresh.append(thresh)
            conv_step.append(step)
        count += 1
    return field, conv_value, conv_thresh, conv_step, conv

# Fetching final energy
def finalE(file):
    count = len(file)-1
    maxcount = file[count]
    z = 0
    u = 0
    h = 0
    g = 0
    cor_z = 0
    cor_u = 0
    cor_h = 0
    cor_g = 0
    mp2E_tot = 0
    mp2E_cor = 0
    mp2if = 'n'
    finalE = 0
    SCF = getSCFEnergies_from_fcon(file)
    steps = []
    while count > 0:
        if file[count].find('Final MP2 energy:') != -1:
            mp2E = float(file[count].split()[:-1])
            mp2if = 'y'
        if file[count].find('Sum of electronic and thermal Free Energies=') != -1:
            g = float(file[count].split()[7])
        elif file[count].find('Sum of electronic and thermal Enthalpies=') != -1:
            h = float(file[count].split()[6])
        elif file[count].find('Sum of electronic and thermal Energies=') != -1:
            u = float(file[count].split()[6])
        elif file[count].find('Sum of electronic and zero-point Energies=') != -1:
            z = float(file[count].split()[6])
        if file[count].find('Thermal correction to Gibbs Free Energy=') != -1:
            cor_g = float(file[count].split()[6])
        elif file[count].find('Thermal correction to Enthalpy=') != -1:
            cor_h = float(file[count].split()[4])
        elif file[count].find('Thermal correction to Energy=') != -1:
            cor_u = float(file[count].split()[4])
        elif file[count].find('Zero-point correction=') != -1:
            cor_z = float(file[count].split()[2])

        elif file[count].find('Optimization completed.') != -1:
            steps.append(1)
        count -= 1
    if mp2if == 'y':
        finalE = SCF[-1] - mp2E_cor 
        return finalE, cor_z, cor_u, cor_h, cor_g, z, u, h, g, mp2E_tot, mp2E_cor, len(steps)
    else:
        finalE = SCF[-1]
        return finalE, cor_z, cor_u, cor_h, cor_g, z, u, h, g, len(steps)

def tail(file, length):
    tail = []
    count = len(file)-1-length
    maxcount = len(file)-1
    while count < maxcount:
        tail.append(file[count].strip())
        count += 1
    return tail

def bulk_output(file):
    o = open(file, 'r').readlines()
    stats = getstats(o)
#    print stats
    if stats[14] == 'n':
        f = color.grn + file + color.end
    elif stats[14] == 'y' and stats[6] == 'opt' and stats[13] == 'n':
        f = color.grn + file + color.end
    elif stats[14] == 'y' and len(getfreq(o)[1]) != 0:
        f = color.red + file + color.end
    else:
        f = color.blk + file + color.end

#    print stats[6]

    # Fetch frequency data
    if stats[6] == 'freq':# and len(getfreq(o)[0]) != 0:
        try:
            freq = sorted(getfreq(o)[0])[0]
            nofq = 'n'
        except:
            print color.grn + file + color.end
            return
    else:
        freq = '-'
        nofq = 'y'

    E = finalE(o)
    if stats[6] == 'freq' and E[1] != 0:
        ZPE = E[5]
        U = E[6]
        H = E[7]
        G = E[8]
    else:
        ZPE = '-'
        U = '-'
        H = '-'
        G = '-'

    if nofq == 'y':
        print "%-99s %13.6f" % (f, E[0])
    else:
        print "%-99s %13.6f %13.5f %13.5f %13.5f %13.6f %13.2f" % (f, E[0], ZPE, U, H, G, freq)

def single_output(file):
    o = open(file, 'r').readlines()
    fail = []
    fail.append([n for n in checkfail(file)])
    if len(fail) == 2 and fail[1] in ('unk', 'uerr'):
        print '-------------------------------------------------------------------------------'
        for n in tail(o, 9):
            print n
        print color.dkcyn + '-------------------------------------------------------------------------------' + color.end
        exit()

    try:
        stats = getstats(o)
    #                 vers, chk, proc, nodes, keywds, job, natoms, Q, M, nbasis, stoi, optcompl, completed, time
    #                 0         1        2            3             4             5        6                  7     8    9                10        11                 12                    13
    except:
        print file, 'error or you checked too soon after job began.'
        print '-------------------------------------------------------------------------------'
        for n in tail(o, 9):
            print n
        print color.dkcyn + '-------------------------------------------------------------------------------' + color.end
        exit()

    E = finalE(o)
    fq = getfreq(o)
    cv = getconv(o)
    cutoff = 100

    print file
    print color.dkcyn + '---------------------------------------------------------------' + color.end
    print "%-8s %9s %13s %-25s" % ('Version:', stats[0], 'Checkpoint:', stats[1])
    print "%-6s %-2s %7s %-2s %6s %-24s" % ('Nodes:', int(stats[3]), 'Cores:', int(stats[2]), 'Date:', stats[15])
    print '---------------------------------------------------------------'
    print "%-105s" % (stats[5])
    print "%10s %-4s %9s %-4s %-15s %-12s" % ('No. Atoms:', stats[7], 'Point Group:', stats[12], 'Stoichiometry:', stats[11])
    print "%-8s %-2s %14s %-2s %12s %-5s" % ('Charge:', stats[8], 'Multiplicity:', stats[9], 'Basis functions:', stats[10])
    print '---------------------------------------------------------------'
    if cv[4] == None and stats[6] == 'opt':
        for n in tail(o, 5):
            print n
    if cv[4] == 'y':
        print 'Step No.:'.ljust(5), (str(E[-1])+'.'+stats[16]).rjust(5),'Max'.rjust(6),'|','value'.ljust(6),'RMS'.rjust(8),'|','value'.ljust(6)
        print 'Force:'.ljust(13), cv[3][-1:][0][0].rjust(17), '|', cv[1][-1:][0][0].ljust(10), cv[3][-1:][0][1].rjust(8), '|', cv[1][-1:][0][1].ljust(10)
        print 'Displacement:'.ljust(13), cv[3][-1:][0][2].rjust(17), '|', cv[1][-1:][0][2].ljust(10), cv[3][-1:][0][3].rjust(8), '|', cv[1][-1:][0][3].ljust(10)
        print '---------------------------------------------------------------'
    if len(E) > 10:
        print "%-22s %-13.6f %-22s %-13.6f" % ('Final SCF Energy:', E[0], 'Final MP2 Energy', E[9])
    else:
        if E[0] == 0:
            print "%-22s %-13.6f %-22s" % ('Final SCF Energy:', 0.000000, 'The first SCF calculation has not finished.')
            print color.prp + '"All human wisdom is summed up in two words; wait and hope." -Alexandre Dumas' + color.end
        else:
            print "%-22s %-13.6f" % ('Final SCF Energy:', E[0])
    if stats[14] == 'n':
        print '---------------------------------------------------------------'
        for n in tail(o, 5):
            print n
    else:
        if stats[6] == 'freq':
            print "%-22s %-13.6f %2s%-8.6f%1s" % ('Zero-point Energy:', E[5], '(', E[1],')')
            print "%-22s %-13.6f %2s%-8.6f%1s" % ('Internal Energy:', E[6], '(', E[2],')')
            print "%-22s %-13.6f %2s%-8.6f%1s" % ('Enthalpy:', E[7], '(', E[3],')')
            print "%-22s %-13.6f %2s%-8.6f%1s" % ('Gibbs Free Energy:', E[8], '(', E[4],')')
            print '---------------------------------------------------------------'

            if len(fq[1]) == 1:
                print color.grn + 'WARNING: There is 1 imaginary vibrational frequency.' + color.end
            elif len(fq[1]) > 1:
                print color.grn + 'WARNING: There are ' + str(len(fq[1])) + ' imaginary vibrational frequencies.' + color.end
            print str(get_lowvib(fq[0], cutoff)[0]) + ' vibrational frequencies below ' + str(cutoff) + ' cm-1:'
            #mat = np.asarray(get_lowvib(fq[0], cutoff)[1][0])
            mat = get_lowvib(fq[0], cutoff)[1]
            #print mat
            print ('\n'.join([''.join(['{0:9}'.format(item) for item in row]) for row in mat]))




def deal_with_argv():
    ## Interpreting user inputs. Setting variables and paths.
    if len(argv) > 4:
        print "Usage: Script.py"
        exit()

    wdir = '.'
    foi = ''
    bulk = 'n'

    ## Search logic
    search = '-a'
    below = 'n'
    if len(argv) == 1:
        wdir = path.dirname('./')
        bulk = 'y'
    elif len(argv) == 2:
        if argv[1] == '-below':
            wdir = path.dirname('./')
            below = 'y'
            bulk = 'y'
        elif path.isfile(argv[1]) and argv[1].split('.')[-1] == 'log':
            if len(path.split(argv[1])) == 1:
                wdir = path.dirname('./')
                foi = argv[1]
                bulk = 'n'
            elif len(path.split(argv[1])) == 2:
                wdir = path.dirname(argv[1])
                foi = path.split(argv[1])[1]
                bulk = 'n'
        elif path.isfile(argv[1]) and 'log' not in argv[1].split('.')[-1]:
            print color.grn + argv[1], 'is not a valid G09 .log file.' + color.end
            exit()
        elif path.isdir(argv[1]):
            wdir = argv[1]
            bulk = 'y'
        elif path.isdir(argv[1]) is False and path.isfile(argv[1]) is False:
            wdir = path.dirname('./')
            search = argv[1]
            bulk = 'y'
    elif len(argv) == 3:
        if path.isdir(argv[1]) and argv[2] != '-below':
            wdir = argv[1]
            search = argv[2]
            bulk = 'y'
        elif path.isdir(argv[1]) and argv[2] == '-below':
            wdir = argv[1]
            bulk = 'y'
            below = 'y'
        elif path.isdir(argv[1]) is False and argv[2] == '-below':
            wdir = path.dirname('./')
            search = argv[1]
            below = 'y'
            bulk = 'y'
    elif len(argv) == 4:
        if path.isdir(argv[1]) and path.isfile(argv[2]) is False and argv[3] == '-below':
            wdir = argv[1]
            search = argv[2]
            bulk = 'y'
            below = 'y'

    return wdir, foi, search, bulk, below

def check_Results(wdir, foi, search, bulk, below):
    ## Searching through directories immediately beneath the working dir for *opt.logs of specified type.
    if bulk == 'y':
        print "%-90s %13s %13s %13s %13s %13s %13s" % ('File', 'Final SCF', 'ZPE', 'U', 'H', 'G', 'Lowest Freq')

        if below == 'y':
            j = 0
            for item in listdir(wdir):
              if path.isdir(path.join(wdir, item)):
                    wd = path.join(wdir, item)
                    j += 1
                    k = 0
                    for file in sorted(glob(path.join(wd, '*.log'))):
                        if search == '-a':
                            if checkfail(file) == 'pass':
                                bulk_output(file)
                            k += 1
                        elif search in file:
                            if checkfail(file) == 'pass':
                                bulk_output(file)
                            k += 1

                    try:
                        k
                    except:
                        k = 0
                    if k == 0:
                        print 'No .log files found within ' + wd + '/'
            if j == 0:
                try:
                    wd
                except:
                    wd = wdir
                print 'No subdirectories found within ' + wd + '/'
        else:
            k = 0
            for file in sorted(glob(path.join(wdir, '*.log'))):
                if search == '-a':
                    if checkfail(file) == 'pass':
                        bulk_output(file)
                    k += 1
                elif search in file:
                    if checkfail(file) == 'pass':
                        bulk_output(file)
                    k += 1

            try:
                k
            except:
                k = 0
            if k == 0:
                print 'No .log files found within ' + wdir + '/'

    if bulk == 'n':
        single_output(path.join(wdir, foi))


if __name__ == '__main__':
    wdir, foi, search, bulk, below = deal_with_argv()
    #print wdir, foi, search, bulk, below
    check_Results(wdir, foi, search, bulk, below)