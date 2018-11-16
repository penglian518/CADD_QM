#!/usr/bin/env python
#
# @purpose
#   to start a new calculation environment.
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 20 2016
#

import os, argparse, logging, json


class LogAnalyze:
    def __init__(self):
        logging.basicConfig(level=logging.INFO)

    def find_imaginary_freq_structures(self, logFile):
        '''
        This function is to find the structures with imaginary freqs

        :param conf:
        :param logFile:
        :return: a list of structures
        '''

        lines = open(logFile).readlines()
        Structures_img = list(set([i.split()[-2].split('/')[-1] for i in lines if i.rfind('imaginary frequency') > 0]))
        Cal_method = list(set([i.split()[-2].split('/')[-2] for i in lines if i.rfind('imaginary frequency') > 0]))

        Structures_unfinished = list(set([i.split()[-1].strip('.') for i in lines if i.rfind('no log files') > 0]))
        Structures_unfinished += list(set([i.split()[-4].split('/')[-1] for i in lines if i.rfind('Cannot find log files') > 0]))
        Structures_unfinished += list(set([i.split()[-1].split('/')[-1] for i in lines if i.rfind('Failed to get xyz') > 0]))

        Mol_not_exist = list(set([i.split()[-4].split('/')[-1] for i in lines if i.rfind('is not exist') > 0]))

        print '''
Calculation model:
    %s

Structures with Imaginary Frequencies:
    %s

Structures with No log files, or Cannot get xyz:
    %s

Molecules not exist:
%s

All improper structures:
%s


''' % (json.dumps(Cal_method), json.dumps(Structures_img), json.dumps(Structures_unfinished),
       json.dumps(Mol_not_exist), json.dumps(Structures_img + Structures_unfinished, indent=4))


        return Cal_method, Structures_img, Structures_unfinished



if __name__ == '__main__':

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Find the structures with imaginary frequencies.')
    parser.add_argument("type", nargs='?', type=str, default="freq", help="Type for the search. Default: freq")
    parser.add_argument("logfile", type=str, default="", help="Input the path to the log file")

    # all arguments
    args = parser.parse_args()
    Analyzer = LogAnalyze()

    # check the folder
    if not os.path.exists(args.logfile):
        print 'Log file %s NOT exist!' % args.folder
        exit()

    # perform the code
    if args.type in ['f', 'freq', 'frequency', 'img']:
        Analyzer.find_imaginary_freq_structures(args.logfile)




