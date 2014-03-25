#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#
# This utility is part of the DecoyDB (working title) project
#
# Author: Stephen Leong Koan
#
# Copyright (C) 2013 Université de Montréal
#

"""
Command-line utility to compute the bp stats for each precursors
"""

import os
import math
import argparse
import subprocess
import shlex
import shutil
import multiprocessing
import cPickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--digested_data',
                        action="store",
                        required=True,
                        dest="digested_data")
    parser.add_argument('--mcfold_output_dir',
                        action="store",
                        required=True,
                        dest='mcfold_output_dir')
    parser.add_argument('--considered_structures',
                        action="store",
                        dest='considered_structures',
                        type=int,
                        default=100)
    parser.add_argument('--out_pk',
                        action="store",
                        required=True,
                        dest="out_pk")

    ns = parser.parse_args()

    digested_data = ns.digested_data
    mcfold_output_dir = ns.mcfold_output_dir
    considered_structures = ns.considered_structures
    threads = ns.threads
    out_pk = ns.out_pk

    digested_data_list = []
    with open(digested_data, 'rb') as dd:
        digested_data_list = cPickle.load(dd)

    list_final_digested = []

    for digested_data_dict in digested_data_list:
        hairpin_name = digested_data_dict['name']
        hairpin_acc = digested_data_dict['accession']
        hairpin_seq = digested_data_dict['sequence']
        mcfold_output = os.path.join(mcfold_output_dir, hairpin_acc)

        if not os.path.exists(mcfold_output):
            print "{dir} not found".format(dir=mcfold_output)
        else:
            with open(mcfold_output, 'rb') as mcf_o:
                if considered_structures == 0:
                    list_struct = [elem.strip().split()[0] for elem in mcf_o.readlines() if elem.strip()]
                else:
                    list_struct = [elem.strip().split()[0] for elem in mcf_o.readlines() if elem.strip()][:considered_structures]

                if list_struct:
                    list_perc = []
                    for i in xrange(len(hairpin_seq)):
                        num_unpaired = len([1 for struct in list_struct if struct[i] == "."])
                        list_perc.append(float(num_unpaired)/float(len(list_struct)))

                    list_final_digested.append(dict(name=hairpin_name,
                                                    accession=hairpin_acc,
                                                    sequence=hairpin_seq,
                                                    list_percentage=list_perc,
                                                    considered_structures=considered_structures))

    with open(out_pk, 'wb') as out_p:
        cPickle.dump(list_final_digested, out_p, -1)