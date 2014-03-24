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
Command-line utility to associate each mature miRNA with its precursor.
The files outputted by this program is used as input for various programs in the DecoyDB (working title) project
"""

import os
import argparse
import cPickle

def extract_seq(input_file):
    """ Extract the fasta and convert it as a dict with key=seq_name, value=sequence """
    result_dict = dict()
    with open(input_file, 'rb') as input:
        list_lines = input.readlines()
        i = 0
        while i < len(list_lines):
            line = list_lines[i]
            if line.startswith(">"):
                splitted = line.strip().split()
                name = splitted[0].replace(">", "")
                acc = splitted[1]
                header = line.strip()
                j = i+1
                seq = ""
                while j < len(list_lines):
                    curr_line = list_lines[j]
                    if curr_line.startswith(">"):
                        break
                    else:
                        seq += curr_line.strip()
                        j += 1
                result_dict[acc] = dict(name=name, accession=acc, sequence=seq, header=header)
                i = j
            else:
                i += 1
    return result_dict

def write_output(hairpin_dir, mature_dir, digested_data):
    """ write the separated output file """
    with open(os.path.join(hairpin_dir, digested_data["accession"]), 'wb') as hairpin_out:
        hairpin_out.write("{header}\n{sequence}".format(header=digested_data["header"],
                                                        accession=digested_data["accession"],
                                                        sequence=digested_data["sequence"]))
    
    to_be_written = []
    for mature in digested_data["matures"]:
        to_be_written.append("{header}\n{sequence}".format(header=mature["header"],
                                                           accession=mature["accession"],
                                                           sequence=mature["sequence"]))
    with open(os.path.join(mature_dir, digested_data["accession"]), 'wb') as hairpin_out:
        hairpin_out.write("\n".join(to_be_written))

def add_5p_3p(dict_mature, hairpin_sequence):
    """ Edit if needed the name of the mature to add 5p or 3p at the end """
    if dict_mature['header'].endswith('3p') or dict_mature['header'].endswith('5p'):
        result_dict = dict_mature
    else:
        mature_seq = dict_mature["sequence"]
        mature_start = hairpin_sequence.find(mature_seq)
#         mature_middle = (mature_start + len(mature_seq)/2)
#         if mature_middle > len(hairpin_sequence)/2:
#             result_dict = dict(name=dict_mature["name"]+"-3p",
#                                accession=dict_mature["accession"],
#                                sequence=dict_mature["sequence"])
#         else:
#             result_dict = dict(name=dict_mature["name"]+"-5p",
#                                accession=dict_mature["accession"],
#                                sequence=dict_mature["sequence"])
        mature_end = mature_start + len(mature_seq)
        remaining_length_on_3p = len(hairpin_sequence) - mature_end
        if mature_start < remaining_length_on_3p:
            result_dict = dict(header=dict_mature["header"]+"-5p",
                               alternative_name=dict_mature["name"]+"-5p",
                               accession=dict_mature["accession"],
                               sequence=dict_mature["sequence"],
                               name=dict_mature["name"])
        else:
            result_dict = dict(header=dict_mature["header"]+"-3p",
                               alternative_name=dict_mature["name"]+"-3p",
                               accession=dict_mature["accession"],
                               sequence=dict_mature["sequence"],
                               name=dict_mature["name"])
    return result_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--hairpin_fa', action="store", required=True, dest="hairpin_fa")
    parser.add_argument('--mature_fa', action="store", required=True, dest="mature_fa")
    parser.add_argument('--mirna_dat', action="store", required=True, dest="mirna_dat")
    parser.add_argument('--out_dir', action="store", required=True, dest="out_dir")

    ns = parser.parse_args()

    mirna_dat = ns.mirna_dat
    hairpin_fa = ns.hairpin_fa
    mature_fa = ns.mature_fa
    out_dir = ns.out_dir

    dict_hairpin = extract_seq(hairpin_fa)
    dict_mature = extract_seq(mature_fa)

    # if no hairping and/or mature dir, create it
    hairpin_dir = os.path.join(out_dir, "hairpin")
    mature_dir = os.path.join(out_dir, "mature")
    if not os.path.exists(hairpin_dir):
        os.mkdir(hairpin_dir)
    if not os.path.exists(mature_dir):
        os.mkdir(mature_dir)

    #the digested_data will be put into a list of dict(accession, name, sequence, matures=[dict(accession, name, sequence)])
    list_digested_data = []

    # associated each pre-microRNA with its mature(s) 
    with open(mirna_dat, 'rb') as mirna_dat_c:
        list_lines = mirna_dat_c.readlines()
        i = 0
        while i < len(list_lines):
            line = list_lines[i]
            # we found the beginning of a section
            if line.startswith("ID"):
                splitted_line = line.split()
                name = splitted_line[1]
                acc = ""
                list_mature = []
                j = i+1
                while j < len(list_lines):
                    curr_line = list_lines[j]
                    if curr_line.startswith("AC"):
                        # we found the accession for the precursor
                        acc = curr_line.split()[1].replace(";", "")
                    elif curr_line.startswith("FT") and "accession" in curr_line:
                        # we found a mature
                        mirna_accession = curr_line.split('"')[1].strip()
                        list_mature.append(mirna_accession)
                    elif curr_line.startswith("ID"):
                        # we're on another pre-miRNA's turf, BAIL OUT!
                        break
                    j += 1

                # only process humans
                if splitted_line[1].startswith("hsa"):
                    digested_data = dict(accession=acc,
                                         name=name,
                                         sequence=dict_hairpin[acc]["sequence"],
                                         header=dict_hairpin[acc]["header"],
                                         matures=[add_5p_3p(dict_mature[elem], dict_hairpin[acc]["sequence"]) for elem in list_mature])

                    write_output(hairpin_dir, mature_dir, digested_data)

                    list_digested_data.append(digested_data)

                i = j
            else:
                i += 1

    # save the data as a pickle
    with open(os.path.join(out_dir, "digested_data.pk"), 'wb') as digested_data_pk:
        cPickle.dump(list_digested_data, digested_data_pk, -1)