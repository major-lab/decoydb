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
Command-line utility to filter MC-Fold result to make sure there is no stem on the complementary side of each mature
"""

import os
import sys
import math
import argparse
import subprocess
import shlex
import shutil
import multiprocessing

THRESHOLD = 2.0
REFINE_THRESHOLD = THRESHOLD + 1.56
RELIEVE_THRESHOLD = THRESHOLD + 2.0
REGULAR_THRESHOLD = THRESHOLD + 3.0
FOUND = False

# The vectors are precomputed by Naim
VECTOR_DICT_MAPPER = {"13": {"AtoB": 36.5306196771968, "AtoC": 25.3202217407352,
                             "AtoD": 14.5687755834181, "BtoC": 15.5526173038495,
                             "BtoD": 48.5117137813127, "CtoD": 35.2541875101384},
                      "14": {"AtoB": 41.0472377267947, "AtoC": 29.9432734015505,
                             "AtoD": 14.3471857170666, "BtoC": 15.5526173038495,
                             "BtoD": 51.9925513415143, "CtoD": 39.5621360141234},
                      "15": {"AtoB": 44.8856673449332, "AtoC": 33.1136380665127,
                             "AtoD": 15.4871359521378, "BtoC": 15.5526173038495,
                             "BtoD": 56.0769255309169, "CtoD": 43.5565930485845},
                      "16": {"AtoB": 48.3460794377372, "AtoC": 35.5524276808209,
                             "AtoD": 15.396983503271, "BtoC": 15.5526173038495,
                             "BtoD": 59.8666426651771, "CtoD": 46.8568560298277},
                      "17": {"AtoB": 51.4500872788375, "AtoC": 37.6623044170162,
                             "AtoD": 15.1277614008154, "BtoC": 15.5526173038495,
                             "BtoD": 63.3616909496582, "CtoD": 49.9966151754296},
                      "18": {"AtoB": 54.2819544968676, "AtoC": 39.8456138740514,
                             "AtoD": 15.755359849905, "BtoC": 15.5526173038495,
                             "BtoD": 66.9310196246852, "CtoD": 52.9056908942696},
                      "19": {"AtoB": 55.3081796211013, "AtoC": 40.2346122884265,
                             "AtoD": 15.7824418896443, "BtoC": 15.5526173038495,
                             "BtoD": 69.3772206779718, "CtoD": 54.7387508991574},
                      "20": {"AtoB": 57.3817541993969, "AtoC": 42.0469492591318,
                             "AtoD": 14.2715699556846, "BtoC": 15.5526173038495,
                             "BtoD": 70.6406870719701, "CtoD": 55.5776207929055},
                      "21": {"AtoB": 57.3962924325256, "AtoC": 42.2683001550808,
                             "AtoD": 15.0380274637334, "BtoC": 15.5526173038495,
                             "BtoD": 71.9972296342019, "CtoD": 56.7566115091449},
                      "22": {"AtoB": 59.2869801389816, "AtoC": 44.7495173940457,
                             "AtoD": 15.7300000635728, "BtoC": 15.5526173038495,
                             "BtoD": 74.1059640244967, "CtoD": 59.0931359888778},
                      "23": {"AtoB": 61.4383611272306, "AtoC": 47.4458354863733,
                             "AtoD": 15.796888997521, "BtoC": 15.5526173038495,
                             "BtoD": 75.3480757551246, "CtoD": 60.7243238990769},
                      "24": {"AtoB": 65.5493830100635, "AtoC": 52.1481783382699,
                             "AtoD": 15.9282791600348, "BtoC": 15.5526173038495,
                             "BtoD": 78.360371170637, "CtoD": 64.1792743570695},
                      "25": {"AtoB": 69.4461850500083, "AtoC": 56.0654800835594,
                             "AtoD": 16.1551173316692, "BtoC": 15.5526173038495,
                             "BtoD": 82.0196073021567, "CtoD": 68.1810333817258},
                      "26": {"AtoB": 73.4726019275213, "AtoC": 59.9833852079057,
                             "AtoD": 16.2752557276376, "BtoC": 15.5526173038495,
                             "BtoD": 86.0297278561312, "CtoD": 72.2933319885036},
                      "27": {"AtoB": 77.2185618682451, "AtoC": 63.3531736695171,
                             "AtoD": 16.2070057382602, "BtoC": 15.5526173038495,
                             "BtoD": 89.7131278799262, "CtoD": 76.026542897859}}

def call_command(command, pipe=None, echo=False):
    if echo:
        print command

    process = subprocess.Popen(shlex.split(command.encode("ascii")),
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    output = process.communicate(input=pipe)
    return output


def extract_pairing(structure):
    dict_open = dict()
    dict_close = dict()

    list_opener = []
    for i, char in enumerate(structure):
        if char == "(":
            list_opener.append(i)
        elif char == ")":
            # we add 1 to correct the fact that the 1st nt's index is 1 and not 0
            opener = list_opener.pop()+1
            closer = i+1
            dict_open[str(opener)] = closer
            dict_close[str(closer)] = opener

    return dict_open, dict_close


def compute_matures_indexes(hairpin_seq, structure, mature5p_seq, mature3p_seq):
    if mature5p_seq and mature3p_seq:
        start_5p = hairpin_seq.find(mature5p_seq) + 1 # correction to start counting from 1 instead of 0
        end_5p = start_5p + len(mature5p_seq) - 1 # correction to be inclusive
        start_3p = hairpin_seq.find(mature3p_seq) + 1 # correction to start counting from 1 instead of 0
        end_3p = start_3p + len(mature3p_seq) - 1 # correction to be inclusive
    else:
        dict_open, dict_close = extract_pairing(structure)
        if not mature5p_seq:
            start_3p = hairpin_seq.find(mature3p_seq) + 1 # correction to start counting from 1 instead of 0
            end_3p = start_3p + len(mature3p_seq) - 1 # correction to be inclusive
            range_3p = range(start_3p, end_3p+1)
            start_5p = None
            end_5p = None
            for pos in range_3p:
                if start_5p is None and str(pos) in dict_close:
                    start_5p = dict_close[str(pos)] + 2
                if str(pos) in dict_close:
                    end_5p = dict_close[str(pos)] + 2
        if not mature3p_seq:
            start_5p = hairpin_seq.find(mature5p_seq) + 1 # correction to start counting from 1 instead of 0
            end_5p = start_5p + len(mature5p_seq) - 1 # correction to be inclusive
            range_5p = range(start_5p, end_5p+1)
            start_3p = None
            end_3p = None
            for pos in range_5p:
                if start_3p is None and str(pos) in dict_open:
                    start_3p = min([dict_open[str(pos)] + 2, len(hairpin_seq)])
                if str(pos) in dict_open:
                    end_3p = min([dict_open[str(pos)] + 2, len(hairpin_seq)])

    return start_5p, end_5p, start_3p, end_3p


def get_vector_dict(mature5p_seq, mature3p_seq):
# def get_vector_dict(hairpin_length):
    if len(mature5p_seq) == 0:
        hairpin_length = len(mature3p_seq)
    elif len(mature3p_seq) == 0:
        hairpin_length = len(mature5p_seq)
    else:
        hairpin_length = min([len(mature3p_seq), len(mature5p_seq)])

    vector_dict = VECTOR_DICT_MAPPER[str(hairpin_length)]

    return vector_dict


def get_pdb_vector_dict(pdb_file, start_5p, end_5p, start_3p, end_3p, process_name):
    list_indexes = [start_5p, end_5p, start_3p, end_3p]
    list_a = []
    list_b = []
    list_c = []
    list_d = []
    list_lines = []
    with open(pdb_file, 'rb') as pdb_c:
        for line in pdb_c:
            if len(line) >= 50:
#                 splitted = line.split()
                # how to parse pdb with python according to http://cupnet.net/pdb_format/
                if line[12:16].strip() == "P" and int(line[22:26].strip()) in list_indexes:
                    index_in_int = int(line[22:26].strip())
                    if index_in_int == start_5p:
                        list_lines.append(line.strip())
                        list_a = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                    elif index_in_int == end_5p:
                        list_lines.append(line.strip())
                        list_b = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                    elif index_in_int == start_3p:
                        list_lines.append(line.strip())
                        list_c = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                    elif index_in_int == end_3p:
                        list_lines.append(line.strip())
                        list_d = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]

#             if list_a and list_b and list_c and list_d:
#                 break

    try:
        vector_dict = {"AtoB": math.sqrt(sum([(list_a[ind]-list_b[ind])**2 for ind in xrange(3)])),
                       "AtoC": math.sqrt(sum([(list_a[ind]-list_c[ind])**2 for ind in xrange(3)])),
                       "AtoD": math.sqrt(sum([(list_a[ind]-list_d[ind])**2 for ind in xrange(3)])),
                       "BtoC": math.sqrt(sum([(list_b[ind]-list_c[ind])**2 for ind in xrange(3)])),
                       "BtoD": math.sqrt(sum([(list_b[ind]-list_d[ind])**2 for ind in xrange(3)])),
                       "CtoD": math.sqrt(sum([(list_c[ind]-list_d[ind])**2 for ind in xrange(3)]))}
    except Exception as e:
        sys.stderr.write('{p} {l}\n'.format(p=process_name, l=list_indexes)) 
        sys.stderr.write(pdb_file + "\n")
        sys.stderr.write("\n".join(list_lines) + "\n")
        sys.stderr.write("\n".join([str(elem) for elem in list_a]) + "\n")
        sys.stderr.write("\n".join([str(elem) for elem in list_b]) + "\n")
        sys.stderr.write("\n".join([str(elem) for elem in list_c]) + "\n")
        sys.stderr.write("\n".join([str(elem) for elem in list_d]) + "\n")
        sys.stderr.write("###############################\n")
        vector_dict = dict()
        copy_to_fail_dir(pdb_file)

    return vector_dict


def compare_vectors(vector_dict, pdb_vector_dict, threshold):
    verdict = True
    for k, val in vector_dict.iteritems():
        diff = math.fabs(val - pdb_vector_dict[k])
        if diff > threshold:
            verdict = False
            break
    return verdict


def copy_to_out_dir(pdb_filepath, pdb, out_dir):
    shutil.copy(pdb_filepath,
                os.path.join(out_dir, pdb))
    call_command("gzip " + os.path.join(out_dir, pdb))


def copy_to_fail_dir(pdb_filepath):
    dest_path = os.path.join("/u/leongs/reproduction_projet_naim/rel20/3D/failed_processed", os.path.basename(pdb_filepath))
    shutil.copy(pdb_filepath,
                dest_path)
    call_command("gzip " + dest_path)


def process_pdb(params_dict):
    decoy_dir = params_dict["decoy_dir"]
    pdb = params_dict["pdb"]
    start_5p = params_dict["start_5p"]
    end_5p = params_dict["end_5p"]
    start_3p = params_dict["start_3p"]
    end_3p = params_dict["end_3p"]
    vector_dict = params_dict["vector_dict"]
    out_dir = params_dict["out_dir"]

    # incremental minimization if we're within the threshold
    pdb_filepath = os.path.join(decoy_dir, pdb)
    pdb_vector_dict = get_pdb_vector_dict(pdb_filepath, start_5p, end_5p, start_3p, end_3p, 'raw')
    if (not pdb_vector_dict is None) and compare_vectors(vector_dict, pdb_vector_dict, REGULAR_THRESHOLD):
        out, err = call_command("{relieve_script} {pdb_filepath}".format(relieve_script=relieve_script,
                                                                         pdb_filepath=pdb_filepath),
                                echo=False)
#         out, err = call_command("python {relieve} --pdb_file {pdb_filepath} --workdir {workdir}".format(relieve=relieve_script,
#                                                                                                         pdb_filepath=pdb_filepath,
#                                                                                                         workdir=decoy_dir),
#                                 echo=False)
        if err:
            sys.stderr.write("####################### relieve ###################\n")
            sys.stderr.write(pdb_filepath + "\n")
            sys.stderr.write(err + "\n")
            sys.stderr.write("#######################\n")
        pdb_vector_dict = get_pdb_vector_dict(pdb_filepath, start_5p, end_5p, start_3p, end_3p, 'relieve')
#         if (not pdb_vector_dict is None) and compare_vectors(vector_dict, pdb_vector_dict, RELIEVE_THRESHOLD):
        if True:
            out, err = call_command("{refine_script} {pdb_filepath}".format(refine_script=refine_script,
                                                                            pdb_filepath=pdb_filepath),
                                    echo=False)
#             out, err = call_command("python {refine_script} --pdb_file {pdb_filepath} --workdir {workdir}".format(refine_script=refine_script,
#                                                                                                                   pdb_filepath=pdb_filepath,
#                                                                                                                   workdir=decoy_dir),
#                                     echo=False)
            if err:
                sys.stderr.write("####################### refine ###################\n")
                sys.stderr.write(pdb_filepath + "\n")
                sys.stderr.write(err + "\n")
                sys.stderr.write("#######################\n")
            pdb_vector_dict = get_pdb_vector_dict(pdb_filepath, start_5p, end_5p, start_3p, end_3p, 'refine')
#             if (not pdb_vector_dict is None) and compare_vectors(vector_dict, pdb_vector_dict, THRESHOLD):
#                 copy_to_out_dir(pdb_filepath, pdb, out_dir)
#                 return True
#             elif (not pdb_vector_dict is None) and compare_vectors(vector_dict, pdb_vector_dict, REFINE_THRESHOLD):
            if True:
                out, err = call_command("{brushup_script} {pdb_filepath}".format(brushup_script=brushup_script,
                                                                                 pdb_filepath=pdb_filepath),
                                        echo=False)
#                 out, err = call_command("python {brushup} --pdb_file {pdb_filepath} --workdir {workdir}".format(brushup=brushup_script,
#                                                                                                                 pdb_filepath=pdb_filepath,
#                                                                                                                 workdir=decoy_dir),
#                                         echo=False)
                if err:
                    sys.stderr.write("####################### brushup ###################\n")
                    sys.stderr.write(pdb_filepath + "\n")
                    sys.stderr.write(err + "\n")
                    sys.stderr.write("#######################\n")
                pdb_vector_dict = get_pdb_vector_dict(pdb_filepath, start_5p, end_5p, start_3p, end_3p, 'brushup')
                if (not pdb_vector_dict is None) and compare_vectors(vector_dict, pdb_vector_dict, THRESHOLD):
                    copy_to_out_dir(pdb_filepath, pdb, out_dir)
                    return True
    return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--hairpin_seq',
                        action="store",
                        required=True,
                        dest="hairpin_seq",
                        help="The hairpin sequence")
    parser.add_argument('--mature5p_seq',
                        action="store",
                        dest="mature5p_seq",
                        help="The 5' mature sequence",
                        default="")
    parser.add_argument('--mature3p_seq',
                        action="store",
                        dest="mature3p_seq",
                        help="The 3' mature sequence",
                        default="")
    parser.add_argument('--structure',
                        action="store",
                        required=True,
                        dest="structure",
                        help="MC-Fold/MC-FlashFold output file")
    parser.add_argument('--decoy_dir',
                        action="store",
                        required=True,
                        dest='decoy_dir',
                        help="The path to the decoy")
    parser.add_argument('--refine_script',
                        action="store",
                        required=True,
                        dest='refine_script',
                        help="The path to the Refine script")
    parser.add_argument('--relieve_script',
                        action="store",
                        required=True,
                        dest='relieve_script',
                        help="The path to the Relieve script")
    parser.add_argument('--brushup_script',
                        action="store",
                        required=True,
                        dest='brushup_script',
                        help="The path to the Brush-Up script")
    parser.add_argument('--out_dir',
                        action="store",
                        required=True,
                        dest='out_dir',
                        help="The path where the valid PDB will be copied to")
    parser.add_argument('--threads',
                        action="store",
                        dest='threads',
                        type=int,
                        default=1,
                        help="Number of threads to run the process in")

    ns = parser.parse_args()

    hairpin_seq = ns.hairpin_seq
    mature5p_seq = ns.mature5p_seq
    mature3p_seq = ns.mature3p_seq
    structure = ns.structure
    decoy_dir = ns.decoy_dir
    refine_script = ns.refine_script
    relieve_script = ns.relieve_script
    brushup_script = ns.brushup_script
    out_dir = ns.out_dir
    threads = ns.threads

    start_5p, end_5p, start_3p, end_3p = compute_matures_indexes(hairpin_seq, structure, mature5p_seq, mature3p_seq)

    # not sure it's the best, we can probably use the shortest distance computed with start_5p, end_5p, start_3p, end_3p
    # instead of just chosing the mature seq
    vector_dict = get_vector_dict(mature5p_seq, mature3p_seq)
#     print start_5p, end_5p, start_3p, end_3p
#     print min([(end_5p-start_5p), (end_3p-start_3p)])
#     vector_dict = get_vector_dict(min([(end_5p-start_5p), (end_3p-start_3p)]))

    list_params_dict = [dict(decoy_dir=decoy_dir,
                             pdb=pdb,
                             start_5p=start_5p,
                             end_5p=end_5p,
                             start_3p=start_3p,
                             end_3p=end_3p,
                             vector_dict=vector_dict,
                             out_dir=out_dir) for pdb in sorted(os.listdir(decoy_dir)) if pdb.endswith(".pdb")]

    # now compute the decoy
    if list_params_dict:
        found = False
        for params_dict in list_params_dict:
            if found:
                break
            else:
                found = process_pdb(params_dict)
#         pool = multiprocessing.Pool(threads)
#         pool.map(process_pdb, list_params_dict, chunksize=1)