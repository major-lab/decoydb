#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# Copyright (C) 2013 Université de Montréal
#

"""
Command-line utility to generate an MC-Sym script
This utility is made to replace mcsymize by Marc Parisien

1.8 (2014-02-14): -add the option to use relative path for the db_path

1.7 (2013-12-04): -add a control to check whether the given db_path exists and can be accessed

1.6 (2013-11-04): -code refactoring
                  -fix an "index out of bounds" error

1.5 (2013-10-31): -fix the "moose structure" bug

1.4 (2013-10-30): -make the MCSYM-DB path same as the one given in the command-line. This way, no need for a symlink

1.3 (2013-10-29): -Another round of bugfixes

1.2 (2013-10-28): -Refactor the make_single_strand method

1.1 (2013-10-22): -Make it so that the distance restraints aren't put if both stems are not placed yet, they are placed otherwise
                  -Fixed a bug in which a link is placed to close a loop when the loop is actually an NCM
"""

import os
import argparse
import sys

__VERSION__ = 1.8

POSSIBLE_NCM_LIST = ["2_2", "2_3", "2_4", "2_5",
                     "3_2", "3_3", "3_4",
                     "4_2", "4_3", "4_4",
                     "6_2"]

CANONICAL_NCM_LIST = ["AAUU", "AUAU", "ACGU", "AGCU",
                      "UAUA", "UUAA", "UCGA", "UGCA",
                      "CAUG", "CUAG", "CCGG", "CGCG",
                      "GAUC", "GUAC", "GCGC", "GGCC"]

def is_canonical_ncm(ncm):
    """
    See if the NCM is made of canonical BPs. Is actually useful only when the
    user set 'use_high_res_ncm' as True
    """
    return ncm in CANONICAL_NCM_LIST

def get_danglings(starting_ind, struct1, struct2, list_unpaired_p):
    """
    Compute the Dangling part(s) of the structure(s)
    """
    def _check_structure_for_dangling(struct, dict_dangling):
        """
        Helper function to find the dangling in the struct
        """
        # check the 5p first
        if "(" in struct:
            curr_ind = 0
            while (curr_ind < len(struct) and \
                   (struct[curr_ind] != "(" or curr_ind+1 in list_unpaired_p)):
                dict_dangling["5p"].append(curr_ind+starting_ind)
                curr_ind += 1
        # now check the 3p
        if ")" in struct:
            curr_ind = len(struct)-1
            while (curr_ind >= 0 and (struct[curr_ind] != ")" or curr_ind+1 in list_unpaired_p)):
                dict_dangling["3p"].append(curr_ind+starting_ind)
                curr_ind -= 1
        return dict_dangling

    dict_dangling1 = {"5p":[], "3p":[]}
    dict_dangling2 = {"5p":[], "3p":[]}

    if struct1:
        dict_dangling1 = _check_structure_for_dangling(struct1, dict_dangling1)

    if struct2:
        dict_dangling2 = _check_structure_for_dangling(struct2, dict_dangling2)

    return dict_dangling1, dict_dangling2

def print_dangling(prefix, dang):
    """
    Print the section that corresponds to the Dangling ends in the MC-Sym script
    """
    list_dangling = []
    if dang["5p"]:
        sorted_dang = sorted(dang["5p"], reverse=True)
        list_dangling.append("(  {prefix}{f}  {remains}  )".format(prefix=prefix,
                                                                   f=sorted_dang[0]+1,
                                                                   remains = "  ".join(["{prefix}{n}".format(prefix=prefix,
                                                                                                             n=n) for n in sorted_dang])))
    if dang["3p"]:
        sorted_dang = sorted(dang["3p"])
        list_dangling.append("(  {prefix}{f}  {remains}  )".format(prefix=prefix,
                                                                   f=sorted_dang[0]-1,
                                                                   remains = "  ".join(["{prefix}{n}".format(prefix=prefix,
                                                                                                             n=n) for n in sorted_dang])))
    return list_dangling

def get_nucleotide_index(seq, pos):
    """
    Convert the nucleotide index to AX or BX so it fits the MC-Sym script formatting
    """
    if pos <= len(seq):
        nuc_index = "A{pos}".format(pos=pos)
    else:
        nuc_index = "B{pos}".format(pos=pos-len(seq))
    return nuc_index


def get_regular_stems(merged_seq, pairing_d, loop_d, external_library_dict_l,
                      stem_type="regular", external_library_r=[]):
    """
    'parse' through the sequence and try to separate into 'stems' and then extract the NCMs
    """
    list_stem = []
    list_ncm = []
    list_opening_keys = sorted([int(elem) for elem in pairing_d[stem_type].keys()])
    previous_key = list_opening_keys[0]
    list_left = []
    list_right = []
    already_inserted_external_library = []
    for op in list_opening_keys[1:]:
        # build the NCM
        left_start = min([previous_key, op])
        left_end = max([previous_key, op])
        right_start = min([pairing_d[stem_type][str(previous_key)], pairing_d[stem_type][str(op)]])
        right_end = max([pairing_d[stem_type][str(previous_key)], pairing_d[stem_type][str(op)]])

        list_left.extend(range(left_start, left_end+1))
        list_right.extend(range(right_start, right_end+1))

        # do not make NCMs for parts that are already in the external_library
        if left_start in external_library_r and left_end in external_library_r and \
           right_start in external_library_r and right_end in external_library_r:
            # determine which library we should load here
            for library_dict in external_library_dict_l:
                if left_start in library_dict["library_range"] and left_end in library_dict["library_range"] and \
                   right_start in library_dict["library_range"] and right_end in library_dict["library_range"] and \
                   not library_dict["library_file"] in already_inserted_external_library:
                    list_ncm.append(dict(sequence=merged_seq[library_dict["left_start"]-1:library_dict["left_end"]]+merged_seq[library_dict["right_start"]-1:library_dict["right_end"]],
                                         ncm="",
                                         left_start=library_dict["left_start"],
                                         left_end=library_dict["left_end"],
                                         right_start=library_dict["right_start"],
                                         right_end=library_dict["right_end"],
                                         is_loop=False,
                                         is_single_strand=False,
                                         library_path=library_dict["library_file"]))
                    already_inserted_external_library.append(library_dict["library_range"])
        else:
            left = left_end - left_start + 1
            right = right_end - right_start + 1

            ncm = "{l}_{r}".format(l=left, r=right)
            sequence = merged_seq[left_start-1:left_end]+merged_seq[right_start-1:right_end]

            # if not in ncm list, or not in MCSYM-DB, we probably have a stem that just reached its end
            if ncm not in POSSIBLE_NCM_LIST or \
               not (os.path.exists(os.path.join(db_path, ncm, sequence)) and \
                    [f for f in os.listdir(os.path.join(db_path, ncm, sequence)) if f.endswith("pdb.gz")]):
                if not list_ncm:
                    list_ncm.append(dict(sequence=merged_seq[left_start-1]+merged_seq[right_end-1],
                                         start=left_start,
                                         end=right_end,
                                         is_loop=False,
                                         is_single_strand=True,
                                         library_path=""))
                list_stem.append(list(list_ncm))
                list_ncm = []
                list_left = []
                list_right = []
            else:
                list_ncm.append(dict(sequence=sequence,
                                     ncm=ncm,
                                     left_start=left_start,
                                     left_end=left_end,
                                     right_start=right_start,
                                     right_end=right_end,
                                     is_loop=False,
                                     is_single_strand=False,
                                     library_path=""))
            if str(left_end) in loop_d:
                sequence = merged_seq[left_end-1:loop_d[str(left_end)]]
                ncm = len(sequence)
                start = left_end
                end = loop_d[str(left_end)]
                list_ncm.append(dict(sequence=sequence,
                                     ncm=ncm,
                                     start=start,
                                     end=end,
                                     is_loop=True,
                                     is_single_strand=False))
        previous_key = op
    if list_ncm:
        list_stem.append(list(list_ncm))
    return list_stem

def get_pairs(dotBracket):
    """
    INPUT dotBracket => OUTPUT dict containing lists of BP for regulars and pseudoknots
    """
    pair_dict = dict(regular=dict(), pseudoknot=dict())
    regular_opener = []
    pseudoknot_opener = []
    for i, element in enumerate(dotBracket):
        if element == "(":
            regular_opener.append(i)
        elif element == ")":
            pair_dict["regular"][str(regular_opener.pop(-1)+1)] = i+1
        elif element == "[":
            pseudoknot_opener.append(i)
        elif element == "]":
            pair_dict["pseudoknot"][str(pseudoknot_opener.pop(-1)+1)] = i+1

    return pair_dict

def find_loops(seq1, struct):
    """
    INPUT struct => OUTPUT list of dict with key=beginning of loop and value=end of loop
    """
    loop_d = dict()
    last_index = struct.find("(")
    i = last_index
    while len(struct)-1 > i >=0:
        i += 1
        if struct[i] == "(":
            last_index = i
        if struct[i] == ")":
            if not (last_index < len(seq1) <= i):
                loop_d[str(last_index+1)] = i+1
                last_index = struct.find("(", i+1)
                i = last_index

    return loop_d

def find_long_ss(seq1, struct, pair_dict,
                 list_placed_in_dang=[], list_unpaired_p=[]):
    """
    INPUT struct => OUTPUT list of tuple for the long single strands
    """
    list_long_ss = []
    begin = 0
    end = None
    # list_broken pairs is used to filter the pairs that can't be accepted (because they're not in an NCM)
    list_broken_pairs = list(list_placed_in_dang) + list_unpaired_p
    for op, cl in pair_dict["regular"].iteritems():
        if int(op) in list_broken_pairs or int(cl) in list_broken_pairs:
            list_broken_pairs.append(int(op))
            list_broken_pairs.append(int(cl))

    for i, char in enumerate(struct):
        if char in ['(', ')'] and not i+1 in list_broken_pairs:
            if begin != None and end != None and \
               not (begin < len(seq1) <= end+1):
                list_long_ss.append((begin+1, end+2))

            end = None
            begin = i
        else:
            end = i

    return list_long_ss


def validate_params(seq1, seq2, struct1, struct2, db_p):
    """
    Make sure the sequence and the structure are valid, i.e. same length, only AUGC, balanced pairs
    """
    merged_seq = seq1+seq2
    merged_struct = struct1+struct2
    if len(merged_seq) != len(merged_struct):
        sys.stderr.write("The sequence(s) and the structure(s) must be of the same length")
        sys.exit()

    for n in merged_seq:
        if n not in ["A", "U", "G", "C"]:
            sys.stderr.write("Invalid character found in sequence: {n}".format(n=n))
            sys.exit()

    if "()" in struct1 or "[]" in struct1 or \
       "()" in struct2 or "[]" in struct2:
        sys.stderr.write("Cannot model () or []")
        sys.exit()

    if (merged_struct.count("(") != merged_struct.count(")")) or \
       (merged_struct.count("[") != merged_struct.count("]")):
        sys.stderr.write("The structure is unbalanced")
        sys.exit()

    if not os.path.exists(db_p):
        sys.stderr.write("{db_path} doesn't exists".format(db_path=db_p))
        sys.exit()


def make_external_library_dict(external_l):
    """
    Parse the external library String to make it into a list of dict containing the following keys:
    """
    external_library_r = []
    external_library_dict_l = []
    if external_l:
        list_library_string = external_l.split(";")
        for library_string in list_library_string:
            library_params = library_string.split(",")
            if len(library_params) != 5:
                sys.stderr.write("Invalid number of comma-separated params in library_string")
                sys.exit()
            else:
                library_file = library_params[0].strip()
                try:
                    list_library_pos = sorted([int(elem.strip()) for elem in library_params[1:]])
                    if len(list_library_pos) != len(set(list_library_pos)):
                        sys.stderr.write("Invalid library positions : {library_string}".format(library_string=library_params[1:]))
                        sys.exit()
                    else:
                        left_start, left_end, right_start, right_end = list_library_pos
                except Exception:
                    sys.stderr.write("Invalid library positions : {library_string}".format(library_string=library_params[1:]))
                    sys.exit()
                if not os.path.exists(library_file):
                    sys.stderr.write("File not found: {library_file}".format(library_file=library_file))
                    sys.exit()
                external_library_dict = dict(library_file=os.path.abspath(library_file),
                                             left_start=left_start,
                                             left_end=left_end,
                                             right_start=right_start,
                                             right_end=right_end,
                                             library_range=range(left_start, left_end+1)+range(right_start, right_end+1))
                external_library_dict_l.append(external_library_dict)
                external_library_r.extend(external_library_dict["library_range"])
    return external_library_dict_l, external_library_r

def make_distance_restraint(seq, s, e, k=0):
    """
    Compute the distance restraints and format it to fit MC-Sym's specs
    """
    if k == 0:
        k = (e - s)
    distance_restraint_string = ""
    if k > 0:
        # the equation is (4.4 * k + 5.8)
        restraint_val = (4.4*k)+5.8
        distance_restraint_string = "distance(  {s_index}:C1'    {e_index}:C1'  0.0  {restraint}  )".format(s_index=get_nucleotide_index(seq=seq,
                                                                                                                                         pos=s),
                                                                                                            e_index=get_nucleotide_index(seq=seq,
                                                                                                                                         pos=e),
                                                                                                            restraint=restraint_val)
    return distance_restraint_string

# the is_loop variable isn't just for loops, it's to close 2 already placed nt (ex: a loop, or n stems that are already joined but not closed
def make_single_strand(start, end, sequence,
                       seq1, is_loop=False,
                       list_pairs_in_ncm=[],
                       list_placed_in_dangling=[],
                       library_diversity=None):
    """
    Build a single strand to either form a loop or to connect stems
    """
    list_library = []
    set_distance_restraints = set()

    rmsd_start_line = "  "
    if not library_diversity:
        ncm_rmsd = 0.5
        if library_diversity != None and float(library_diversity) == 0.0:
            rmsd_start_line = "//"
    else:
        ncm_rmsd = library_diversity

    placed_connecting = False
    for i in xrange(0,len(sequence)):
        j = i+1
        temp_seq = sequence[i:j+1]
        # no meaning linking a nucleotide to itself
        if len(temp_seq) > 1:
            temp_first = 1
            temp_second = 2
            temp_start = start + i
            temp_end = start + j
            if temp_start not in list_placed_in_dangling and temp_end not in list_placed_in_dangling:
                if not (placed_connecting and \
                        (temp_end in list_pairs_in_ncm or temp_start in list_pairs_in_ncm)) and \
                   (list_pairs_in_ncm.count(temp_start) < 2 and list_pairs_in_ncm.count(temp_end) < 2):

                    if use_relative_path:
                        pdb_path = os.path.join(os.path.basename(db_path), "ss2", temp_seq, "*.pdb.gz")
                    else:
                        pdb_path = os.path.join(db_path, "ss2", temp_seq, "*.pdb.gz")

                    list_library.append(('lnk = library(\n'
                                         '    pdb( "{pdb_path}" )\n'
                                         '    #{first}:#{second} <- {start}:{end}\n'
                                         '{rmsd_start_line}  rmsd( {rmsd} sidechain && !( pse || lp || hydrogen ) ) \n'
                                         ')').format(pdb_path=pdb_path,
                                                     first=temp_first,
                                                     second=temp_second,
                                                     start=get_nucleotide_index(seq=seq1,
                                                                                pos=temp_start),
                                                     end=get_nucleotide_index(seq=seq1,
                                                                              pos=temp_end),
                                                     rmsd_start_line=rmsd_start_line,
                                                     rmsd=ncm_rmsd))

                    if (temp_end in list_pairs_in_ncm or temp_start in list_pairs_in_ncm) and is_loop:
                        placed_connecting = True

                # make the distance restraints
                if temp_start != start:
                    dr_against_begin = make_distance_restraint(seq=seq1, s=start, e=temp_start)
                    if dr_against_begin:
                        set_distance_restraints.add(dr_against_begin)
                    dr_against_last = make_distance_restraint(seq=seq1, s=temp_start, e=end)
                    if dr_against_last:
                        set_distance_restraints.add(dr_against_last)
                if temp_end != end:
                    dr_against_begin = make_distance_restraint(seq=seq1, s=start, e=temp_end)
                    if dr_against_begin:
                        set_distance_restraints.add(dr_against_begin)
                    dr_against_last = make_distance_restraint(seq=seq1, s=temp_end, e=end)
                    if dr_against_last:
                        set_distance_restraints.add(dr_against_last)

    return list_library, sorted(list(set_distance_restraints))


def make_loop_library(start, end, sequence, db_path, ncm,
                      sequence1,
                      sequence2="",
                      list_pairs_in_ncm=[],
                      library_diversity=None):
    """
    Try to build the loop as an NCM. If not possible, build it manually by a combination of single strands
    """
    first = 1
    second = 1 + (end - start)
    list_library = []
    list_distance_restraints = []

    rmsd_start_line = "  "
    if not library_diversity:
        ncm_rmsd = 0.5
        if library_diversity != None and float(library_diversity) == 0.0:
            rmsd_start_line = "//"
    else:
        ncm_rmsd = library_diversity

    # see if we can find the loop directly in the NCMs
    if os.path.exists(os.path.join(db_path, str(len(sequence)), sequence)) and \
       [f for f in os.listdir(os.path.join(db_path, str(len(sequence)), sequence)) if f.endswith("pdb.gz")]:

        if use_relative_path:
            pdb_path = os.path.join(os.path.basename(db_path), str(len(sequence)), sequence, "*.pdb.gz")
        else:
            pdb_path = os.path.join(db_path, str(len(sequence)), sequence, "*.pdb.gz")

        list_library.append(('ncm = library(\n'
                             '    pdb( "{pdb_path}" )\n'
                             '    #{first}:#{second} <- {start}:{end}\n'
                             '{rmsd_start_line}  rmsd( {rmsd} sidechain && !( pse || lp || hydrogen ) ) \n'
                             ')').format(pdb_path=pdb_path,
                                         first=first,
                                         second=second,
                                         start=get_nucleotide_index(seq=sequence1,
                                                                    pos=ncm['start']),
                                         end=get_nucleotide_index(seq=sequence1,
                                                                  pos=ncm['end']),
                                         rmsd_start_line=rmsd_start_line,
                                         rmsd=ncm_rmsd))
    else:
        # since both ends are already placed, we don't put the last lnk
        # cut it to small single strands
        ss_list_library, list_distance_restraints = make_single_strand(start=start,
                                                                       end=end,
                                                                       sequence=sequence,
                                                                       seq1=sequence1,
                                                                       is_loop=True,
                                                                       list_pairs_in_ncm=list_pairs_in_ncm,
                                                                       list_placed_in_dangling=[],
                                                                       library_diversity=library_diversity)
        list_library.extend(ss_list_library)
    list_pairs_in_ncm.append(ncm['start'])
    list_pairs_in_ncm.append(ncm['end'])
    return list_library, list_distance_restraints


def make_ncm_library(ncm, use_high_res_ncm, sequence1, sequence2, library_diversity, library_path=""):
    """
    Format the NCM so it fits the format in MC-Sym script specs
    """
    first = 1
    second = first + (ncm['left_end'] - ncm['left_start'])
    third = second + 1
    fourth = third + (ncm['right_end'] - ncm['right_start'])
    rmsd_start_line = "  "
    if not library_diversity:
        ncm_rmsd = 0.5
        if is_canonical_ncm(ncm['sequence']) and use_high_res_ncm:
            ncm_rmsd = 0.1
        if library_diversity != None and float(library_diversity) == 0.0:
            rmsd_start_line = "//"
    else:
        ncm_rmsd = library_diversity

    if library_path:
        pdb_path = library_path
        type="lib"
    else:
        if use_relative_path:
            pdb_path = os.path.join(os.path.basename(db_path), ncm['ncm'], ncm['sequence'], '*.pdb.gz')
        else:
            pdb_path = os.path.join(db_path, ncm['ncm'], ncm['sequence'], '*.pdb.gz')
        type="ncm"

    library_st = ('{type} = library(\n'
                  '    pdb( "{pdb_path}" )\n'
                  '    #{first}:#{second}, #{third}:#{fourth} <- {l_s}:{l_e}, {r_s}:{r_e}\n'
                  '{rmsd_start_line}  rmsd( {rmsd} sidechain && !( pse || lp || hydrogen ) ) \n'
                  ')').format(type=type,
                              pdb_path=pdb_path,
                              first=first,
                              second=second,
                              third=third,
                              fourth=fourth,
                              l_s=get_nucleotide_index(seq=sequence1,
                                                       pos=ncm['left_start']),
                              l_e=get_nucleotide_index(seq=sequence1,
                                                       pos=ncm['left_end']),
                              r_s=get_nucleotide_index(seq=sequence1,
                                                       pos=ncm['right_start']),
                              r_e=get_nucleotide_index(seq=sequence1,
                                                       pos=ncm['right_end']),
                              rmsd_start_line=rmsd_start_line,
                              rmsd=ncm_rmsd)
    return library_st


def assemble_stems(list_regular_stem, library_diversity):
    """
    Get the stem list and tries to assemble each part of the stems
    """
    list_stem_library = []
    list_distance_restraints = []
    list_extremities = []
    list_pairs_in_ncm = []
    list_ss_extremities = []
    for i, stem in enumerate(list_regular_stem):
        list_left = []
        list_right = []
        list_loop_extremity = []
        list_library = []
        for j, ncm in enumerate(stem):
            if ncm['is_single_strand']:
                list_ss_extremities.extend((ncm['start'], ncm["end"]))
                list_ss_extremities = list(set(list_ss_extremities))
            else:
                if ncm['is_loop']:
                    # in a loop, there's no need to link the last 2 nts
                    loop_list_library, loop_list_distance_restraints = make_loop_library(sequence = ncm['sequence'],
                                                                                         start = ncm["start"],
                                                                                         end = ncm["end"],
                                                                                         ncm=ncm,
                                                                                         db_path=db_path,
                                                                                         sequence1=sequence1,
                                                                                         sequence2=sequence2,
                                                                                         list_pairs_in_ncm=list_pairs_in_ncm,
                                                                                         library_diversity=library_diversity)
                    list_library.extend(loop_list_library)
                    list_distance_restraints.extend(loop_list_distance_restraints)

                    list_loop_extremity.extend((ncm['start'], ncm['end']))
                    list_left.append(ncm['start'])
                    list_right.append(ncm["end"])
                else:
                    list_library.append(make_ncm_library(ncm, use_high_res_ncm, sequence1, sequence2,
                                                         library_diversity,
                                                         library_path=ncm['library_path']))

                    list_left.append(ncm['left_start'])
                    list_left.append(ncm['left_end'])
                    list_right.append(ncm["right_start"])
                    list_right.append(ncm["right_end"])

                    list_pairs_in_ncm.append(ncm['left_start'])
                    list_pairs_in_ncm.append(ncm['left_end'])
                    list_pairs_in_ncm.append(ncm['right_start'])
                    list_pairs_in_ncm.append(ncm['right_end'])

        if list_left and list_right:
            list_stem_extremities = [max(list_left),
                                     min(list_left),
                                     max(list_right),
                                     min(list_right)]

            if len(set(list_stem_extremities)) > len(set(list_loop_extremity)):
               list_stem_extremities = [elem for elem in list_stem_extremities if elem not in list_loop_extremity]

            list_stem_library.append(dict(list_index=list_left+list_right,
                                          list_library=list(list_library),
                                          list_extremities=list(set(list_stem_extremities))))

            list_extremities.extend(set(list_stem_extremities))
    return list_stem_library, list_distance_restraints, list_extremities, list_pairs_in_ncm, list_ss_extremities

def make_danglings(structure1, structure2, list_unpaired_pairs):
    """
    Make the Dangling part
    """
    list_relation = []
    dangling1 = dict()
    dangling2 = dict()
    list_placed_in_dangling = []
    dangling1, dangling2 = get_danglings(starting_ind=1,
                                         struct1=structure1,
                                         struct2=structure2,
                                         list_unpaired_p=list_unpaired_pairs)

    # make the relation part
    if dangling1["5p"] == dangling1["3p"][::-1]:
        dangling1["3p"] = []
        if dangling1["5p"]:
            dangling1["5p"].remove(max(dangling1["5p"]))
    for k, v in dangling1.iteritems():
        if v:
            v.sort()
            if k == "5p":
                one = v[0]
                two = v[-1]+1
            else:
                one = v[0]-1
                two = v[-1]
            list_relation.append('A{one}:  A{two} {{ file( "helixA_RNA" ) stack }} 1'.format(one=one,
                                                                                             two=two))
            list_placed_in_dangling.extend(v)
    if sequence2:
        if dangling2["5p"] == dangling2["3p"][::-1]:
            dangling2["3p"] = []
            if dangling2["5p"]:
                dangling2["5p"].remove(max(dangling2["5p"]))
        for k, v in dangling2.iteritems():
            if v:
                v.sort()
                if k == "5p":
                    one = v[0]
                    two = v[-1]+1
                else:
                    one = v[0]-1
                    two = v[-1]
                list_relation.append('B{one}:  B{two} {{ file( "helixA_RNA" ) stack }} 1'.format(one=one,
                                                                                                 two=two))
                list_placed_in_dangling.extend(v)

    return list_relation, dangling1, dangling2, list_placed_in_dangling


def merge_libraries(list_stem_library, list_extremities, list_ss_extremities,
                    sequence1, sequence2, merged_sequence, long_ss_list,
                    db_path, list_pairs_in_ncm, list_placed_in_dangling,
                    library_diversity):
    i = 1
    list_library = ["//----- Fragment {i} -----".format(i=i), ]
    list_library.extend(list_stem_library[0]["list_library"])
    if len(list_stem_library) > 1:
        sorted_set_extremities = sorted(list(set(list_extremities)))

        list_merged_stem = ["0", ]

        list_placed_unavailable_extremities = [elem for elem in [min(sorted_set_extremities), max(sorted_set_extremities)] if elem not in list_ss_extremities]
        list_placed_available_extremities = [extremity for extremity in sorted_set_extremities \
                                             if (extremity in list_stem_library[0]['list_index'] \
                                                              and extremity not in list_placed_unavailable_extremities)]

        index = 1
        while True:
            best_prev = None
            best_curr = None
            dict_distance = dict()
            for i_prev, prev in enumerate(sorted(list_placed_available_extremities)):
                remaining_extremities = sorted([ext for ext in sorted_set_extremities if \
                                                ext not in list_placed_unavailable_extremities and \
                                                ext != prev])

                for curr in remaining_extremities:
                    distance = abs(curr - prev)
                    if str(distance) not in dict_distance:
                        dict_distance[str(distance)] = []
                    dict_distance[str(distance)].append((prev, curr))

            back_stem = None
            forward_stem = None
            for dist in xrange(1, len(merged_sequence) + 1):
                if str(dist) in dict_distance:
                    for best_pairs in dict_distance[str(dist)]:
                        best_prev = min(best_pairs)
                        best_curr = max(best_pairs)
                        if int(dist) > 1 and not [(a,b) for (a,b) in long_ss_list if a==best_prev and b==best_curr]:
                            continue
                        # find the corresponding STEMs
                        for stem_index, stem in enumerate(list_stem_library):
                            if best_prev in stem['list_index']:
                                back_stem = stem_index
                            if best_curr in stem['list_index']:
                                forward_stem = stem_index
                            if back_stem and forward_stem and back_stem != forward_stem:
                                break
                        if back_stem != forward_stem:
                            break
                        else:
                            back_stem = None
                            forward_stem = None
                if back_stem != forward_stem:
                    break
                else:
                    back_stem = None
                    forward_stem = None

            if back_stem!=None and forward_stem!=None and back_stem != forward_stem:
                prev = best_prev
                curr = best_curr
                is_loop = False

                if str(back_stem) in list_merged_stem and str(forward_stem) in list_merged_stem:
                    is_loop = True

                ss_list_library = []
                if not (is_loop and len(merged_sequence[prev-1:curr]) <= 2):
                    # make and put the link
                    ss_list_library, ss_list_distance_restraints = make_single_strand(start=prev,
                                                                                      end=curr,
                                                                                      sequence=merged_sequence[prev-1:curr],
                                                                                      is_loop=is_loop,
                                                                                      seq1=sequence1,
                                                                                      list_pairs_in_ncm=list_pairs_in_ncm,
                                                                                      list_placed_in_dangling=list_placed_in_dangling,
                                                                                      library_diversity=library_diversity)

                    if str(back_stem) in list_merged_stem and str(forward_stem) in list_merged_stem:
                        list_distance_restraints.extend(ss_list_distance_restraints)

                # make sure the back stem is already in list_library before linking
                if str(back_stem) not in list_merged_stem:
                    i += 1
                    list_library.append("//----- Fragment {i} -----".format(i=i))
                    list_library.extend(ss_list_library[::-1])

                    # if the stem's list library is actually a loop formed by only lnk, decide whether to reverse it or not
                    list_lib_has_ncm = [lib for lib in list_stem_library[back_stem]["list_library"] if lib.startswith("ncm")]
                    if not list_lib_has_ncm:
                        list_library.extend(list_stem_library[back_stem]["list_library"][::-1])
                    else:
                        list_library.extend(list_stem_library[back_stem]["list_library"])

                    list_merged_stem.append(str(back_stem))
                    list_placed_available_extremities.extend([stem_extremity for stem_extremity \
                                                              in list_stem_library[back_stem]["list_extremities"] \
                                                              if stem_extremity not in list_placed_unavailable_extremities+list_placed_available_extremities])
                else:
                    list_library.extend(ss_list_library)

                # put the stem on if not merged yet
                if str(forward_stem) not in list_merged_stem:
                    i += 1
                    list_library.append("//----- Fragment {i} -----".format(i=i+1))
                    list_library.extend(list_stem_library[forward_stem]["list_library"])
                    list_merged_stem.append(str(forward_stem))

                    list_placed_available_extremities.extend([stem_extremity for stem_extremity \
                                                              in list_stem_library[forward_stem]["list_extremities"] \
                                                              if stem_extremity not in list_placed_unavailable_extremities+list_placed_available_extremities])

                list_placed_available_extremities = list(set(list_placed_available_extremities))

                if best_prev in list_ss_extremities:
                    list_ss_extremities.remove(best_prev)
                else:
                    if best_prev in list_placed_available_extremities:
                        list_placed_available_extremities.remove(best_prev)
                    list_placed_unavailable_extremities.append(best_prev)

                if best_curr in list_ss_extremities:
                    list_ss_extremities.remove(best_curr)
                else:
                    if best_curr in list_placed_available_extremities:
                        list_placed_available_extremities.remove(best_curr)
                    list_placed_unavailable_extremities.append(best_curr)

            else:
                break
    return list_library


def print_script(sequence1, structure1,
                 sequence2, structure2,
                 list_relation,
                 new_list_library,
                 list_merge,
                 list_distance_restraints,
                 clash_threshold,
                 exploration_method,
                 construction_method,
                 bond_threshold,
                 max_number,
                 timeout,
                 model_diversity,
                 name,
                 unzipped,
                 print_header):
    if print_header:
        print ('// MC-Sym 4.2 script generated by mcsymizer.py {version}\n'
               '// (c) Stephen Leong Koan & Francois Major, University of Montreal\n'
               '// Please cite: Parisien M, Major F. Nature. (2008) 452:51-55.\n'
               '//\n'
               '// web site FAQ: www.major.iric.ca/MC-Sym/faq.html\n'
               '//\n'
               '\n'
               '\n').format(version=__VERSION__)

    print '// ==================== Sequence ====================\n'

    if sequence1:
        print ('sequence( r A1  {seq} )\n'
               '//              {struct}\n').format(seq=sequence1,
                                                    struct=structure1)
    if sequence2:
        print ('sequence( r B1  {seq} )\n'
               '//              {struct}\n').format(seq=sequence2,
                                                    struct=structure2)

    if list_relation:
        print "\n".join(["// ==================== Relations ====================",
                         "relation",
                         "(",
                         "\n".join(list_relation),
                         ")"])


    # print library section
    print "// ==================== Library ====================\n"
    print "\n".join(new_list_library)

    # print backtrack/merge section
    print ('// ===================== Backtrack =====================\n'
           'structure = backtrack\n'
           '(')
    for elem in list_merge:
        print "    {elem}".format(elem=elem)
    print ')\n'

    # if applicable, make the distance restraints part
    if list_distance_restraints:
        print "// =================== Distance Restraints ==================="
        print "\n".join(list_distance_restraints)

    if unzipped:
        zipped = ""
    else:
        zipped = "zipped"


    # print remaining
    print ('// =================== Backtrack Restraints ===================\n'
           'clash\n'
           '(\n'
           '    structure\n'
           '    {ct} !( pse || lp || hydrogen )\n'
           ')\n'
           'backtrack_rst\n'
           '(\n'
           '    structure\n'
           '    width_limit  = 25%,\n'
           '    height_limit = 33%,\n'
           '    method       = {em}\n'
           ')\n'
           '// =================== Ribose Restraints ===================\n'
           'ribose_rst\n'
           '(\n'
           '    structure\n'
           '    method    = {cm},\n'
           '    pucker    = C3p_endo,\n'
           '    threshold = {bt}\n'
           ')\n'
           '// =================== Exploration Initialization =========\n'
           'explore\n'
           '(\n'
           '    structure\n'
           '    option(\n'
           '    model_limit = {mn},\n'
           '    time_limit  = {t}m,\n'
           '    seed        = 3210 )\n'
           '    rmsd( {md} sidechain && !( pse || lp || hydrogen ) )\n'
           '    pdb( "{name}" {zipped} )\n'
           ')').format(ct=clash_threshold,
                       em=exploration_method,
                       cm=construction_method,
                       bt=bond_threshold,
                       mn=max_number,
                       t=timeout,
                       md=model_diversity,
                       name=name,
                       zipped=zipped)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Generate a MC-Sym script",
                                     epilog="Copyright(c) Stephen Leong Koan & Francois Major, University of Montreal")

    group = parser.add_mutually_exclusive_group()

    parser.add_argument('-D', '--db_path', action="store", dest="db_path",
                        required=True,
                        help='MCSYM-DB directory')

    parser.add_argument('-s1', '--sequence1', action="store", dest="sequence1",
                        required=True,
                        help='RNA sequence (mandatory)')

    parser.add_argument('-S1', '--structure1', action="store", dest="structure1",
                        required=True,
                        help='RNA structure in dot-brackets notation (mandatory)')

    parser.add_argument('-s2', '--sequence2', action="store", dest="sequence2",
                        default="",
                        help='Annealed RNA sequence')

    parser.add_argument('-S2', '--structure2', action="store", dest="structure2",
                        default="",
                        help='Structure of the annealed RNA in dot-brackets notation')

    parser.add_argument('-n', '--name', action="store", dest="name",
                        default="structure",
                        help='Name')

    parser.add_argument('-mr', '--merge_rmsd', action="store", dest="merge_rmsd",
                        type=float, default=1.5,
                        help='Merge RMSD threshold')

    parser.add_argument('-md', '--model_diversity', action="store", dest="model_diversity",
                        type=float, default=3.0,
                        help='Model diversity threshold in Angstrom')

    group.add_argument('-ld', '--library_diversity', action="store", dest="library_diversity",
                        type=float, default=None,
                        help='Library diversity threshold in Angstrom')

    group.add_argument('--use_high_res_ncm', action="store_true",
                        dest="use_high_res_ncm",
                        help='Use high-resolution NCMs for canonical stacks of base pairs')

    parser.add_argument('-ct', '--clash_threshold', action="store", dest="clash_threshold",
                        type=float, default=1.5,
                        help='Clash threshold. Prevents non-bonded atoms to be too close to one another')

    parser.add_argument('-cm', '--construction_method', action="store", dest="construction_method",
                        choices=["ccm", "estimate"], default="ccm",
                        help='ccm=Cyclic Coordinate Minimization, estimate=Interpolation Estimation')

    parser.add_argument('--bond_threshold', '-bt', action="store", dest="bond_threshold",
                        type=float, default=2.0,
                        help='Bond threshold for covalent bonds in the backbone')

    parser.add_argument('-em', '--exploration_method', action="store", dest="exploration_method",
                        choices=["probabilistic", "exhaustive"], default="probabilistic",
                        help='probabilistic=allow back-jumps and random domain assignments, exhaustive=classic back-track algorithm')

    parser.add_argument('-mn', '--max_number', action="store", dest="max_number",
                        type=int, default=1000,
                        help='Maximum number of models generated')

    parser.add_argument('-t', '--timeout', action="store", dest="timeout",
                        type=int, default=30,
                        help='Time limit in minutes')

    parser.add_argument('-nd', '--no_dangling', action="store_false", dest="no_dangling",
                        help='Do not include dangling ends in model')

    parser.add_argument('-u', '--unzipped', action="store_true", dest="unzipped",
                        help='Do not zip the model(s)')

    parser.add_argument('-el', '--external_library', action="store", dest="external_library",
                        default='',
                        help=('Use the pdb.gz fragment as library, at which position. '
                              '(e.g. fragment1.pdb.gz,5,7,21,23;fragment2.pdb.gz,9,12,17,19). '
                              'N.B: we recommend protecting the input string using quotes (").'))

    parser.add_argument('-ur', '--use_relative_path', action="store_true",
                        dest="use_relative_path",
                        help='Use relative_paths in the NCM paths')

    parser.add_argument('-nh', '--no_header', action="store_false",
                        dest="no_header",
                        help='Do not print the header in the script')

    ns = parser.parse_args()

    sequence1 = ns.sequence1
    sequence2 = ns.sequence2
    structure1 = ns.structure1
    structure2 = ns.structure2
    db_path = os.path.abspath(ns.db_path)
    name = ns.name
    merge_rmsd = ns.merge_rmsd
    model_diversity = ns.model_diversity
    clash_threshold = ns.clash_threshold
    construction_method = ns.construction_method
    bond_threshold = ns.bond_threshold
    exploration_method = ns.exploration_method
    max_number = ns.max_number
    timeout = ns.timeout
    use_high_res_ncm = ns.use_high_res_ncm
    dangling = ns.no_dangling
    library_diversity = ns.library_diversity
    external_library = ns.external_library
    unzipped = ns.unzipped
    use_relative_path = ns.use_relative_path
    print_header = ns.no_header

    # convert the sequences
    sequence1 = sequence1.upper().replace("T", "U")
    sequence2 = sequence2.upper().replace("T", "U")

    merged_sequence = sequence1+sequence2
    merged_structure = structure1+structure2
    # validate to make sure everything is ok before making computations
    validate_params(seq1=sequence1,
                    seq2=sequence2,
                    struct1=structure1,
                    struct2=structure2,
                    db_p=db_path)

    # make the library_dict
    external_library_dict_list, external_library_range = make_external_library_dict(external_l=external_library)

    # get each pairs
    pairing_dict = get_pairs(dotBracket=merged_structure)

    # get the position of each loop
    loop_dict = find_loops(seq1=sequence1,
                           struct=merged_structure)

    # build the NCMS within stems
    list_regular_stem = get_regular_stems(merged_seq=merged_sequence,
                                          pairing_d=pairing_dict,
                                          loop_d=loop_dict,
                                          stem_type="regular",
                                          external_library_dict_l=external_library_dict_list,
                                          external_library_r=external_library_range)

    # build the ncms and lnk for stems
    list_stem_library, list_distance_restraints, \
    list_extremities, list_pairs_in_ncm, list_ss_extremities = assemble_stems(list_regular_stem, library_diversity)

    #make the distance for the pairs that are not in a x_y NCM
    list_unpaired_pairs = []
    for opening, closing in pairing_dict["regular"].iteritems():
        if int(opening) not in list_pairs_in_ncm:
            list_unpaired_pairs.append(int(opening))
            list_unpaired_pairs.append(int(closing))
            restraint = 3.0
            list_distance_restraints.append("distance(  {s_index}:C1'    {e_index}:C1'  0.0  {restraint}  )".format(s_index=get_nucleotide_index(seq=sequence1,
                                                                                                                                                 pos=int(opening)),
                                                                                                                    e_index=get_nucleotide_index(seq=sequence1,
                                                                                                                                                 pos=int(closing)),
                                                                                                                    restraint=restraint))

    # do the dangling
    if dangling:
        list_relation, dangling1, dangling2, \
        list_placed_in_dangling = make_danglings(structure1, structure2, list_unpaired_pairs)
    else:
        list_relation = []
        dangling1 = dict()
        dangling2 = dict()
        list_placed_in_dangling = []

    # get the long ss
    long_ss_list = find_long_ss(seq1=sequence1, struct=merged_structure,
                                list_placed_in_dang=list_placed_in_dangling,
                                pair_dict=pairing_dict,
                                list_unpaired_p=list_unpaired_pairs)

    # make the merge part and add remaining links
    # start by putting the first stem
    list_library = merge_libraries(list_stem_library, list_extremities, list_ss_extremities,
                                   sequence1, sequence2, merged_sequence, long_ss_list,
                                   db_path, list_pairs_in_ncm, list_placed_in_dangling,
                                   library_diversity)

    list_passed = []
    list_merge = []
    new_list_library = []
    curr_ncm_index = 1
    curr_lnk_index = 1
    curr_lib_index = 1
    for lib in list_library:
        curr_line = lib
        if curr_line in list_passed:
            continue
        else:
            list_passed.append(curr_line)
        if curr_line.startswith("ncm"):
            ncm_index = "ncm_{index:02d}".format(index=curr_ncm_index)
            curr_line = curr_line.replace("ncm =", ncm_index + " =")
            if not list_merge:
                text_merge = ncm_index
            else:
                text_merge = 'merge( {ncm_index}  {merge_rmsd} )'.format(ncm_index=ncm_index,
                                                                         merge_rmsd=merge_rmsd)
            list_merge.append(text_merge)
            curr_ncm_index += 1

        elif curr_line.startswith("lnk"):
            lnk_index = "lnk_{index:02d}".format(index=curr_lnk_index)
            curr_line = curr_line.replace("lnk =", lnk_index + " =")
            if not list_merge:
                text_merge = lnk_index
            else:
                text_merge = 'merge( {lnk_index}  {merge_rmsd} )'.format(lnk_index=lnk_index,
                                                                         merge_rmsd=merge_rmsd)
            list_merge.append(text_merge)
            curr_lnk_index += 1
        elif curr_line.startswith("lib"):
            lib_index = "lib_{index:02d}".format(index=curr_lib_index)
            curr_line = curr_line.replace("lib =", lib_index + " =")
            if not list_merge:
                text_merge = lib_index
            else:
                text_merge = 'merge( {lib_index}  {merge_rmsd} )'.format(lib_index=lib_index,
                                                                         merge_rmsd=merge_rmsd)
            list_merge.append(text_merge)
            curr_lib_index += 1

        new_list_library.append(curr_line)


    # add the danglings
    if dangling and (dangling1 or dangling2):
        list_merge.append("//------ dangling ends --------")
        list_merge.extend(print_dangling(prefix="A", dang=dangling1))
        list_merge.extend(print_dangling(prefix="B", dang=dangling2))

    # compute the distance restraints for pseudoknots
    for pseudo_open, pseudo_close in pairing_dict["pseudoknot"].iteritems():
        list_distance_restraints.append(make_distance_restraint(sequence1, int(pseudo_open), int(pseudo_close), k=1))

    print_script(sequence1, structure1,
                 sequence2, structure2,
                 list_relation,
                 new_list_library,
                 list_merge,
                 list_distance_restraints,
                 clash_threshold,
                 exploration_method,
                 construction_method,
                 bond_threshold,
                 max_number,
                 timeout,
                 model_diversity,
                 name,
                 unzipped,
                 print_header)