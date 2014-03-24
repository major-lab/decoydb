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
Command-line utility to refine a pdb file
"""

import os
import argparse
import subprocess
import shlex
import shutil

def call_command(command, pipe=None, echo = False):
    if echo:
        print command
    
    process = subprocess.Popen(shlex.split(command.encode("ascii")),
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    
    output=process.communicate(input=pipe)
    return output

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--pdb_file',
                        action="store",
                        required=True,
                        dest="pdb_file",
                        help="The pdb_file")

    parser.add_argument('--workdir',
                        action="store",
                        required=True,
                        dest="workdir",
                        help="The workdir")

    ns = parser.parse_args()

    pdb_file = ns.pdb_file
    workdir = ns.workdir

    # prepare pdb_file to be minimized
    call_command('{addions} "{pdb_file}"'.format(addions=os.path.join('/u/major/bin/Tinker', 'addions.exe'), 
                                                 pdb_file=pdb_file))
    
    call_command('{get_sched} "{pdb_file}"'.format(get_sched=os.path.join('/u/major/bin/Tinker', 'get_sched.exe'), 
                                                   pdb_file=pdb_file))
    
    call_command('{prepchains} "{pdb_file}"'.format(prepchains=os.path.join('/u/major/bin/Tinker', 'prepchains.exe'), 
                                                    pdb_file=pdb_file))
    
    call_command('{pdbxyz} "{pdb_file}" ALL {amber99}'.format(pdbxyz=os.path.join('/u/major/bin/Tinker', 'pdbxyz'), 
                                                              pdb_file=pdb_file, 
                                                              amber99=os.path.join('/u/major/bin/Tinker', 'amber99')))

    xyz_file = pdb_file.replace('pdb', 'xyz') 
           
    minkey, err = call_command('{prepkey} "{xyz_file}" {steps}'.format(prepkey=os.path.join('/u/major/bin/Tinker', 'prepkey.exe'),
                                                                       xyz_file=xyz_file, 
                                                                       steps=200))

    with open(os.path.join(workdir, 'min.key'), 'w') as keyfile:
        keyfile.write(minkey.replace('/var/www/html/major_f/cgi-bin/Tinker/amber99', os.path.join('/u/major/bin/Tinker', 'amber99')))

    minimize_out, minimize_err = call_command('{minimize} "{xyz_file}" -k min.key {GRMS}'.format(minimize=os.path.join('/u/major/bin/Tinker', 'minimize'), 
                                                                                                xyz_file=xyz_file, 
                                                                                                GRMS=5.0))
    
    call_command('{xyzpdb} "{xyz_2}" {amber99}'.format(xyzpdb=os.path.join('/u/major/bin/Tinker', 'xyzpdb'), 
                                                       xyz_2=pdb_file.replace('pdb', 'xyz_2'), 
                                                       amber99=os.path.join('/u/major/bin/Tinker', 'amber99')))
    
    for f in [elem for elem in sorted(os.listdir(workdir)) if 'pdb_2' in elem]:
        os.rename(os.path.join(workdir, f), os.path.join(workdir, f.replace('pdb_2', 'pdb.min')))
            

    mcrmsd_out, mcsrmsd_err = call_command('{MCRMSD_exec_path} -n "{pdb_file}" "{pdb_file}.min"'.format(MCRMSD_exec_path='/soft/bioinfo/linux/mcrms-dev/bin/MC-RMSD', 
                                                                                                        pdb_file=pdb_file))

    # resulting minimized structure replaces original one
    for f in [elem for elem in sorted(os.listdir(workdir)) if 'pdb.min' in elem]:
        os.remove(os.path.join(workdir, f.replace("pdb.min", "pdb")))
        os.rename(os.path.join(workdir, f), os.path.join(workdir, f.replace('pdb.min', 'pdb')))

    call_command('{three2one} "{pdb_file}"'.format(three2one=os.path.join('/u/major/bin/Tinker', 'three2one.exe'), 
                                                   pdb_file=pdb_file))

    call_command('{set_sched} "{pdb_file}"'.format(set_sched=os.path.join('/u/major/bin/Tinker', 'set_sched.exe'), pdb_file=pdb_file))

    # delete the 'seq' and 'xyz'files
    to_be_rm = [elem for elem in sorted(os.listdir(workdir)) if ((elem.endswith('seq')) or (elem.endswith('xyz')) or elem.endswith('xyz_2'))]
    for file_rm in to_be_rm:
        try:
            os.remove(file_rm)
        except Exception as e:
            pass

    #delete the min.key file
    try:
        os.remove('min.key')
    except Exception as e:
        pass