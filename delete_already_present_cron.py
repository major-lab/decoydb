#!/usr/bin/env python

#PBS -N cron_delete                
#PBS -o /u/leongs/reproduction_projet_naim/rel20/3D/cron_job_for_deleting_already_present/out.txt        
#PBS -e /u/leongs/reproduction_projet_naim/rel20/3D/cron_job_for_deleting_already_present/err.txt
#PBS -l walltime=10000:00:00

import os
import shutil
import cPickle
import subprocess
import shlex
import time
import datetime


def call_command(command, pipe=None, echo=False):
    if echo:
        print command

    process = subprocess.Popen(shlex.split(command.encode("ascii")),
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    output = process.communicate(input=pipe)
    return output

b = "/u/leongs/reproduction_projet_naim/rel20/3D/decoys_500"

already_present_pickle = "/u/leongs/reproduction_projet_naim/rel20/3D/cron_job_for_deleting_already_present/already_exist.pk"

already_present_dict = dict()
with open(already_present_pickle, "rb") as pk:
    already_present_dict = cPickle.load(pk)

last_status = False
while True:
    print "working", datetime.datetime.now()
    list_dir = os.listdir(b)

    number_removed = 0
    for acc in sorted(list_dir):
        bmir_dir = os.path.join(b, acc)

        bsubmir_list = os.listdir(bmir_dir)
        for sub in sorted(bsubmir_list):
            bsubmir_dir = os.path.join(bmir_dir, sub)

            list_pdb = [elem for elem in os.listdir(bsubmir_dir) if elem.endswith(".pdb.gz")]

            for pdb in list_pdb:
                bpdb = os.path.join(bsubmir_dir, pdb)
                if pdb in already_present_dict.get(sub, []):
                    number_removed += 1
                    os.remove(bpdb)

    print number_removed

    qsub_exec = os.path.join("/opt/torque-2.5.9/bin", "qstat")
    stdout, stderr = call_command(qsub_exec)

    monitoring_job_id = ("599357.binsrv3", '599358.binsrv3', '599359.binsrv3', '599361.binsrv3', '599361.binsrv3')
    all_finished_running = len([splitted for splitted in \
                                (elem.split() for elem in stdout.splitlines()) \
                                if splitted[2] == "leongs" and splitted[0] in monitoring_job_id]) == 0

    if not all_finished_running:
        print "sleeping", datetime.datetime.now()
        time.sleep(300)
    elif last_status and all_still_running:
        break
    last_status = all_finished_running