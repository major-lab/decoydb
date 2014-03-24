import os
import cPickle
import subprocess
import shlex
import shutil
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

workdir = "/u/leongs/reproduction_projet_naim/test_process_3d"

decoy_dir = os.path.join(workdir, "MI0000060")
out_dir = os.path.join(workdir, "out")

digested_list = []

with open("/u/leongs/reproduction_projet_naim/rel20/2D/digested_data.pk", 'rb') as pk_o:
    digested_list = cPickle.load(pk_o)

sequence = ""
mat_5p = ""
mat_3p = ""
for digested_data in digested_list:
    if digested_data["accession"] == "MI0000060":
        sequence = digested_data["sequence"]
        for mature in digested_data["matures"]:
            if mature["header"].endswith("5p"):
                mat_5p = mature["sequence"]
            else:
                mat_3p = mature["sequence"]
        break

with open("/u/leongs/reproduction_projet_naim/rel20/2D/filtered/MI0000060", 'rb') as mir_c:
    for i, line in enumerate(mir_c):

        name = "MI0000060_" + str(i)
        if not os.path.exists(os.path.join(decoy_dir, name)):
            continue
        struct = line.split()[0]

        print "start Analysis on MI0000060_" + str(i)

        start = datetime.datetime.now()
        out, err = call_command(('python /u/leongs/git/various-codes/DecoyDB/3_structure_verificator.py '
                                 '--hairpin_seq "{seq}" '
                                 '--mature5p_seq "{mature_5p}" '
                                 '--mature3p_seq "{mature_3p}" '
                                 '--structure "{struct}" '
                                 '--decoy_dir {subdecoy_dir} '
                                 '--refine_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/refine.bash '
                                 '--relieve_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/relieve.bash '
                                 '--brushup_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/brushup.bash '
                                 '--out_dir {out} --threads 15 ').format(seq=sequence,
                                                                         struct=struct,
                                                                         mature_5p=mat_5p,
                                                                         mature_3p=mat_3p,
                                                                         out=out_dir,
                                                                         subdecoy_dir=os.path.join(decoy_dir, name)),
                                echo=False)
        end = datetime.datetime.now()
        print end - start

        if err:
            print err