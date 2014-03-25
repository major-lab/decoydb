import os
import cPickle
import shutil


flashfold_dir = "/u/leongs/reproduction_projet_naim/rel20/2D/filtered"
tmp_dir = "/tmp"
workdir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_10_mar"
pbs_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_10_mar/pbs"
digested_data_pk = "/u/leongs/reproduction_projet_naim/rel20/2D/digested_data.pk"
raw_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/completed_decoy_10_mar"
job_name = "10_mar_{index}"

def write_pbs_script(list_command, list_hairpins, index):
    curr_index = index%3 + 1
    pbs_subdir = os.path.join(pbs_dir, str(curr_index))
    if not os.path.exists(pbs_subdir):
        os.mkdir(pbs_subdir)
    pbs_folder = os.path.join(pbs_subdir, "pbs")
    if not os.path.exists(pbs_folder):
        os.mkdir(pbs_folder)
    if list_command:
        name = list_hairpins[0]
        commands = ('#!/bin/sh\n'
                    '#PBS -N {name}\n'
                    '#PBS -l walltime=10000:00:00\n'
                    '{commands}').format(name=name,
                                         commands="\n".join(list_command))

        with open(os.path.join(pbs_folder, name + ".pbs"), 'wb') as pbs_f:
            pbs_f.write(commands)

        with open(os.path.join(pbs_subdir, "list.txt"), 'a') as list:
            list.write(os.path.join(pbs_folder, name + ".pbs") + '\n')

        with open(os.path.join(pbs_subdir, "list_acc"), 'a') as l_acc:
            l_acc.write("{d}\n".format(d=name))

        with open(os.path.join(pbs_subdir, "execute.bash"), 'wb') as exec_bash:
                exec_bash.write("#!/bin/bash\n$*")

        with open(os.path.join(pbs_subdir, "launch.bash"), 'wb') as pbs_f:
            pbs_f.write("#!/bin/bash\ncat {list} | xargs --max-args=1 --max-procs=$* {execute_bash} \n".format(list=os.path.join(pbs_subdir, "list.txt"),
                                                                                                               execute_bash=os.path.join(pbs_subdir, "execute.bash")))

        with open(os.path.join(pbs_subdir, "launch_pbs.pbs"), 'wb') as pbs_f:
            pbs_f.write(("#!/bin/bash\n"
                         "#PBS -N {job_name}\n"
                         "#PBS -l nodes=1:ppn=12,walltime=10000:00:00\n"
                         "#PBS -o {stdout_filename}\n"
                         "#PBS -e {stderr_filename}\n"
                         "rm -rf /tmp/MI*\n"
                         "{pbs_dir}/{name}/launch.bash 12").format(stdout_filename=os.path.join(pbs_subdir, "out.txt"),
                                                                  stderr_filename=os.path.join(pbs_subdir, "err.txt"),
                                                                  pbs_dir=pbs_dir,
                                                                  job_name=job_name.format(index=str(curr_index)),
                                                                  name=str(curr_index)))

list_digested_data = []
with open(digested_data_pk, 'rb') as dd:
    list_digested_data = cPickle.load(dd)

list_command = []
list_hairpins = []
for index, hairpin_dict in enumerate(list_digested_data):
    list_command = []
    list_hairpins = []

    hairpin_name = hairpin_dict['name']
    hairpin_acc = hairpin_dict['accession']
    hairpin_seq = hairpin_dict['sequence']

    mat_5p = ""
    mat_3p = ""
    for mature in hairpin_dict["matures"]:
        if mature["header"].endswith("5p"):
            mat_5p = mature["sequence"]
        else:
            mat_3p = mature["sequence"]

    mcfold_output = os.path.join(flashfold_dir, hairpin_acc)

    list_hairpins.append(hairpin_acc)

    #hairpin_main_dir = os.path.join(workdir, 'decoy', hairpin_acc)
    hairpin_main_dir = os.path.join(tmp_dir, hairpin_acc)
    hairpin_out_dir = os.path.join(workdir, 'decoy', hairpin_acc, "out")


    local_list_command = []
    with open(mcfold_output, 'rb') as mcfold_o:
        for i, line in enumerate(mcfold_o):
#             if i >= 50:
#                 continue
            if line.strip():
                structure = line.split()[0]
                raw_decoy_dir = os.path.join(raw_dir, hairpin_acc, hairpin_acc + "_" + str(i))
                if os.path.exists(raw_decoy_dir):
                    print raw_decoy_dir
                    list_pdb = [elem for elem in os.listdir(raw_decoy_dir) if elem.endswith(".pdb.gz")]
                    if len(list_pdb) == 0:
                        continue
                    destination_dir = os.path.join(hairpin_main_dir, hairpin_acc + "_" + str(i))
                    local_list_command.append("cp -r {source} {dest}".format(source=raw_decoy_dir,
                                                                             dest=destination_dir))

                    local_list_command.append("cd {dest}".format(dest=destination_dir))
                    local_list_command.append("gunzip " + os.path.join(destination_dir, "*.pdb.gz"))

                    local_list_command.append(('python /u/leongs/git/various-codes/DecoyDB/3_structure_verificator.py '
                                               '--hairpin_seq "{seq}" '
                                               '--mature5p_seq "{mature_5p}" '
                                               '--mature3p_seq "{mature_3p}" '
                                               '--structure "{struct}" '
                                               '--decoy_dir {subdecoy_dir} '
                                               '--refine_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/refine.bash '
                                               '--relieve_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/relieve.bash '
                                               '--brushup_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/brushup.bash '
                                               '--out_dir {out} --threads 1 ').format(seq=hairpin_seq,
                                                                                      struct=structure,
                                                                                      mature_5p=mat_5p,
                                                                                      mature_3p=mat_3p,
                                                                                      out=hairpin_out_dir,
                                                                                      subdecoy_dir=destination_dir))

#                     local_list_command.append(('python /u/leongs/git/various-codes/DecoyDB/3_structure_verificator.py '
#                                                '--hairpin_seq "{seq}" '
#                                                '--mature5p_seq "{mature_5p}" '
#                                                '--mature3p_seq "{mature_3p}" '
#                                                '--structure "{struct}" '
#                                                '--decoy_dir {subdecoy_dir} '
#                                                '--refine_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/refine.py '
#                                                '--relieve_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/relieve.py '
#                                                '--brushup_script /u/leongs/git/various-codes/DecoyDB/minimization_scripts/brushup.py '
#                                                '--out_dir {out} --threads 1 ').format(seq=hairpin_seq,
#                                                                                       struct=structure,
#                                                                                       mature_5p=mat_5p,
#                                                                                       mature_3p=mat_3p,
#                                                                                       out=hairpin_out_dir,
#                                                                                       subdecoy_dir=destination_dir))

                    local_list_command.append("rm -rf " + destination_dir)

    if local_list_command:
        list_command.append("mkdir {decoy_dir}".format(decoy_dir=hairpin_main_dir))
        list_command.append("mkdir {decoy_dir}".format(decoy_dir=os.path.join(workdir, 'decoy', hairpin_acc)))
        list_command.append("mkdir {hairpin_out_dir}".format(hairpin_out_dir=hairpin_out_dir))
        list_command.extend(local_list_command)

        write_pbs_script(list_command, list_hairpins, index)