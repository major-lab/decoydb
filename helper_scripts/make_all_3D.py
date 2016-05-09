import os
import cPickle


flashfold_dir = "/u/mailhoto/workdir_decoy/flashfold_filtered"
workdir = "/u/mailhoto/workdir_decoy/3D"
pbs_dir = "/u/mailhoto/workdir_decoy/3D/pbs"
digested_data_pk = "/u/mailhoto/workdir_decoy/digested_data.pk"
path_mcsymizer = "/u/mailhoto/decoydb/mcsymizer.py"

# list_to_process = []
# with open("/u/leongs/reproduction_projet_naim/rel20/3D/failed_3D.txt", 'rb') as failed:
#     list_to_process = [elem.split()[0] for elem in failed]


def write_pbs_script(list_command, name, index):
    modulo_val = 3
    pbs_subdir = os.path.join(pbs_dir, str(index%modulo_val + 1))
    if (index%modulo_val + 1):
        if not os.path.exists(pbs_subdir):
            os.mkdir(pbs_subdir)
        pbs_folder = os.path.join(pbs_subdir, "pbs")
        if not os.path.exists(pbs_folder):
            os.mkdir(pbs_folder)
        if list_command:
            commands = ('#!/bin/sh\n'
                        '#PBS -N {name}\n'
                        '#PBS -l walltime=10000:00:00\n'
                        'export MCSYM_DB=/soft/bioinfo/share/mcsym/db/mcsymdb-4.2.1.bin.gz\n'
                        '{commands}').format(name=name,
                                             commands="\n".join(list_command))

            with open(os.path.join(pbs_folder, name + ".pbs"), 'wb') as pbs_f:
                pbs_f.write(commands)

            with open(os.path.join(pbs_subdir, "list.txt"), 'a') as list:
                list.write(os.path.join(pbs_folder, name + ".pbs") + '\n')

            with open(os.path.join(pbs_subdir, "execute.bash"), 'wb') as exec_bash:
                exec_bash.write("#!/bin/bash\n$* > /dev/null")

            with open(os.path.join(pbs_subdir, "launch.bash"), 'wb') as pbs_f:
                pbs_f.write("#!/bin/bash\ncat {list} | xargs --max-args=1 --max-procs=$* {execute_bash} \n".format(list=os.path.join(pbs_subdir, "list.txt"),
                                                                                                                   execute_bash=os.path.join(pbs_subdir, "execute.bash")))
    
            with open(os.path.join(pbs_subdir, "launch_pbs.pbs"), 'wb') as pbs_f:
                pbs_f.write(("#!/bin/bash\n"
                             "#PBS -N mcsym_{name}\n"
                             "#PBS -l nodes=1:ppn=12,walltime=10000:00:00\n"
                             "#PBS -o {stdout_filename}\n"
                             "#PBS -e {stderr_filename}\n"
                             "/u/mailhoto/workdir_decoy/3D/pbs/{name}/launch.bash 12").format(stdout_filename=os.path.join(pbs_subdir, "out.txt"),
                                                                                                            stderr_filename=os.path.join(pbs_subdir, "err.txt"),
                                                                                                            name=str(index%modulo_val + 1)))

list_digested_data = []
with open(digested_data_pk, 'rb') as dd:
        list_digested_data = cPickle.load(dd)

list_already_computed_decoy = []
# dec_dir = "/u/mailhoto/workdir_decoy/3D"
# for dec in ["completed_decoy_11_feb", "completed_decoy_17_feb", "completed_decoy_19_feb", "completed_decoy_21_feb", "completed_decoy_24_feb"]:
#     for acc in os.listdir(os.path.join(dec_dir, dec)):
#         list_already_computed_decoy.extend(os.listdir(os.path.join(dec_dir, dec, acc)))

dict_sum_content = dict()
# dec_dir = "/u/mailhoto/workdir_decoy/3D"
# for dec in ["processed_11_feb", "processed_17_feb", "processed_19_feb", "processed_21_feb"]:
#     list_acc = os.listdir(os.path.join(dec_dir, dec, "decoy"))
#     for acc in list_acc:
#         dict_sum_content[acc] = dict_sum_content.get(acc, 0) + len(os.listdir(os.path.join(dec_dir, dec, "decoy", acc, "out")))

for k, v in dict_sum_content.iteritems():
    if v > 0:
        list_already_computed_decoy.extend([k + "_" + str(i) for i in xrange(100)])


list_command = []
for index, hairpin_dict in enumerate(list_digested_data):
    list_command = []

    hairpin_name = hairpin_dict['name']
    hairpin_acc = hairpin_dict['accession']
    hairpin_seq = hairpin_dict['sequence']
    mcfold_output = os.path.join(flashfold_dir, hairpin_acc)

    hairpin_main_dir = os.path.join("/u/mailhoto/workdir_decoy/3D/decoys_100", hairpin_acc)

    local_list_command = []
    with open(mcfold_output, 'rb') as mcfold_o:
        for i, line in enumerate(mcfold_o):
            if line.strip():
                structure = line.split()[0]
            if i >= 100:
                break
            else:
                workdir = os.path.join(hairpin_main_dir, hairpin_acc + "_" + str(i))
                
                if "{name}_{i}".format(name=hairpin_acc, i=i) in list_already_computed_decoy:
                    continue
                
                local_list_command.append("mkdir "+workdir)
                local_list_command.append("cd "+ workdir)
                # make mcsym_script
                cmd = ('python {mcsymizer} --db_path {db_path} '
                       '--sequence1 {seq} --structure1 "{struct}" '
                       '--use_high_res_ncm '
                       '--max_number 100 '
                       '--name {acc}_{i} '
                       '--timeout 30 > {script_file}').format(mcsymizer=path_mcsymizer,
                                                              db_path="/soft/bioinfo/share/mcsym/db/mcsymdb-4.2.1.bin.gz",
                                                              seq=hairpin_seq,
                                                              struct=structure,
                                                              acc=hairpin_acc,
                                                              i=i,
                                                              script_file=os.path.join(workdir, "script.mcc"))

                local_list_command.append(cmd)

                cmd = "/soft/bioinfo/linux_RH5/mcsym-4.2.2/bin/mcsym script.mcc > /dev/null"
                local_list_command.append(cmd)
                local_list_command.append('echo "{workdir}"'.format(workdir=workdir))
    if local_list_command:
        list_command.append("mkdir {decoy_dir}".format(decoy_dir=hairpin_main_dir))
        list_command.extend(local_list_command)
        write_pbs_script(list_command, hairpin_acc, index)
        
        
        