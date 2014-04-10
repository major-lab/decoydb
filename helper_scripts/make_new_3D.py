import os
import cPickle


pbs_dir = "/u/leongs/reproduction_projet_naim/rel20/3D_bis/pbs"
decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D_bis/tmp"
best_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/best_struct_bis"
path_mcsymizer = "/u/leongs/workspace/mcweb/script/mcsymizer.py"
libaries_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/libraries"


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
                             "{exec_f} 12").format(exec_f=os.path.join(pbs_subdir, "launch.bash"),
                                                   stdout_filename=os.path.join(pbs_subdir, "out.txt"),
                                                   stderr_filename=os.path.join(pbs_subdir, "err.txt"),
                                                   name=str(index%modulo_val + 1)))


for index, acc in enumerate(sorted(os.listdir(best_dir))):
    list_struct = []
    with open(os.path.join(best_dir, acc, acc + ".2d"), 'rb') as st_f:
        list_struct = [elem.split()[0] for elem in st_f.readlines() if elem.strip()]

    info_dict = dict()
    with open(os.path.join(best_dir, acc, acc + ".pk"), 'rb') as pk_f:
        info_dict = cPickle.load(pk_f)

    library_pos = ""
    with open(os.path.join(libaries_dir, acc, "positions.txt"), 'rb') as p_f:
        library_pos = p_f.read().strip()

    list_command = []
    local_list_command = []
    for i, struct in enumerate(list_struct):
        workdir = os.path.join(decoy_dir, acc, "{acc}_{i}".format(acc=acc, i=i))
        local_list_command.append("mkdir " + workdir)
        local_list_command.append("cd "+ workdir)
        # make mcsym_script
        cmd = ('python {mcsymizer} --db_path {db_path} '
               '--sequence1 {seq} --structure1 "{struct}" '
               '--use_high_res_ncm '
               '--max_number 500 '
               '--name {acc}_{i} '
               '--external_library "{f},{pos}" '
               '--timeout 120 > {script_file}').format(mcsymizer=path_mcsymizer,
                                                      db_path="/u/leongs/Downloads/MC-Sym_Linux/MCSYM-DB",
                                                      seq=info_dict["seq"],
                                                      struct=struct,
                                                      acc=acc,
                                                      i=i+1,
                                                      f=os.path.join(libaries_dir, acc, acc+".pdb.gz"),
                                                      pos=library_pos,
                                                      script_file=os.path.join(workdir, "script.mcc"))

        local_list_command.append(cmd)

        cmd = "/soft/bioinfo/linux_RH5/mcsym-4.2.2/bin/mcsym script.mcc > /dev/null"
        local_list_command.append(cmd)
        local_list_command.append('echo "{workdir}"'.format(workdir=workdir))
    if local_list_command:
        list_command.append("mkdir {decoy_dir}".format(decoy_dir=os.path.join(decoy_dir, acc)))
        list_command.extend(local_list_command)
        write_pbs_script(list_command, acc, index)
