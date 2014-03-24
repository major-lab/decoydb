import os

list_completed = ["completed_decoy_05_mar", "completed_decoy_10_mar"]

main_dir = "/u/leongs/reproduction_projet_naim/rel20/3D"
merge_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/merged_completed"

list_commands = []
for c_decoy in list_completed:
    completed_dir = os.path.join(main_dir, c_decoy)
    for acc in sorted(os.listdir(completed_dir)):
        mkdir_cmd = "mkdir "+os.path.join(merge_dir, acc)
        if not mkdir_cmd in list_commands:
            list_commands.append(mkdir_cmd)

        for sub_acc in sorted(os.listdir(os.path.join(completed_dir, acc))):
#             print sub_acc
            list_commands.append("mv {source} {dest}".format(source=os.path.join(completed_dir, acc, sub_acc),
                                                             dest=os.path.join(merge_dir, acc, sub_acc)))

print "\n".join(list_commands)