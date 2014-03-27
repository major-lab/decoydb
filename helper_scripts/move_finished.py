import os

not_finished_file = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_10_mar/pbs/2/list_acc"

not_finished_list = []

with open(not_finished_file, 'rb') as nff:
    not_finished_list = [elem.strip() for elem in nff if elem.strip()]

old_decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_10_mar/decoy"
new_decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/finished_processing_26_mar"

for acc in sorted(os.listdir(old_decoy_dir)):
    if acc not in not_finished_list:
        out_dir = os.path.join(old_decoy_dir, acc, 'out')
        if os.listdir(out_dir):
            print "mkdir " + os.path.join(new_decoy_dir, acc)
            for pdb in sorted(os.listdir(out_dir)):
                print "cp {pdb} {dest}".format(pdb=os.path.join(out_dir, pdb),
                                               dest=os.path.join(new_decoy_dir, acc, pdb))