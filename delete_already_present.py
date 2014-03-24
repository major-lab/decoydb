import os
import shutil
import cPickle

# a = "/u/leongs/reproduction_projet_naim/rel20/3D/decoy"
b = "/u/leongs/reproduction_projet_naim/rel20/3D/completed_decoy_28_jan"
# b = "/u/leongs/reproduction_projet_naim/rel20/3D/decoys_500"

already_present_pickle = "/u/leongs/reproduction_projet_naim/rel20/3D/already_exist.pk"

already_present_dict = dict()
with open(already_present_pickle, "rb") as pk:
    already_present_dict = cPickle.load(pk)

list_dir = os.listdir(b)

for acc in sorted(list_dir):
    bmir_dir = os.path.join(b, acc)

    bsubmir_list = os.listdir(bmir_dir)
    for sub in sorted(bsubmir_list):
        bsubmir_dir = os.path.join(bmir_dir, sub)

        list_pdb = [elem for elem in os.listdir(bsubmir_dir) if elem.endswith(".pdb.gz")]

        for pdb in list_pdb:
            bpdb = os.path.join(bsubmir_dir, pdb)
            if pdb in already_present_dict.get(sub, []):
                print bpdb
                os.remove(bpdb)