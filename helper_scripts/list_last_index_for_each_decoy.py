import os
import cPickle
import tarfile

decoy_dir_list = ["/u/leongs/reproduction_projet_naim/rel20/3D/decoy",
                  "/u/leongs/reproduction_projet_naim/rel20/3D/completed_decoy_23_Jan",
                  "/u/leongs/reproduction_projet_naim/rel20/3D/completed_decoy_28_jan",
                  "/u/leongs/reproduction_projet_naim/rel20/3D/decoys_500"]


dict_failed_mir = dict()

for decoy_dir in decoy_dir_list:
    for acc in sorted(os.listdir(decoy_dir)):
        for sub_acc in sorted(os.listdir(os.path.join(decoy_dir, acc))):
            if not sub_acc in dict_failed_mir:
                dict_failed_mir[sub_acc] = [50,]

            for pdb in os.listdir(os.path.join(decoy_dir, acc, sub_acc)):
                if pdb.endswith("pdb.gz"):
                    splitted = pdb.split("-")
                    subdecoy_index = int(splitted[1].replace(".pdb.gz", ""))
                    dict_failed_mir[sub_acc].append(subdecoy_index)

for sub_decoy in sorted(dict_failed_mir.keys()):
    list_index = dict_failed_mir[sub_decoy]
    print sub_decoy, max(list_index)