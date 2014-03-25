import os
import cPickle
import tarfile

decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/completed_decoy_10_mar"
processed_decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_10_mar/decoy"
processed_pbs_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_10_mar/pbs"


dict_failed_mir = dict()

for pbs_dir in os.listdir(processed_pbs_dir):
    list_failed_mir = []
    err_file = os.path.join(processed_pbs_dir, pbs_dir, "err.txt")
    with open(err_file, 'rb') as err_f:
        for line in err_f:
            if line.startswith("/tmp"):
                decoy = line.strip().split("/")[-1]
                list_failed_mir.append(decoy)

    for elem in list_failed_mir:
        splitted = elem.split("-")
        subdecoy = splitted[0]
        subdecoy_index = int(splitted[1].replace(".pdb", "").replace(".gz", ""))
        if not subdecoy in dict_failed_mir:
            dict_failed_mir[subdecoy] = []
        dict_failed_mir[subdecoy].append(subdecoy_index)


for acc in os.listdir(processed_decoy_dir):
    out_dir = os.path.join(processed_decoy_dir, acc, "out")
    for pdb in os.listdir(out_dir):
        splitted = pdb.split("-")
        subdecoy = splitted[0]
        subdecoy_index = int(splitted[1].replace(".pdb.gz", ""))
        if not subdecoy in dict_failed_mir:
            dict_failed_mir[subdecoy] = []
        dict_failed_mir[subdecoy].append(subdecoy_index)

# dict_max_index = dict()
# for subdecoy, list_index in dict_failed_mir.iteritems():
#     if len(list_index) > 1:
#         print subdecoy
#     starting_ind = 10
#     while starting_ind < max(list_index):
#         starting_ind += 10
#     dict_max_index[subdecoy] = starting_ind

for acc in sorted(os.listdir(decoy_dir)):
    for sub_acc in sorted(os.listdir(os.path.join(decoy_dir, acc))):
        for pdb in os.listdir(os.path.join(decoy_dir, acc, sub_acc)):
            if pdb.endswith("pdb.gz"):
                splitted = pdb.split("-")
                subdecoy_index = int(splitted[1].replace(".pdb.gz", ""))
                if not subdecoy_index in dict_failed_mir.get(sub_acc, []):
                    print "rm -rf {p}".format(p=os.path.join(decoy_dir, acc, sub_acc, pdb))