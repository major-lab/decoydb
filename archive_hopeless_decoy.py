import os
import cPickle
import tarfile

decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/completed_decoy"
processed_decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_23_jan/decoy"
processed_pbs_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed_23_jan/pbs"

list_failed_mir = []

for pbs_dir in os.listdir(processed_pbs_dir):
    err_file = os.path.join(processed_pbs_dir, pbs_dir, "err.txt")
    with open(err_file, 'rb') as err_f:
        for line in err_f:
            if line.startswith(processed_decoy_dir):
                decoy = line.strip().split("/")[-2]
                list_failed_mir.append(decoy)

list_found = []
for acc in os.listdir(processed_decoy_dir):
    out_dir = os.path.join(processed_decoy_dir, acc, "out")
    for pdb in os.listdir(out_dir):
        list_found.append(pdb.split("-")[0])


list_already_exists = []
for acc in sorted(os.listdir(decoy_dir)):
    list_sub_decoy = os.listdir(os.path.join(decoy_dir, acc))
    for sub_decoy in list_sub_decoy:
        list_content = os.listdir(os.path.join(decoy_dir, acc, sub_decoy))
        list_already_exists.extend(list_content)
        if not (sub_decoy in list_failed_mir or sub_decoy in list_found):
            print "rm -rf {a}".format(a=os.path.join(decoy_dir, acc, sub_decoy))

already_exists_pickle = "/u/leongs/reproduction_projet_naim/rel20/3D/already_exist.pk"

final_list = []
if os.path.exists(already_exists_pickle):
    dict_pickled_content = []
    with open(already_exists_pickle, 'rb') as pk:
        dict_pickled_content = cPickle.load(pk)
    pickle_content = []
    for k, v in dict_pickled_content.iteritems():
        pickle_content.extend(v)
    final_list = list(set(pickle_content+list_already_exists))
else:
    final_list = list(set(list_already_exists))

dict_already_exists = dict()
for elem in final_list:
    subdecoy = elem.split("-")[0]
    if not subdecoy in dict_already_exists:
        dict_already_exists[subdecoy] = []
    dict_already_exists[subdecoy].append(elem)

new_pickle = "/u/leongs/reproduction_projet_naim/rel20/3D/already_exist_new.pk"
with open(new_pickle, 'wb') as pk:
    cPickle.dump(dict_already_exists, pk, -1)