import os
import cPickle


best_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/best_struct_bis"
other_pv = "/u/leongs/reproduction_projet_naim/rel20/3D/for_pseudoviewer_bis"


for index, acc in enumerate(sorted(os.listdir(best_dir))):
    list_struct = []
    with open(os.path.join(best_dir, acc, acc + ".2d"), 'rb') as st_f:
        list_struct = [elem.split()[0] for elem in st_f.readlines() if elem.strip()]

    info_dict = dict()
    with open(os.path.join(best_dir, acc, acc + ".pk"), 'rb') as pk_f:
        info_dict = cPickle.load(pk_f)

    for i, struct in enumerate(list_struct):
        name = "{acc}_{i}".format(acc=acc, i=i)
        with open(os.path.join(other_pv, name + ".txt"), "wb") as pv_f:
            pv_f.write(">{name}\n{seq}\n{struct}".format(name=name,
                                                         seq=info_dict["seq"],
                                                         struct=struct))