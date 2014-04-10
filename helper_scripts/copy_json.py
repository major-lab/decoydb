import os
import shutil

for acc in sorted(os.listdir("/u/leongs/git/mirdb/mirdecoy/data/2d")):
    if not os.path.exists(os.path.join("/u/leongs/git/mirdb/mirdecoy/data/additional_info",
                                       acc)):
        os.makedirs(os.path.join("/u/leongs/git/mirdb/mirdecoy/data/additional_info",
                                 acc))

    shutil.copy(os.path.join("/u/leongs/reproduction_projet_naim/rel20/3D/best_struct_bis",
                             acc,
                             acc+".json"),
                os.path.join("/u/leongs/git/mirdb/mirdecoy/data/additional_info",
                             acc,
                             acc+".json"))